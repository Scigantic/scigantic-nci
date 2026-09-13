# scigantic-nci: build brief for the module authors

You are writing one module of a small public Python package (MIT-0, `Scigantic/scigantic-nci` on GitHub, published to PyPI). The package gives one-import access to two public Scigantic mirrors of NCI data. Read this whole file before writing code.

## What already exists (do not rewrite)

- `pyproject.toml`: name, deps (pandas, pyarrow, numpy, boto3; `pydicom` only in the `dicom` extra), `scigantic-nci` console script pointing at `scigantic_nci.cli:main`.
- `src/scigantic_nci/_store.py`: `resolve(bucket, key, root=None) -> Store`, `Store.read_parquet(rel, columns=None, filters=None)`, `read_bytes`, `read_json`, `exists`, `listdir`, `walk_files`, `parquet_schema`, and `list_keys(bucket)`. It picks local dir / `$SCIGANTIC_NCI_ROOT` / the notebook mount / anonymous S3. Use it for every read. Never call boto3 directly for the mirror buckets. `NciError` and `NotMirroredError` live there.
- `src/scigantic_nci/__init__.py` imports `gdc` and `idc` modules. Your module must import cleanly with only the base dependencies (pydicom imported lazily inside functions).
- `.venv/` has Python + pandas 2.x + pyarrow + boto3 + pydicom + pytest + mypy. Use `./.venv/bin/python` and `./.venv/bin/pytest`.
- `docs/table_schemas.md`: the exact parquet schemas of one GDC project and one IDC collection.

## Where the data is

- GDC: `s3://scigantic-gdc-open/<PROJECT>/` (57 projects: TCGA-*, TARGET-*, CPTAC-2/3, HCMI-CMDC, MMRF-COMMPASS, BEATAML1.0-COHORT, CGCI-*, ...). Not all projects are uploaded yet; a complete local copy of the ones built so far is at `/private/tmp/claude-501/-Users-aaronkanzer-projects-scigantic/c826c167-1893-4b78-a1d1-5610f06f1463/scratchpad/gdc_mirror/<PROJECT>/`. Point `SCIGANTIC_NCI_ROOT` at `.../scratchpad/gdc_mirror` for tests (the resolver accepts `<root>/<key>`).
- IDC: `s3://scigantic-idc-open/<collection_id>/` (176 collections, most already uploaded; local copy at `.../scratchpad/idc_mirror/<collection_id>/`). Same env var trick with `idc_mirror`.
- Both buckets are public-read; `Store` reads them anonymously. The raw NCI buckets (`tcga-2-open`, `idc-open-data`, ...) are public too; newer GDC objects are SSE-KMS and refuse `botocore.UNSIGNED`, so a *signed* default boto3 client is needed for GDC raw files (any principal works), while IDC raw objects accept unsigned reads.

### GDC per-project layout (every table documented in the project's README.md)

```
MANIFEST.parquet                 every open file: file_id, file_name, data_category, data_type, data_format,
                                 experimental_strategy, file_size, md5sum, project_id, case_id, case_submitter_id,
                                 sample_submitter_id, sample_type, workflow_type, platform, aliquot_submitter_id,
                                 s3_bucket, s3_key, gdc_download_url
clinical/cases.parquet           one row per case; demographic_* columns flattened; demographic_vital_status,
                                 demographic_days_to_death, demographic_age_at_index (years), demographic_sex_at_birth
clinical/diagnoses.parquet       MORE rows than cases (metastasis/recurrence rows); diagnosis_is_primary_disease flags
                                 the primary; days_to_last_follow_up; stage column varies (ajcc_pathologic_stage,
                                 ensat_pathologic_stage, figo_stage, ...); age_at_diagnosis in DAYS
clinical/treatments|exposures|follow_ups|samples|aliquots.parquet
rna_seq/genes.parquet            gene_id (Ensembl with version), gene_name, gene_type; 60,660 rows GENCODE v36;
                                 45 ids end in _PAR_Y and duplicate X-chromosome names; 6 more names are duplicated
rna_seq/samples.parquet          column -> aliquot_submitter_id, file_id, case_submitter_id, sample_submitter_id,
                                 sample_type, tissue_type, tumor_descriptor
rna_seq/counts_unstranded.parquet   gene_id + one int32 column per aliquot (raw STAR counts)
rna_seq/tpm_unstranded.parquet      gene_id + one float32 column per aliquot
rna_seq/star_qc.parquet
somatic_mutations/masked_somatic_mutations.parquet   all MAFs stacked; file_id, aliquot_submitter_id + 140 MAF
                                 columns (Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_depth, ...)
copy_number/gene_level_<workflow>.parquet   gene_id, gene_name, chromosome, start, end + one float32 column per
                                 tumor aliquot; workflows: ascat2, ascat3, absolute_liftover, ascatngs (separate matrices)
copy_number/gene_level_<workflow>_samples.parquet, copy_number/segments_<wf>.parquet, masked_segments_dnacopy,
                                 allele_specific_segments_<wf>
mirna/read_count.parquet, mirna/rpm.parquet (miRNA_ID + columns), mirna/samples.parquet, mirna/isoforms.parquet
protein/rppa_protein_expression.parquet (AGID, peptide_target + columns), protein/antibodies.parquet, protein/samples.parquet
methylation/betas_<platform>.parquet (probe_id + float32 columns, 50k-probe row groups), *_samples.parquet   [optional, may be absent]
README.md, BUILD_REPORT.json ({project, gdc_data_release, tables: {name: {rows, columns, bytes}}, ...})
```

Column ids of the wide matrices are aliquot barcodes; join to cases through the sibling `*samples.parquet` (`column` -> `case_submitter_id`, `sample_type`). Some projects lack some table groups (FM-AD has only clinical). Tables can be big: read wide matrices with `columns=` (gene_id plus the aliquots you want) or by row-group; never load methylation whole.

### IDC per-collection layout

```
series.parquet       one row per DICOM series: collection_id, PatientID, StudyInstanceUID, SeriesInstanceUID, Modality,
                     BodyPartExamined, StudyDate, SeriesDescription, instanceCount, series_size_MB, license_short_name,
                     source_DOI, aws_bucket, crdc_series_uuid, series_aws_url, PatientSex, PatientAge, Manufacturer, ...
studies.parquet      one row per study: PatientID, StudyDate, modalities, n_series, instances, size_MB, body_parts
patients.parquet     one row per patient: PatientSex, PatientAge, n_studies, n_series, modalities, size_MB
sm_series.parquet    slide microscopy attributes (staining_usingSubstance_CodeMeaning is LIST-valued, ObjectiveLensPower, ...)  [when SM present]
analysis_results.parquet  [when present]
clinical/<table>.parquet + clinical/dictionary.parquet (short_table_name, column, column_label, values)  [59 collections]
sample/<Modality>_<SeriesInstanceUID>/*.dcm   whole small radiology series, or thumbnail + 2 lowest levels of one slide
README.md, BUILD_REPORT.json ({collection_id, collection_name, samples: [{modality, folder, kind: series|sm_levels, ...}],
                     modalities, licenses, source_DOIs, clinical_tables, n_patients, n_series, n_studies, total_GB})
```

Raw IDC objects: `s3://<aws_bucket>/<crdc_series_uuid>/<instance_uuid>.dcm` plus a folder-marker object ending in `/` (skip it). Anonymous reads work. Slide (SM) series are DICOM WSI pyramids: one instance per level, tiled frames; assemble a level row-major with `ncol = ceil(TotalPixelMatrixColumns / Columns)`, canvas shape `(rows, cols) + frames.shape[3:]` (RGB slides have a channel axis, fluorescence slides do not); a level with NumberOfFrames == 1 is a single tile.

## Rules for the code

- Small, typed, documented functions that return pandas DataFrames; no classes beyond what is needed. `from __future__ import annotations`; mypy strict must pass (`./.venv/bin/mypy src`).
- Every public function takes the project/collection id first and an optional `root: str | None = None` last; pass it to `resolve`.
- Validate ids: unknown project/collection raises `NotMirroredError` with the list of valid ids in the message (cache `list_keys` per process).
- Do not print. Do not download whole matrices when a column subset is asked for. Do not use duckdb, polars, s3fs, or anything outside the declared dependencies.
- Tests in `tests/test_<module>.py` run against the LOCAL mirror via `SCIGANTIC_NCI_ROOT` (set it in a `conftest.py` fixture from the scratchpad paths above, skip if absent) AND one small live-S3 test per module marked with `@pytest.mark.network` (IDC: `acrin_6698` is uploaded; GDC: skip the live test if `s3://scigantic-gdc-open/TCGA-CHOL/BUILD_REPORT.json` is absent, using `Store.exists`). Assert real measured numbers (e.g. TCGA-BLCA has 412 cases, 431 RNA aliquots, 19 Solid Tissue Normal; tcga_lihc has 3,457 series, 377 patients).
- Measure before you claim: any number in a docstring or README section must come from the data you read.
- No emoji, no em-dashes, plain prose in docstrings.

## Deliverables

Write your module, its tests, and a `docs/<module>.md` with a usage walkthrough whose every snippet you actually ran (paste the real output). Finish with `./.venv/bin/pytest -q` and `./.venv/bin/mypy src` clean, and report: functions implemented, test counts, anything in the data you found surprising.
