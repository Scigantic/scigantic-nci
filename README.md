<h1 align="center">scigantic-nci</h1>

<p align="center">
    <a href="https://github.com/Scigantic/scigantic-nci/actions/workflows/ci.yml">
        <img alt="CI" src="https://github.com/Scigantic/scigantic-nci/actions/workflows/ci.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/scigantic-nci/">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/scigantic-nci" /></a>
    <a href="https://pypi.org/project/scigantic-nci/">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/scigantic-nci" /></a>
    <a href="https://github.com/Scigantic/scigantic-nci/blob/main/LICENSE">
        <img alt="License" src="https://img.shields.io/github/license/Scigantic/scigantic-nci" /></a>
</p>

Analysis-ready NCI cancer data in one import: per-project [Genomic Data Commons](https://gdc.cancer.gov/) tables (clinical, RNA-Seq, mutations, copy number, miRNA, RPPA, methylation) and per-collection [Imaging Data Commons](https://imaging.datacommons.cancer.gov/) indexes with sample DICOM, read as pandas DataFrames from public S3 or a mounted copy. No account, no download step, no BigQuery.

```python
import scigantic_nci as nci

nci.gdc.projects()                                                  # one row per mirrored GDC project
clin = nci.gdc.clinical("TCGA-BLCA")                                # 412 cases: os_time_days, os_event, stage
tpm = nci.gdc.expression("TCGA-BLCA", genes=["TP53", "ESR1"], sample_type="Primary Tumor")   # (2, 412)
maf = nci.gdc.mutations("TCGA-BLCA", genes=["TP53"], non_silent=True)                        # MAF rows
freq = nci.gdc.mutation_frequency("TCGA-BLCA")                      # TP53 0.488, TTN 0.425, KMT2D 0.263, ...

nci.idc.collections()                                               # 176 collections, 1,032,911 series
s = nci.idc.series("tcga_lihc", modality="CT")                      # 777 rows, one per DICOM series
vol = nci.idc.read_sample("tcga_lihc", "CT")                        # Volume (36, 512, 512) in HU, from the mirror
nci.idc.pull_series(s.iloc[0], "/tmp/ct")                           # the full series from idc-open-data, anonymous
```

## Installation

```console
$ pip install scigantic-nci
$ pip install "scigantic-nci[dicom]"     # pydicom and the pylibjpeg decoders, for read_sample() and pulled series
```

Python 3.10 or newer. Dependencies are pandas, pyarrow, numpy and boto3.

## Why this exists

Both NCI data commons are public on S3, and both are stored in a way that no one can use without an index they do not ship.

The GDC's open buckets (`tcga-2-open`, `gdc-cptac-phs001287-2-open`, ...) hold one folder per file, named by the file's UUID: every key is `<file_id>/<file_name>`, with nothing above it for program, project, case or data type. Across the 73 projects with open files that is 390,418 folders; the 33 TCGA projects alone are 284,498. Which of those folders belong to TCGA-BLCA is only knowable through the GDC API, and every folder holds one per-sample file: TCGA-BRCA has 1,231 `Gene Expression Quantification` TSVs, one per aliquot, each listing the same 60,660 genes, so its expression matrix is 1,231 downloads followed by a column-wise stack. The newer GDC objects are also encrypted with SSE-KMS and refuse anonymous reads.

The IDC's buckets (`idc-open-data`, `idc-open-data-cr`, `idc-open-data-two`) hold 1,032,911 DICOM series across 176 collections, 96.9 TB, at `s3://<bucket>/<crdc_series_uuid>/<instance_uuid>.dcm`. The key carries no collection, patient or modality; that mapping is the IDC index, served from BigQuery or the `idc-index` package's downloadable tables.

This package reads two Scigantic mirrors that hold those indexes and the analysis-ready tables:

- `s3://scigantic-gdc-open/<PROJECT>/`: every open-access file of a project stacked into parquet tables, one prefix per project. The build so far covers 36 projects (the 33 TCGA projects and TARGET-ALL-P1, -P2, -P3): 13,230 cases, 1,178 tables, 4.97 GB of parquet, built from 168,844 source files (178.4 GB) read from the GDC buckets, at Data Release 46.0 (August 10, 2026). More projects are being added; `gdc.projects()` lists what is live.
- `s3://scigantic-idc-open/<collection_id>/`: for each of the 176 collections a `series`, `studies` and `patients` table (`sm_series` for the 73 with slides, `analysis_results` for 55, clinical tables for 59), plus sample DICOM that can be read without touching the raw buckets: 286 whole small radiology series and the thumbnail and two lowest pyramid levels of 47 slides, 20,897 files, 5.23 GB. Built from idc-index 0.12.5, IDC data v24.

Each prefix carries a `README.md` describing its tables and a `BUILD_REPORT.json`. For GDC the report records the source file count, files and bytes read, md5 checks against the GDC manifest, duplicates and any re-alignment (`realigned_files`); the tables are a re-arrangement of the GDC files with values unchanged from the source. For IDC the report records the series, patient and study counts, modalities, licenses and source DOIs of the collection and of every sample. Both mirrors were built on 2026-09-13.

What the functions do about size: column subsets and row filters are pushed down to parquet, so `expression("TCGA-BRCA", genes=[5 symbols])` reads 5 x 1,231 values in 0.10 s instead of the 296 MB full matrix (0.37 s locally), and `methylation(..., probes=[...])` reads only the 50,000-probe row groups whose statistics can contain those ids. The numbers are in [docs/gdc.md](docs/gdc.md).

Note on names: `gdc.samples()` is the GDC biospecimen table (one row per sample of a case); `idc.samples()` lists the sample DICOM series shipped with a collection. The two modules use the same word for different things because their sources do.

## Data provenance and citation

**GDC.** Only the open-access tier is mirrored; it carries no access restrictions. Cite the program that generated the data (TCGA, TARGET, ...) and the GDC: Grossman et al., *Toward a Shared Vision for Cancer Genomic Data*, NEJM 2016, [doi:10.1056/NEJMp1607591](https://doi.org/10.1056/NEJMp1607591). TCGA and TARGET publication guidelines: <https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing>.

**IDC.** Every series carries its own license in `series.parquet` (`license_short_name`) and its source DOI (`source_DOI`); across the mirror 865,935 series are CC BY 4.0, 132,303 CC BY 3.0, 28,783 CC BY-NC 4.0, 5,851 CC BY-NC 3.0 and 39 under the National Library of Medicine terms. Collections can mix licenses, so check the row, not the collection. The CC BY-NC series (28,783 + 5,851) may not be used for commercial purposes; `license_short_name` is authoritative and the only place the license is recorded, so filter on it before using a series in commercial work, and on the same column of `idc.samples(c)` before using a mirrored sample series. Cite the source DOIs (`idc.report(c)["source_DOIs"]`, 237 distinct across the mirror) and the IDC: Fedorov et al., *NCI Imaging Data Commons*, Radiographics 2021, [doi:10.1148/rg.2021210020](https://doi.org/10.1148/rg.2021210020); Fedorov et al., *National Cancer Institute Imaging Data Commons: Toward Transparency, Reproducibility, and Scalability in Imaging Artificial Intelligence*, Nature Methods 2023, [doi:10.1038/s41592-023-01893-2](https://doi.org/10.1038/s41592-023-01893-2).

This package's code is MIT-0. That does not extend to the data.

## Public API

Every per-project or per-collection function takes the id first and `root: str | None = None` last. Walkthroughs with real output: [docs/gdc.md](docs/gdc.md), [docs/idc.md](docs/idc.md).

### `scigantic_nci.gdc`

| function | returns |
|---|---|
| `projects()` | one row per mirrored project (release, tables, cases, files and bytes read), cached |
| `report(p)`, `tables(p)`, `readme(p)` | build report dict, table sizes, README text |
| `manifest(p, data_type=None)` | every open-access file with bucket, key, md5 and HTTPS URL |
| `clinical(p, primary_only=True)` | cases joined to the primary diagnosis, with `os_time_days`, `os_event`, `age_at_diagnosis_years`, `stage` |
| `diagnoses`, `treatments`, `exposures`, `follow_ups`, `samples`, `aliquots` | the clinical tables as stored |
| `expression(p, genes, sample_type, kind, protein_coding, drop_par_y)` | genes x aliquots (TPM or counts), `attrs["gene_id"]` |
| `expression_samples(p)`, `star_qc(p)` | column to case mapping, STAR summary rows |
| `mutations(p, genes, non_silent, columns)` | MAF rows, filters pushed down |
| `mutation_frequency(p, non_silent=True)` | fraction of tumor barcodes mutated per gene |
| `copy_number(p, workflow, genes)` | genes x tumor aliquots, `attrs["workflow"]`, `attrs["gene_location"]` |
| `copy_number_workflows(p, kind)`, `copy_number_samples(p, workflow)` | available workflows, column mapping |
| `segments(p, kind, workflow)` | long segment tables |
| `mirna(p, kind, sample_type)`, `mirna_samples`, `mirna_isoforms` | miRNA matrices and tables |
| `rppa(p, sample_type)`, `rppa_samples`, `rppa_antibodies` | RPPA matrix (`attrs["agid"]`) and tables |
| `methylation(p, platform, probes, columns, sample_type)` | probes x aliquots betas |
| `methylation_platforms(p)`, `methylation_samples(p, platform)` | platforms, column mapping |
| `case_for_barcode(barcode)` | case id for TCGA and TARGET barcodes |
| `fetch_raw(row_or_project, file_id, dest)` | one raw GDC file, md5-verified |

### `scigantic_nci.idc`

| function | returns |
|---|---|
| `collections()`, `collection_ids()` | one row per collection (cached), sorted ids |
| `report(c)`, `tables(c)`, `readme(c)` | build report dict, table sizes, README text |
| `series(c, modality, patient, max_size_mb, columns)` | series.parquet, filters pushed down |
| `studies(c)`, `patients(c)` | studies.parquet, patients.parquet |
| `slides(c)` | sm_series joined to series, list columns flattened, `stain` column |
| `analysis_results(c)` | analysis_results.parquet or an empty frame |
| `clinical_tables(c)`, `clinical(c, table)`, `clinical_dictionary(c)` | clinical tables and their dictionary |
| `samples(c)`, `sample_files(c, modality)` | sample series metadata with `license_short_name`, .dcm paths |
| `read_sample(c, modality)` | `Volume` or `SlideLevels` (needs the `dicom` extra) |
| `assemble_level(path)` | one WSI level as an array |
| `pull_series(row_or_uuid, dest, bucket)` | local .dcm paths from the raw bucket |
| `series_size(row_or_uuid)`, `viewer_url(row)` | MB, IDC viewer link |

`NciError` is the base exception; `NotMirroredError` (a subclass) means the project, collection or table is not in the mirror and names what is.

## Where the bytes come from

A project or collection is resolved to a place to read from, in this order:

1. an explicit `root` argument (a local directory holding that prefix's files),
2. `$SCIGANTIC_NCI_ROOT/<bucket>/<key>` or `$SCIGANTIC_NCI_ROOT/<key>` if that directory exists (local mirrors, tests),
3. the Scigantic notebook mount at `$SCIGANTIC_MOUNT_PATH` (default `/mnt/archive`) when its `BUILD_REPORT.json` names the same project or collection,
4. the public bucket, read anonymously through pyarrow's S3 filesystem.

`gdc.projects()` and `idc.collections()` list the ids the same way (bucket listing plus any local copy), then read their rows from one index file at the bucket root, `s3://scigantic-gdc-open/PROJECTS.parquet` and `s3://scigantic-idc-open/COLLECTIONS.parquet`, instead of one `BUILD_REPORT.json` per id (and, for GDC, one `clinical/cases.parquet` per project). An id the index does not cover (a project added after the index was written, or one read from `$SCIGANTIC_NCI_ROOT` or the mount) is built from its own files, and when the index is absent, unreadable or lacks a column every row is built that way, 16 at a time. A local copy of the index at `$SCIGANTIC_NCI_ROOT/<bucket>/PROJECTS.parquet` or `$SCIGANTIC_NCI_ROOT/PROJECTS.parquet` (likewise `COLLECTIONS.parquet`) is used in place of the bucket's. The index files are objects at the bucket root, not prefixes, so older clients, which take the bucket's top-level directories as ids, do not see them. The bucket owners regenerate the index after any rebuild with `scigantic-nci build-index gdc --out PROJECTS.parquet` (or `idc --out COLLECTIONS.parquet`), which ignores any existing index, and upload it to the bucket root; until they do, `projects()` and `collections()` can show a rebuilt project's or collection's previous `built_at` and counts. The ids are always current because they come from the listing.

Anonymous S3 works because both mirror buckets are public-read; no credentials are required, and none are used even if present. The two exceptions are the raw NCI buckets: `idc.pull_series` reads `idc-open-data*` unsigned, and `gdc.fetch_raw` uses a signed default boto3 client because the newer GDC objects refuse unsigned reads, falling back to the GDC HTTPS download URL when S3 answers AccessDenied, 403 or 404 or no credentials are configured.

## Limits

- **No controlled-access data.** Only the GDC open-access tier is mirrored; there are no BAMs, no germline calls and no controlled clinical supplements.
- **Slides are not mirrored.** GDC `Slide Image` files appear in `manifest()` but not as tables (TCGA-BLCA alone lists 926 of them, 833 GB). IDC slide microscopy is in `sm_series` and `slides()`, but the mirror keeps only the thumbnail and the two lowest pyramid levels of one slide per collection; full pyramids come through `pull_series`.
- **Methylation is optional and large.** The beta matrices are built only on request, so `methylation_platforms()` is empty for most projects and `methylation()` raises `NotMirroredError`. Where present, a 450k matrix is 486,427 probes wide (486,427 x 45 for TCGA-CHOL); a 1,200-aliquot project read whole is about 2.3 GB in pandas, so pass `probes`, `columns` or `sample_type`.
- **Copy-number samples table.** In `copy_number_samples()` the `sample_type` describes the matrix column's own aliquot, and every column is the tumor member of the tumor/normal pair the file was called on (in TCGA-BLCA all 392 ASCAT3 columns are primary tumor barcodes and all 392 rows say Primary Tumor). The paired workflows (ASCAT2, ASCAT3, AscatNGS) keep both aliquots of the pair pipe-joined in `aliquot_submitter_id`, in no fixed order. Tables built before 2026-09-13 reported the first aliquot of that pair instead, so 195 BLCA rows said Blood Derived Normal; the current tables do not. Every matrix column is a single tumor aliquot.
- **Gene names are not unique.** GENCODE v36 gives 67 names to more than one gene id (plus 44 `_PAR_Y` duplicates, dropped by default); `expression()` and `copy_number()` keep the first and list the rest in `attrs["dropped_duplicate_names"]`.
- **Compressed DICOM needs a decoder.** Slide levels are JPEG baseline, and series pulled from IDC can be JPEG 2000 (tomosynthesis) or JPEG-LS; without a pydicom plugin `read_sample()` raises `NciError` with the install hint. The `dicom` extra installs them.
- **Eight collections have no sample.** `b_mode_and_ceus_liver` (120 US series) has none small enough to copy, and seven are entirely CC BY-NC, which the mirror never redistributes (`breast_cancer_screening_dbt`, `midrc_ricord_1a`, `midrc_ricord_1b`, `midrc_ricord_1c`, `nsclc_radiomics_genomics`, `nsclc_radiomics_interobserver1`, `phantom_fda`). `samples()` is empty and `read_sample()` raises `ValueError`; the index tables are complete and `pull_series` fetches any series under its license. Mirrored samples are only ever CC BY or NLM-terms series: 321 samples in 168 collections (190 CC BY 4.0, 129 CC BY 3.0, 2 NLM terms).
- **Which GDC projects are mirrored.** All 57 open-access projects of GDC Data Release 46.0 are built and live; `gdc.projects()` is the source of truth, and an unknown id raises `NotMirroredError` listing the mirrored ones.

## What is not here

- Controlled-access GDC data of any kind.
- Raw sequencing reads (BAM, FASTQ) and the per-file TSVs themselves; `fetch_raw()` gets one when needed.
- Stranded RNA-Seq counts and FPKM (only `counts_unstranded` and `tpm_unstranded` are stacked; the rest stays in the source TSVs listed in the manifest).
- GDC slide images, and IDC pixel data beyond the per-collection samples.
- Methylation for projects built without it.
- Anything derived: no normalisation, batch correction, survival models or segmentations that are not in the source data.

## Command line

`scigantic-nci` wraps the same functions; tables go to stdout as CSV, or to `--out FILE.csv`. `--root` maps to the `root=` argument. An id with a lowercase letter is an IDC collection, otherwise a GDC project. Library errors print to stderr with exit code 2.

```console
$ scigantic-nci projects
$ scigantic-nci collections
$ scigantic-nci tables TCGA-BLCA
$ scigantic-nci readme tcga_lihc
$ scigantic-nci clinical TCGA-BLCA --out blca_clinical.csv
$ scigantic-nci expression TCGA-BLCA --genes TP53,ESR1 --sample-type "Primary Tumor" --counts
$ scigantic-nci mutations TCGA-BLCA --genes TP53,RB1 --non-silent
$ scigantic-nci mutation-frequency TCGA-BLCA --top 20
$ scigantic-nci series tcga_lihc --modality CT --max-size-mb 50
$ scigantic-nci sample tcga_lihc --modality CT
$ scigantic-nci viewer-url tcga_lihc --series-uid 1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335
$ scigantic-nci pull-series tcga_lihc --series-uid 1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335 --dest /tmp/ct
```

Real output of each command is in [docs/cli.md](docs/cli.md).

## Testing

The suite runs against a complete local copy of the mirrors when `SCIGANTIC_NCI_ROOT` points at one (`tests/conftest.py` skips otherwise) and asserts measured numbers: TCGA-BLCA has 412 cases, 431 RNA aliquots and 19 solid tissue normals; tcga_lihc has 3,457 series and 377 patients. Tests marked `network` read the public buckets (`pytest -m "not network"` to skip them). CI runs Python 3.10 through 3.14 plus `mypy --strict`.

## License

MIT-0 for the code. See [Data provenance and citation](#data-provenance-and-citation) for the data.
