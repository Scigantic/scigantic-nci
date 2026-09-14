# Changelog

## 0.1.1 (2026-09-14)

- Documentation fix for `copy_number_samples`. The `*_samples.parquet` tables of the paired copy-number workflows (ASCAT2, ASCAT3, AscatNGS) were rebuilt on 2026-09-13 so that `sample_type`, `sample_submitter_id` and `case_submitter_id` describe the matrix column's own aliquot, the tumor member of the pair; the pair itself stays pipe-joined in `aliquot_submitter_id`. The 0.1.0 README, docs and docstring still described the earlier tables, where `sample_type` was the pair's first aliquot (195 of 392 TCGA-BLCA ASCAT3 rows said Blood Derived Normal), and told users to distrust the column. Checked against the mirror: 138 tables in 42 projects, 37,826 rows, no row whose `sample_type` differs from `clinical/aliquots` for the column's aliquot. The 82 AscatNGS columns in CGCI-BLGSP and HCMI-CMDC that are themselves pipe-joined pairs are noted as the remaining exception.
- The test that pinned the old counts now asserts the corrected invariant.
- `gdc.clinical` now takes `case_submitter_id` from the cases table for every row. It used to come from the diagnoses side of the left join, so a case with no primary diagnosis row had `case_submitter_id` NaN although its `case_id` and `submitter_id` were filled (TCGA-OV 21 of 608 rows, TCGA-LUAD 63 of 585, TCGA-BRCA 2 of 1098, TCGA-BLCA 1 of 412, CPTAC-3 24 of 1866; 1,497 rows across the 57 projects). `clinical()` is now one row per case with a non-null, unique `case_submitter_id`; `os_time_days`, `os_event`, `age_at_diagnosis_years` and `stage` are unchanged.
- The GDC mirror is complete: the remaining 21 projects finished building later on 2026-09-13, so all 57 open-access projects of GDC Data Release 46.0 are live and `gdc.projects()` returns 57 rows. The README bullet that said 36 of 57 were built is replaced.

## 0.1.0 (2026-09-13)

First release.

- `scigantic_nci.gdc`: per-project GDC tables as pandas DataFrames: `projects`, `report`, `tables`, `readme`, `manifest`, `clinical` (cases joined to the primary diagnosis with `os_time_days`, `os_event`, `age_at_diagnosis_years`, `stage`), the six clinical tables, `expression` and `expression_samples`, `star_qc`, `mutations` and `mutation_frequency`, `copy_number`, `copy_number_workflows`, `copy_number_samples`, `segments`, `mirna`, `mirna_samples`, `mirna_isoforms`, `rppa`, `rppa_samples`, `rppa_antibodies`, `methylation`, `methylation_platforms`, `methylation_samples`, `case_for_barcode`, `fetch_raw`. Gene, probe and sample subsets are pushed down to parquet.
- `scigantic_nci.idc`: per-collection IDC indexes and samples: `collections`, `collection_ids`, `report`, `tables`, `readme`, `series`, `studies`, `patients`, `slides`, `analysis_results`, `clinical_tables`, `clinical`, `clinical_dictionary`, `samples`, `sample_files`, `read_sample` (`Volume` for radiology, `SlideLevels` for slides), `assemble_level`, `pull_series`, `series_size`, `viewer_url`. Reading DICOM needs the `dicom` extra.
- `scigantic_nci._store`: one resolver for an explicit `root`, `$SCIGANTIC_NCI_ROOT`, the notebook mount, or anonymous S3.
- `scigantic-nci` console script: `projects`, `collections`, `tables`, `readme`, `clinical`, `expression`, `mutations`, `mutation-frequency`, `series`, `pull-series`, `sample`, `viewer-url`; CSV to stdout or `--out`.
- Mirrors at this release: 36 GDC projects (GDC Data Release 46.0) and 176 IDC collections (IDC data v24), both built 2026-09-13.
