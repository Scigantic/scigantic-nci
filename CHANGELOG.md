# Changelog

## 0.1.0 (2026-09-13)

First release.

- `scigantic_nci.gdc`: per-project GDC tables as pandas DataFrames: `projects`, `report`, `tables`, `readme`, `manifest`, `clinical` (cases joined to the primary diagnosis with `os_time_days`, `os_event`, `age_at_diagnosis_years`, `stage`), the six clinical tables, `expression` and `expression_samples`, `star_qc`, `mutations` and `mutation_frequency`, `copy_number`, `copy_number_workflows`, `copy_number_samples`, `segments`, `mirna`, `mirna_samples`, `mirna_isoforms`, `rppa`, `rppa_samples`, `rppa_antibodies`, `methylation`, `methylation_platforms`, `methylation_samples`, `case_for_barcode`, `fetch_raw`. Gene, probe and sample subsets are pushed down to parquet.
- `scigantic_nci.idc`: per-collection IDC indexes and samples: `collections`, `collection_ids`, `report`, `tables`, `readme`, `series`, `studies`, `patients`, `slides`, `analysis_results`, `clinical_tables`, `clinical`, `clinical_dictionary`, `samples`, `sample_files`, `read_sample` (`Volume` for radiology, `SlideLevels` for slides), `assemble_level`, `pull_series`, `series_size`, `viewer_url`. Reading DICOM needs the `dicom` extra.
- `scigantic_nci._store`: one resolver for an explicit `root`, `$SCIGANTIC_NCI_ROOT`, the notebook mount, or anonymous S3.
- `scigantic-nci` console script: `projects`, `collections`, `tables`, `readme`, `clinical`, `expression`, `mutations`, `mutation-frequency`, `series`, `pull-series`, `sample`, `viewer-url`; CSV to stdout or `--out`.
- Mirrors at this release: 36 GDC projects (GDC Data Release 46.0) and 176 IDC collections (IDC data v24), both built 2026-09-13.
