# scigantic-nci: the console script

`scigantic-nci` wraps the functions in `scigantic_nci.gdc` and `scigantic_nci.idc`. Tables are written as CSV to stdout, or to `--out FILE.csv`; `readme`, `sample`, `viewer-url` and `pull-series` print text, and `build-index` writes a parquet file. The global `--root DIR` is passed as the `root=` argument (a local directory holding one project's or collection's files); without it the resolver uses `$SCIGANTIC_NCI_ROOT`, the notebook mount, or the public bucket. An id containing a lowercase letter (`tcga_lihc`) is an IDC collection, anything else (`TCGA-LIHC`) a GDC project.

Errors raised by the library (`NciError`, `NotMirroredError`, `ValueError`) are printed to stderr prefixed with `scigantic-nci:` and the exit code is 2.

Every output below was produced on 2026-09-13 against the local mirrors (`SCIGANTIC_NCI_ROOT` set to the GDC or IDC mirror); long lines are cut.

```
$ scigantic-nci --help
usage: scigantic-nci [-h] [--version] [--root ROOT] command ...

GDC projects and IDC collections from the Scigantic mirrors, as CSV.

positional arguments:
  command
    projects            one row per mirrored GDC project
    collections         one row per mirrored IDC collection
    build-index         write the catalog index file (PROJECTS.parquet or
                        COLLECTIONS.parquet) for a mirror bucket
    tables              tables of a GDC project or IDC collection with rows,
                        columns, bytes
    readme              print the README.md of a project or collection
    clinical            GDC cases joined to the primary diagnosis, with
                        survival columns
    expression          GDC STAR expression, genes x aliquots
    mutations           GDC masked somatic mutations, one row per call
    mutation-frequency  fraction of tumor samples mutated per gene
    series              IDC series.parquet, optionally one modality
    pull-series         download one series from the raw IDC bucket
                        (anonymous)
    sample              describe the sample series read_sample() would load
    viewer-url          IDC web viewer link for one series
```

## GDC

```
$ scigantic-nci projects | head -4
project,gdc_data_release,built_at,n_tables,cases,source_files_in_project,files_read,bytes_read
TARGET-ALL-P1,"Data Release 46.0 - August 10, 2026",2026-09-13T20:57:15.218280+00:00,13,24,18,13,50731894
TARGET-ALL-P2,"Data Release 46.0 - August 10, 2026",2026-09-13T20:58:03.683302+00:00,21,1587,2408,2397,3364435733
TARGET-ALL-P3,"Data Release 46.0 - August 10, 2026",2026-09-13T21:01:39.074492+00:00,17,191,820,394,597928155
```

```
$ scigantic-nci tables TCGA-BLCA | head -5
table,rows,columns,bytes
MANIFEST,10603,19,1474026
clinical/aliquots,4836,7,145554
clinical/cases,412,25,56581
clinical/diagnoses,1014,34,78641
```

`clinical` writes one row per case with the helper columns of `gdc.clinical()` (`os_time_days`, `os_event`, `age_at_diagnosis_years`, `stage`); `--all-diagnoses` keeps every diagnosis row.

```
$ scigantic-nci clinical TCGA-BLCA --out blca_clinical.csv      # 412 rows, 62 columns, nothing on stdout
```

`expression` needs `--genes`; `--sample-type` takes one type or a comma-separated list; `--counts` switches from TPM to raw STAR counts. Columns are aliquot barcodes.

```
$ scigantic-nci expression TCGA-BLCA --genes TP53,ESR1 --sample-type "Solid Tissue Normal" | cut -c1-110
gene_name,TCGA-BT-A20Q-11A-11R-A14Y-07,TCGA-GC-A3WC-11A-11R-A22U-07,TCGA-BL-A13J-11A-13R-A10U-07,TCGA-CU-A0YR-
ESR1,3.4256,2.814,1.3847,3.1993,0.9256,1.468,4.0321,28.9746,4.9208,5.8203,2.6042,2.2363,14.606,2.151,3.749,0.8
TP53,26.7296,28.1733,19.8971,34.7217,30.0243,23.9229,53.0271,43.7137,62.368,30.1476,32.3026,22.002,33.3175,18.
```

`mutations` prints the 17 default MAF columns (`gdc.DEFAULT_MAF_COLUMNS`); `--genes` and `--non-silent` are parquet row filters.

```
$ scigantic-nci mutations TCGA-BLCA --genes TP53,RB1 --non-silent | cut -d, -f3,7,13 | head -4
Hugo_Symbol,Variant_Classification,HGVSp_Short
TP53,Missense_Mutation,p.R158H
TP53,Frame_Shift_Ins,p.S240Kfs*24
RB1,Nonsense_Mutation,p.R320*
$ scigantic-nci mutations TCGA-BLCA --genes TP53,RB1 --non-silent | wc -l
     307
```

`mutation-frequency` prints the top 20 genes by default (`--top 0` for all); `--include-silent` counts silent calls too. The denominator is every tumor barcode in the MAF (414 for TCGA-BLCA).

```
$ scigantic-nci mutation-frequency TCGA-BLCA --top 5
Hugo_Symbol,fraction_mutated
TP53,0.48792270531400966
TTN,0.4251207729468599
KMT2D,0.2632850241545894
ARID1A,0.25120772946859904
KDM6A,0.25120772946859904
```

An unknown id lists the mirrored ones and exits 2:

```
$ scigantic-nci clinical TCGA-NOPE
scigantic-nci: 'TCGA-NOPE' is not a mirrored GDC project. Mirrored (36): TARGET-ALL-P1, TARGET-ALL-P2, TARGET-ALL-P3, TCGA-ACC, TCGA-BLCA, TCGA-BRCA, TCGA-CESC, ...
$ echo $?
2
```

## Catalog index

`build-index gdc|idc --out PATH` builds the `projects()` or `collections()` frame from every project's or collection's `BUILD_REPORT.json`, ignoring any existing index, and writes it as parquet. It is for the bucket owners: after any rebuild the result is uploaded to `s3://scigantic-gdc-open/PROJECTS.parquet` or `s3://scigantic-idc-open/COLLECTIONS.parquet`. Run on 2026-09-15 against the public buckets (`SCIGANTIC_NCI_ROOT` unset, no mount):

```
$ scigantic-nci build-index gdc --out PROJECTS.parquet
PROJECTS.parquet: 57 rows, 5772 bytes
$ scigantic-nci build-index idc --out COLLECTIONS.parquet
COLLECTIONS.parquet: 176 rows, 15648 bytes
```

## IDC

```
$ scigantic-nci collections | head -3
collection_id,collection_name,n_patients,n_studies,n_series,total_GB,modalities,licenses,clinical_tables,idc_data_version,built_at
4d_lung,4D-Lung,20,589,6690,178.75,"CT, RTSTRUCT",CC BY 3.0,0,v24,2026-09-13T17:44:13.574883+00:00
acrin_6698,ACRIN-6698,385,1123,18747,822.22,"MR, SEG",CC BY 4.0,1,v24,2026-09-13T17:44:34.185940+00:00
```

```
$ scigantic-nci tables tcga_lihc
table,rows,columns,bytes
analysis_results,2,13,16812
clinical/bamf_aimi_annotations_liver_ct_qa_results,98,15,16479
clinical/bamf_aimi_annotations_liver_mr_qa_results,72,15,13378
clinical/dictionary,30,4,4561
patients,377,7,9782
series,3457,31,308879
sm_series,870,20,43293
studies,614,9,29530
```

`series` prints series.parquet (31 columns); `--modality` is validated against the collection's modalities and `--max-size-mb` keeps series at or below that size.

```
$ scigantic-nci series tcga_lihc --modality CT --max-size-mb 1 | cut -d, -f3,4,12,24,31 | head -4
PatientID,SeriesInstanceUID,Modality,instanceCount,series_size_MB
TCGA-BC-A10W,1.3.6.1.4.1.14519.5.2.1.8421.4008.141230740284190533687726815345,CT,1,0.527958
TCGA-BC-A10W,1.3.6.1.4.1.14519.5.2.1.8421.4008.322726584798806176092185592280,CT,1,0.527966
TCGA-BC-A10W,1.3.6.1.4.1.14519.5.2.1.8421.4008.731035386785425769837043483583,CT,1,0.52796
$ scigantic-nci series tcga_lihc --modality XX
scigantic-nci: tcga_lihc has no 'XX' series; modalities: MR, SEG, SM, CT, PT
```

`sample` describes the sample series that `idc.read_sample(collection, modality)` would load, without needing pydicom: the BUILD_REPORT entry and the `.dcm` paths under the collection prefix.

```
$ scigantic-nci sample tcga_lihc --modality CT | head -16
modality: CT
kind: series
files: 36
bytes: 18957986
instanceCount_full_series: 36
series_size_MB_full: 18.957986
PatientID: TCGA-DD-A3A6
StudyInstanceUID: 1.3.6.1.4.1.14519.5.2.1.3344.4008.325199116162170045267705852537
SeriesInstanceUID: 1.3.6.1.4.1.14519.5.2.1.3344.4008.1590978269182606855431778873.8
series_aws_url: s3://idc-open-data/54e10b05-fd18-4886-81c9-765faeafb7a4/*
license: CC BY 3.0
license_short_name: CC BY 3.0
source_DOI: 10.7937/k9/tcia.2016.immqw8uq
folder: sample/CT_1.3.6.1.4.1.14519.5.2.1.3344.4008.1590978269182606855431778873.8
dcm_files: 36
  sample/CT_1.3.6.1.4.1.14519.5.2.1.3344.4008.1590978269182606855431778873.8/0ddfe116-0d02-4fc3-bc34-26744e742684.dcm
  sample/CT_1.3.6.1.4.1.14519.5.2.1.3344.4008.1590978269182606855431778873.8/1630ca3b-f591-4ae2-a24e-749e1d51673d.dcm
```

`viewer-url` and `pull-series` name a series by `--series-uid` (SeriesInstanceUID) or `--uuid` (crdc_series_uuid); the row is looked up in the collection's series.parquet so the right raw bucket is used.

```
$ scigantic-nci viewer-url tcga_lihc --series-uid 1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335
https://viewer.imaging.datacommons.cancer.gov/viewer/1.3.6.1.4.1.14519.5.2.1.3344.4008.164131928981776564121389657757?seriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335
```

```
$ scigantic-nci pull-series tcga_lihc --series-uid 1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335 --dest /tmp/ct
/tmp/ct/48c50a43-3e94-4b84-b228-041635474c53.dcm
$ ls -l /tmp/ct
-rw-r--r--  1 user  wheel  409608 Sep 13 17:11 48c50a43-3e94-4b84-b228-041635474c53.dcm
```

That is the smallest CT series of tcga_lihc (one 0.41 MB scout), read anonymously from `idc-open-data`. The destination directory and user name in this listing are shortened; the file name and size are as returned.
