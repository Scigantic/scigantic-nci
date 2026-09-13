# scigantic_nci.idc: Imaging Data Commons collections

Every IDC collection is mirrored at `s3://scigantic-idc-open/<collection_id>/` as a handful of
parquet indexes (series, studies, patients, slide attributes, analysis results, clinical
tables) plus one or more small sample DICOM series. `scigantic_nci.idc` reads those with
pandas; pixels come from the samples in the mirror or, on demand, from the raw IDC buckets
anonymously. Nothing here needs an account or BigQuery.

All snippets below were run against the mirror built 2026-09-13 (IDC data v24, 176
collections). Every number is what the code returned.

```python
import scigantic_nci as nci
```

## What is mirrored

```python
>>> c = nci.idc.collections(); len(c)
176
>>> c.sort_values('n_series', ascending=False).head(8)[['collection_id','n_patients','n_series','total_GB','modalities','clinical_tables']]
                   collection_id  n_patients  n_series  total_GB                         modalities  clinical_tables
90                          nlst       26410    590572  25737.06                    SR, CT, SEG, SM                5
74                         ispy2         719     32411   1705.31                            MR, SEG                1
62                        eay131        2813     30293    778.58  RTSTRUCT, CT, SEG, MR, PT, NM, XA                0
64                          gtex         971     25503   8354.19                                 SM                0
10   breast_cancer_screening_dbt        5060     22032   1599.43                                 MG                0
114                    prostatex         346     19177     16.17                        MR, SEG, SR                6
1                     acrin_6698       385     18747    822.22                            MR, SEG                1
38               covid_19_ny_sbu        1384     17950    499.49     CR, CT, DX, MR, NM, PT, OT, SR                1
>>> c[c.modalities.str.contains('SM')].shape[0], c.total_GB.sum().round(1)
73 96940.4
```

`collections()` reads the 176 `BUILD_REPORT.json` files (16 threads) once per process and
returns a copy of the cached frame. The `modalities` string lists a collection's modalities
by series count, most first; `clinical_tables` is the number of clinical tables (59
collections have at least one). `collection_ids()` gives just the sorted ids, and
`report(collection)` / `readme(collection)` return one collection's BUILD_REPORT dict and
README text.

An unknown id raises `NotMirroredError` listing the valid ones:

```python
>>> nci.idc.series('nope')
NotMirroredError: 'nope' is not a mirrored IDC collection (s3://scigantic-idc-open/nope/ is not in the mirror (no BUILD_REPORT.json)). Valid ids (176): 4d_lung, acrin_6698, acrin_contralateral_breast_mr, acrin_flt_brea ...
```

## Series, studies, patients

```python
>>> s = nci.idc.series('tcga_lihc'); s.shape
(3457, 31)
>>> s.Modality.value_counts()
Modality
MR     910
SEG    899
SM     870
CT     777
PT       1
```

`series()` filters are pushed down to pyarrow, so a filtered read of a large collection does
not build the full frame. `modality` is validated against the collection's modalities,
`patient` accepts one id or a list, `max_size_mb` caps `series_size_MB`, and `columns`
reads a subset of the 31 columns.

```python
>>> ct = nci.idc.series('tcga_lihc', modality='CT', max_size_mb=50)
>>> ct[['PatientID','SeriesDescription','instanceCount','series_size_MB','crdc_series_uuid']].head(4)
      PatientID SeriesDescription  instanceCount  series_size_MB                      crdc_series_uuid
0  TCGA-BC-A10W         LOCALIZER              1        0.527958  3b5c3246-2f5f-423d-931e-64576cf73df3
1  TCGA-BC-A10W  150cc OMNIPAQUE)             34       17.951720  2a6025a1-0941-48f8-88ee-ff7c05be29df
2  TCGA-BC-A10W  150cc OMNIPAQUE)              1        0.527966  e909b395-f9fc-461d-86c3-a64608e4ccf4
3  TCGA-BC-A10W  150cc OMNIPAQUE)              60       31.679832  0f467e8d-0d55-4b09-aa2a-905575968a1e
>>> nci.idc.series('tcga_lihc', modality='XX')
ValueError: tcga_lihc has no 'XX' series; modalities: MR, SEG, SM, CT, PT
```

```python
>>> nci.idc.studies('tcga_lihc').head(3)
                                                   StudyInstanceUID     PatientID   StudyDate                        StudyDescription modalities  n_series  instances     size_MB body_parts
0  1.3.6.1.4.1.14519.5.2.1.1079.4008.108664975487956635580635879165  TCGA-G3-A25T  2001-06-26          MR ABDOMEN ENHANCED & NON-BODY         MR        21       1087  188.479384      LIVER
1  1.3.6.1.4.1.14519.5.2.1.1079.4008.110146226557974939388520158381  TCGA-G3-AAV1  2007-01-25  MR ABDOMEN NONENHANCED & ENHANCED-BODY    MR, SEG        24       1573  322.447554      LIVER
2  1.3.6.1.4.1.14519.5.2.1.1079.4008.116712770482771075241626231125  TCGA-G3-A3CH  2004-05-07       CT Abdomen Nonenh & Enhanced-BODY    CT, SEG         9        664  355.920054      LIVER
>>> nci.idc.patients('tcga_lihc').head(3)
      PatientID PatientSex PatientAge  n_studies  n_series modalities     size_MB
0  TCGA-2V-A95S                  None          1         4    SEG, SM  492.470990
1  TCGA-2Y-A9GS                  None          1         4    SEG, SM  947.391868
2  TCGA-2Y-A9GT                  None          1         4    SEG, SM  630.389374
```

tcga_lihc has 614 studies and 377 patients. Pathology-only patients (SEG + SM) carry an
empty PatientSex and no PatientAge; the DICOM slides do not record them.

## Slides

`slides()` joins `sm_series.parquet` to `series.parquet` on SeriesInstanceUID. The four
list-valued code columns (embedding medium, fixative, stain, and their code strings) are
flattened to comma-joined strings, and `stain` holds the flattened
`staining_usingSubstance_CodeMeaning`.

```python
>>> sl = nci.idc.slides('tcga_lihc'); sl.shape
(870, 51)
>>> sl[['ContainerIdentifier','stain','ObjectiveLensPower','max_TotalPixelMatrixColumns','series_size_MB']].head(3)
       ContainerIdentifier                                         stain  ObjectiveLensPower  max_TotalPixelMatrixColumns  series_size_MB
0  TCGA-ED-A459-01Z-00-DX1  hematoxylin stain, water soluble eosin stain                  40                       114240     1613.738424
1  TCGA-CC-5260-01Z-00-DX1  hematoxylin stain, water soluble eosin stain                  40                        52360      674.268904
2  TCGA-DD-AACL-01Z-00-DX1  hematoxylin stain, water soluble eosin stain                  40                        59760      466.882722
```

All 870 tcga_lihc slides are H&E; the widest is 209,439 px. A collection without slide
microscopy raises `ValueError` naming its modalities.

## Analysis results and clinical tables

```python
>>> nci.idc.analysis_results('tcga_lihc')[['analysis_result_id','subjects','modalities','license_short_name']]
      analysis_result_id  subjects modalities license_short_name
0  bamf_aimi_annotations      4226        SEG          CC BY 4.0
1      tcga_sbu_til_maps      7600        SEG          CC BY 4.0
```

`analysis_results()` returns an empty 13-column frame for the 121 collections without one.
The `analysis_result_id` column of `series.parquet` links each SEG series back to these rows.

```python
>>> nci.idc.clinical_tables('acrin_6698')
['acrin_6698_clinical']
>>> clin = nci.idc.clinical('acrin_6698'); clin.shape
(385, 32)
>>> clin.columns.tolist()[:12]
['dicom_patient_id', 'source_batch', 'i_spy_2_research_id', 'tcia_patient_id', 't0', 't1', 't2', 't3', 'bmmr2_train', 'bmmr2_test', 'primary_aim_t0', 'primary_aim_t1']
>>> nci.idc.clinical_dictionary('acrin_6698').iloc[4:7]
      short_table_name column                                         column_label                                                                                                    values
4  acrin_6698_clinical     t0       T0 (baseline) MRI study included in collection  [{'option_code': '0', 'option_description': ' No'}\n {'option_code': '1', 'option_description': ' Yes'}]
5  acrin_6698_clinical     t1       T1 (early-Tx) MRI study included in collection  [{'option_code': '0', 'option_description': ' No'}\n {'option_code': '1', 'option_description': ' Yes'}]
6  acrin_6698_clinical     t2  T2 (inter-regimen) MRI study included in collection  [{'option_code': '0', 'option_description': ' No'}\n {'option_code': '1', 'option_description': ' Yes'}]
```

`clinical()` defaults to the first table in BUILD_REPORT; the dictionary explains each
column and its coded values (the `values` column is a string, not parsed). Clinical tables
join to `patients.parquet` through `dicom_patient_id` (the DICOM PatientID).

## Sample DICOM in the mirror

Each collection keeps whole small radiology series and, for slides, the thumbnail plus the
two lowest pyramid levels of one slide, under `sample/<Modality>_<SeriesInstanceUID>/`.

```python
>>> nci.idc.samples('tcga_lihc')[['modality','kind','files','bytes','instanceCount_full_series','PatientID']]
  modality    kind  files     bytes  instanceCount_full_series     PatientID
0       CT  series     36  18957986                         36  TCGA-DD-A3A6
1       MR  series    111   6561042                        111  TCGA-G3-A3CJ
2       PT  series    215  12890310                        215  TCGA-DD-A4NG
>>> nci.idc.sample_files('tcga_lihc', 'CT')[:2]
['sample/CT_1.3.6.1.4.1.14519.5.2.1.3344.4008.1590978269182606855431778873.8/0ddfe116-0d02-4fc3-bc34-26744e742684.dcm', 'sample/CT_1.3.6.1.4.1.14519.5.2.1.3344.4008.1590978269182606855431778873.8/1630ca3b-f591-4ae2-a24e-749e1d51673d.dcm']
```

`read_sample()` needs pydicom (`pip install "scigantic-nci[dicom]"`). Radiology samples
come back as a `Volume`: float32 array `(n_instances, rows, cols)`, instances sorted by
ImagePositionPatient z then InstanceNumber, CT rescaled to Hounsfield units, and
`spacing` as (slice, row, column) in mm. A multi-frame instance contributes its middle
frame. If the instances have different shapes, `array` is None, `ragged` is True and
`images` holds the 2D arrays.

```python
>>> vol = nci.idc.read_sample('tcga_lihc', 'CT'); vol
Volume(modality='CT', description='AXIAL', shape=(36, 512, 512), spacing=(5.0, 0.585938, 0.585938), n_instances=36)
>>> vol.array.shape, vol.array.dtype, vol.array.min(), vol.array.max()
(36, 512, 512) float32 -1000.0 1323.0
>>> [float(ds.ImagePositionPatient[2]) for ds in vol.datasets[:4]]
[-174.900009, -169.900009, -164.900009, -159.900009]
>>> nci.idc.read_sample('tcga_lihc', 'MR')
Volume(modality='MR', description='DWI LIVER   B50, B500', shape=(111, 144, 192), spacing=(6.5, 1.9791666269302, 1.9791666269302), n_instances=111)
```

The CT slice spacing (5.0 mm from the positions) is smaller than its SliceThickness
(6.0 mm): the slices overlap. The MR sample is a two b-value DWI, so pairs of instances
share a position and are ordered by InstanceNumber within it.

Slide samples come back as `SlideLevels`, lowest resolution first; `level(i)` assembles
that level's tiles row-major (`ncol = ceil(TotalPixelMatrixColumns / Columns)`) and crops
to the full pixel matrix. `assemble_level(path)` does the same for any WSI instance on disk.

```python
>>> slide = nci.idc.read_sample('tcga_chol'); slide.description, slide.widths
Frozen HE TP TSA [708, 1992, 7968]
>>> for lv in slide.levels: print(lv.index, lv.image_type, lv.width, lv.height, lv.rows, lv.cols, lv.n_frames, lv.um_per_px)
0 THUMBNAIL 708 768 768 708 1 11.366779661017
1 VOLUME 1992 2158 240 240 81 4.04
2 VOLUME 7968 8633 240 240 1224 1.01
>>> slide.level(0).shape
(768, 708, 3)
>>> slide.level(1).shape
NciError: cannot decode pixel data (transfer syntax 1.2.840.10008.1.2.4.50): Unable to decompress 'JPEG Baseline (Process 1)' pixel data because all plugins are missing dependencies:
	gdcm - requires gdcm>=3.0.10
	pylibjpeg - requires pylib ...
```

The thumbnail is stored uncompressed; the VOLUME levels are JPEG baseline and need a
pydicom decoder plugin (`pip install "pylibjpeg[all]"`, or Pillow). With one installed,
`slide.level(1)` is `(2158, 1992, 3)` uint8.

Fluorescence slides have no channel axis. The htan_ohsu sample is a multiplex IHC slide
whose three "levels" are three single-tile instances of the same 143 x 103 px size, one
per channel, not a pyramid:

```python
>>> fl = nci.idc.read_sample('htan_ohsu'); fl.description, [(lv.height, lv.width, lv.n_frames) for lv in fl.levels]
mIHC [(103, 143, 1), (103, 143, 1), (103, 143, 1)]
>>> fl.level(0).shape, fl.level(0).dtype
(103, 143) uint8
```

Digital breast tomosynthesis stores a whole reconstruction in one multi-frame instance
(21 frames of 2457 x 1996 here); `read_sample` takes the middle frame, and JPEG 2000 frames
also need a decoder plugin:

```python
>>> nci.idc.samples('breast_cancer_screening_dbt')[['modality','files','bytes','series_aws_url']]
  modality  files    bytes                                                series_aws_url
0       MG      1  2322508  s3://idc-open-data-cr/8217b99a-d12f-4de7-914d-1f762dc38006/*
>>> nci.idc.read_sample('breast_cancer_screening_dbt')
NciError: cannot decode pixel data (transfer syntax 1.2.840.10008.1.2.4.90): Unable to decompress 'JPEG 2000 Image Compression (Lossless Only)' pixel data because all plugins are missing dependencies:
	gdcm - requires gdcm>=3.0.10
	pylibjpe ...
```

One collection, `b_mode_and_ceus_liver` (120 US series), has no sample; `samples()` is
empty and `read_sample()` raises `ValueError`.

## Pulling a full series from the raw IDC bucket

Raw objects live at `s3://<aws_bucket>/<crdc_series_uuid>/<instance>.dcm` and accept
anonymous reads. `pull_series` takes a `series.parquet` row (or a uuid string, bucket
defaulting to `idc-open-data`), lists the prefix, skips the folder marker, downloads with
8 threads and returns the local paths. `series_size` and `viewer_url` work from the row.

```python
>>> row = ct.sort_values('series_size_MB').iloc[0]; nci.idc.series_size(row), nci.idc.viewer_url(row)
0.409608
https://viewer.imaging.datacommons.cancer.gov/viewer/1.3.6.1.4.1.14519.5.2.1.3344.4008.164131928981776564121389657757?seriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335
>>> files = nci.idc.pull_series(row, dest); files
['48c50a43-3e94-4b84-b228-041635474c53.dcm'] (3.5 s)
>>> ds = pydicom.dcmread(files[0]); ds.Modality, ds.SeriesDescription, ds.pixel_array.shape
CT SCOUT (395, 512)
```

Note the bucket varies by collection: tcga_lihc is in `idc-open-data`, while
breast_cancer_screening_dbt (CC BY-NC) is in `idc-open-data-cr`. Always pass the row so
`aws_bucket` is used.

## Reading from S3 instead of a local copy

With no local mirror or mount, the same calls read the public bucket anonymously.
A column subset of `acrin_6698` series, filtered to SEG:

```python
>>> s = nci.idc.series('acrin_6698', modality='SEG', columns=['SeriesInstanceUID','PatientID','series_size_MB']); s.shape, s.PatientID.nunique()
(2213, 3) 385      # 2.3 s from S3
```

## Function reference

| Function | Returns |
| --- | --- |
| `collections()` | one row per collection (cached) |
| `collection_ids()` | sorted list of ids |
| `report(c)`, `tables(c)`, `readme(c)` | BUILD_REPORT dict, table sizes, README text |
| `series(c, modality, patient, max_size_mb, columns)` | series.parquet, filtered |
| `studies(c)`, `patients(c)` | studies.parquet, patients.parquet |
| `slides(c)` | sm_series joined to series, lists flattened, `stain` column |
| `analysis_results(c)` | analysis_results.parquet or empty frame |
| `clinical_tables(c)`, `clinical(c, table)`, `clinical_dictionary(c)` | clinical tables |
| `samples(c)`, `sample_files(c, modality)` | sample series metadata, .dcm paths |
| `read_sample(c, modality)` | `Volume` or `SlideLevels` (pydicom) |
| `assemble_level(path)` | one WSI level as an array |
| `pull_series(row_or_uuid, dest, bucket)` | local .dcm paths from the raw bucket |
| `series_size(row_or_uuid)`, `viewer_url(row)` | MB, IDC viewer link |

Every per-collection function takes the collection id first and `root: str | None = None`
last (a local directory holding that collection's files).
