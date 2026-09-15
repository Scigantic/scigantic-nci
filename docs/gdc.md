# scigantic_nci.gdc: GDC projects as pandas tables

Every function takes the GDC project id first (`"TCGA-BLCA"`) and an optional `root=` last. Without `root` the resolver reads `$SCIGANTIC_NCI_ROOT`, a mounted Scigantic archive, or `s3://scigantic-gdc-open/<PROJECT>/` anonymously. Every snippet below was run against the Data Release 46.0 build of the mirror; the outputs are pasted as printed.

```python
from scigantic_nci import gdc
```

## What is mirrored

```python
gdc.projects()[["project", "cases", "n_tables", "files_read", "bytes_read"]].head(6)
```
```
      project  cases  n_tables  files_read  bytes_read
TARGET-ALL-P1     24        13          13    50731894
     TCGA-ACC     92        34        1387  1549135889
    TCGA-BLCA    412        35        6529  6426110741
    TCGA-BRCA   1098        35       16612 17431911368
    TCGA-CESC    307        35        4374  4635989112
    TCGA-CHOL     51        34         698   730495995
```

`projects()` lists the mirrored projects, reads their rows from the bucket's `PROJECTS.parquet` index (written from each project's `BUILD_REPORT.json` plus the row count of `clinical/cases`), builds any row the index lacks from those two files, and caches the result for the process. `report(project)` returns the raw build report, `readme(project)` the README text, and `tables(project)` the per-table sizes:

```python
gdc.tables("TCGA-BLCA").head(8)
```
```
              table  rows  columns   bytes
           MANIFEST 10603       19 1474026
  clinical/aliquots  4836        7  145554
     clinical/cases   412       25   56581
 clinical/diagnoses  1014       34   78641
 clinical/exposures   677       17   47273
clinical/follow_ups  1967       17   86965
   clinical/samples  1240       16   70033
clinical/treatments  2306       29  125282
```

An unknown project lists the valid ids:

```
>>> gdc.clinical("TCGA-NOPE")
NotMirroredError: 'TCGA-NOPE' is not a mirrored GDC project. Mirrored (34): TARGET-ALL-P1, TCGA-ACC, TCGA-BLCA, TCGA-BRCA, TCGA-CESC, TCGA-CHOL, TCGA-COAD, TCGA-DLBC, ...
```

## Clinical

`clinical(project)` joins `clinical/cases` to the one diagnosis row flagged `diagnosis_is_primary_disease` (columns that exist in both tables get `_dx` on the diagnosis side) and appends four helper columns:

- `os_time_days`: `demographic_days_to_death` when vital status is Dead, else `days_to_last_follow_up`
- `os_event`: True when vital status is Dead
- `age_at_diagnosis_years`: `age_at_diagnosis` (GDC stores days) / 365.25
- `stage`: a copy of the best-covered `*_stage` column, falling back to `tumor_grade`; the source is in `df.attrs["stage_column"]`

```python
c = gdc.clinical("TCGA-BLCA")
c[["submitter_id", "demographic_sex_at_birth", "demographic_vital_status",
   "os_time_days", "os_event", "age_at_diagnosis_years", "stage"]].head()
```
```
submitter_id demographic_sex_at_birth demographic_vital_status  os_time_days  os_event  age_at_diagnosis_years     stage
TCGA-DK-A3IS                     male                    Alive        1529.0     False               68.714579  Stage II
TCGA-CF-A1HS                   female                    Alive         382.0     False               75.468857 Stage III
TCGA-BT-A20V                   female                     Dead         154.0      True               59.003422  Stage IV
TCGA-XF-AAN8                   female                     Dead         118.0      True               74.767967 Stage III
TCGA-GC-A3WC                   female                    Alive         540.0     False               80.205339 Stage III
```
```
>>> c.shape, c.attrs["stage_column"], c.attrs["stage_coverage"]
((412, 62), 'ajcc_pathologic_stage', {'ajcc_pathologic_stage': 409})
>>> c["stage"].value_counts().to_dict()
{'Stage III': 141, 'Stage IV': 136, 'Stage II': 131, 'Stage I': 1}
```

412 cases, 182 events, 50 cases with no survival time (no death date and no follow-up on the primary diagnosis row). The stage column is whatever the disease uses:

```
TCGA-ACC:  stage_column='ensat_pathologic_stage' coverage={'ajcc_pathologic_stage': 0, 'ensat_pathologic_stage': 90, 'figo_stage': 0}
TCGA-OV:   stage_column='figo_stage'             coverage={'figo_stage': 582}
TCGA-THYM: stage_column='masaoka_stage'          coverage={'masaoka_stage': 122, 'ajcc_pathologic_stage': 0}
TCGA-LGG:  stage_column='tumor_grade'            coverage={}
TCGA-GBM:  stage_column=None                     coverage={}
```

`clinical(project, primary_only=False)` keeps every diagnosis row (1,014 rows for BLCA's 412 cases). The other clinical tables are thin wrappers: `diagnoses`, `treatments`, `exposures`, `follow_ups`, `samples`, `aliquots`. A table the project does not have raises `NotMirroredError` naming what is there (TCGA-CHOL has no `clinical/exposures`).

## RNA-Seq expression

Matrix columns are aliquot barcodes; `expression_samples` maps them to cases and sample types:

```python
s = gdc.expression_samples("TCGA-BLCA")
s["sample_type"].value_counts()
```
```
Primary Tumor          412
Solid Tissue Normal     19
```

`expression` picks the columns from that table first and reads only those, so a `sample_type` subset never touches the other aliquots. `genes` are symbols matched on `gene_name`; missing symbols are reported rather than raised (all missing raises `ValueError`).

```python
e = gdc.expression("TCGA-BLCA", genes=["TP53", "ESR1", "ERBB2", "NOPE"], sample_type="Primary Tumor")
e.iloc[:, :3].round(2)
```
```
           TCGA-FD-A3NA-01A-11R-A21D-07  TCGA-DK-A6B2-01A-11R-A30C-07  TCGA-GV-A3QG-01A-11R-A220-07
gene_name                                                                                          
ESR1                           0.210000                      0.640000                      0.970000
TP53                          13.460000                     89.970001                     53.820000
ERBB2                          97.309998                    136.850006                     39.869999
```
```
>>> e.shape, e.attrs["missing_genes"], dict(e.attrs["gene_id"])
((3, 412), ['NOPE'], {'ESR1': 'ENSG00000091831.24', 'TP53': 'ENSG00000141510.18', 'ERBB2': 'ENSG00000141736.14'})
```

`kind="counts"` gives raw STAR counts (int32) instead of TPM (float32); `protein_coding=True` keeps the 19,938 protein-coding names. The 44 `_PAR_Y` ids are dropped by default. GENCODE v36 also gives 67 further names to more than one gene id (1,189 extra rows: `Y_RNA`, `U6`, `SNORD116`, `5S_rRNA`, ...); the first occurrence is kept and the names are listed in `attrs["dropped_duplicate_names"]`, so the full matrix is 59,427 x 431 for BLCA.

```
>>> full = gdc.expression("TCGA-BLCA")
>>> full.shape, len(full.attrs["dropped_duplicate_names"]), full.attrs["dropped_duplicate_names"][:4]
((59427, 431), 67, ['5S_rRNA', '5_8S_rRNA', '7SK', 'ACTL10'])
```

Timings on TCGA-BRCA (1,231 aliquots), local mirror, median of 3:

| call | shape | time | pandas memory |
|---|---|---|---|
| `expression("TCGA-BRCA", genes=[5 symbols])` | (5, 1231) | 0.10 s | 0.05 MB |
| `expression("TCGA-BRCA")` | (59427, 1231) | 0.37 s | 296 MB |

The gene subset is pushed down as a `gene_id in [...]` parquet filter, so the 5-gene call never builds the 296 MB frame. `attrs` values are read-only dicts (`gene_id`, `gene_type`) that survive slicing and `concat` at no cost.

## Somatic mutations

`mutations` reads the stacked masked MAFs. `genes` becomes a `Hugo_Symbol in [...]` parquet filter, `non_silent=True` keeps the `gdc.NON_SILENT` classes (the maftools nonSyn set), and `columns=None` reads 17 useful columns (`gdc.DEFAULT_MAF_COLUMNS`); pass `"all"` for the 142.

```python
m = gdc.mutations("TCGA-BLCA", genes=["TP53", "RB1"], non_silent=True)
m[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "t_depth", "t_alt_count"]].head()
```
```
        Tumor_Sample_Barcode Hugo_Symbol Variant_Classification  HGVSp_Short  t_depth  t_alt_count
TCGA-GC-A3YS-01A-11D-A23M-08        TP53      Missense_Mutation      p.R158H       45           18
TCGA-DK-A3IS-01A-21D-A21A-08        TP53        Frame_Shift_Ins p.S240Kfs*24       96           57
TCGA-FD-A43P-01A-31D-A23U-08         RB1      Nonsense_Mutation      p.R320*       61           22
TCGA-FD-A43P-01A-31D-A23U-08        TP53      Missense_Mutation      p.D228N       28           10
TCGA-FD-A6TC-01A-21D-A339-08         RB1      Nonsense_Mutation      p.E315*       36           16
```

`Tumor_Sample_Barcode` is the tumor aliquot; the MAF's `aliquot_submitter_id` holds normal and tumor pipe-joined, as in the manifest. The whole BLCA MAF is 117,053 rows from 414 tumor aliquots (415 MAF files; one is empty).

```python
gdc.mutation_frequency("TCGA-BLCA").head(8).round(3)
```
```
Hugo_Symbol
TP53      0.488
TTN       0.425
KMT2D     0.263
ARID1A    0.251
KDM6A     0.251
MUC16     0.249
PIK3CA    0.205
SYNE1     0.191
```

The denominator (`attrs["n_samples"]`, 414) is every barcode in the MAF, before the non-silent filter. With `non_silent=False` TTN comes first (0.522) because of its silent load.

### Joining tables through the case id

For TCGA and TARGET, `case_for_barcode` truncates a barcode to the case (12 or 16 characters). Other programs (CPTAC, HCMI, MMRF, ...) have no positional structure: join through `aliquots()` (`aliquot_submitter_id` to `case_submitter_id`) or the `*_samples` tables instead.

```python
import numpy as np, pandas as pd
tp53 = set(gdc.mutations("TCGA-BLCA", genes="TP53", non_silent=True, columns=["Tumor_Sample_Barcode"])
           ["Tumor_Sample_Barcode"].map(gdc.case_for_barcode))
mdm2 = gdc.expression("TCGA-BLCA", genes="MDM2", sample_type="Primary Tumor").loc["MDM2"]
mut = mdm2.index.map(gdc.case_for_barcode).isin(tp53)
pd.Series({"TP53 mutant": np.median(mdm2[mut]), "TP53 wild-type": np.median(mdm2[~mut])}).round(1)
```
```
TP53 mutant       11.7
TP53 wild-type    22.7
```

198 mutant and 214 wild-type primary tumors; the lower MDM2 in TP53-mutant tumors is the expected loss of p53-driven MDM2 transcription.

## Copy number

Four gene-level workflows are separate pipelines and are not interchangeable; `copy_number` defaults to ascat3 when present and records the choice in `attrs["workflow"]`.

```
>>> gdc.copy_number_workflows("TCGA-BLCA")
['absolute_liftover', 'ascat2', 'ascat3', 'ascatngs']
>>> gdc.copy_number_workflows("TCGA-BLCA", "segments"), gdc.copy_number_workflows("TCGA-BLCA", "allele_specific")
(['dnacopy', 'gatk4_cnv'], ['ascat2', 'ascat3', 'ascatngs'])
```

```python
cn = gdc.copy_number("TCGA-BLCA", genes=["CDKN2A", "ERBB2", "TP53"])
cn.iloc[:, :4]
```
```
           TCGA-DK-A6B2-01A-11D-A30D-01  TCGA-GV-A3QG-01A-11D-A21Y-01  TCGA-C4-A0EZ-01A-21D-A10T-01  TCGA-DK-A3IU-01A-11D-A20B-01
gene_name
CDKN2A                              2.0                           2.0                           3.0                           2.0
TP53                                3.0                           2.0                           3.0                           2.0
ERBB2                               2.0                           2.0                           5.0                           2.0
```
```
>>> cn.shape, cn.attrs["workflow"], dict(cn.attrs["gene_location"])
((3, 392), 'ascat3', {'CDKN2A': ('chr9', 21967752, 21995301), 'TP53': ('chr17', 7661779, 7687538), 'ERBB2': ('chr17', 39687914, 39730426)})
>>> cn.loc["CDKN2A"].value_counts().sort_index()
0.0    112
1.0     32
2.0    132
3.0     51
4.0     37
...
```

112 of 392 BLCA tumors have homozygous CDKN2A deletion by ASCAT3. Values are total copy number (float, NaN when not called).

The sibling samples table, `copy_number_samples`, maps each column to its own aliquot: `sample_type`, `sample_submitter_id` and `case_submitter_id` describe the column. Every column is the tumor aliquot (389 of 392 end in `-01A`, the other three are `-01B`/`-01C` vials), so the table says:

```
>>> gdc.copy_number_samples("TCGA-BLCA")["sample_type"].value_counts().to_dict()
{'Primary Tumor': 392}
```

The paired workflows (ASCAT2, ASCAT3, AscatNGS) keep the tumor/normal pair the file was called on pipe-joined in `aliquot_submitter_id`, in no fixed order (in BLCA ASCAT3 the tumor is first in 186 rows and second in 206), the same convention as the MAF. Tables built before 2026-09-13 took `sample_type` from the first aliquot of that pair, which is why 195 BLCA rows used to say Blood Derived Normal. Every matrix column is a single tumor aliquot.

Long segment tables come from `segments(project, kind, workflow)`, with `kind` in `"segments"` (dnacopy 588,674 rows or gatk4_cnv 2,386,633 rows for BLCA), `"masked"` (266,389) and `"allele_specific"`:

```python
gdc.segments("TCGA-BLCA", "allele_specific").head(3)
```
```
                             file_id         aliquot_submitter_id                          GDC_Aliquot Chromosome    Start      End  Copy_Number  Major_Copy_Number  Minor_Copy_Number
3de11494-5753-4b5d-8fe6-565eb3337a70 TCGA-FD-A3NA-01A-11D-A219-01 93129d0d-3817-4aa2-b886-51726cd264ff       chr1    61735 24072270            4                  3                  1
3de11494-5753-4b5d-8fe6-565eb3337a70 TCGA-FD-A3NA-01A-11D-A219-01 93129d0d-3817-4aa2-b886-51726cd264ff       chr1 24074709 24291499            6                  5                  1
3de11494-5753-4b5d-8fe6-565eb3337a70 TCGA-FD-A3NA-01A-11D-A219-01 93129d0d-3817-4aa2-b886-51726cd264ff       chr1 24294311 92802098            4                  3                  1
```

## miRNA and RPPA

```python
gdc.mirna("TCGA-BLCA").iloc[:3, :3].round(1)          # kind="rpm" (default) or "counts"
```
```
              TCGA-DK-A3IU-01A-11R-A20E-13  TCGA-CU-A0YR-01A-12R-A10V-13  TCGA-FD-A43X-01A-11R-A23X-13
miRNA_ID
hsa-let-7a-1                   4545.899902                  10222.900391                   5481.200195
hsa-let-7a-2                   4555.700195                  10233.200195                   5391.100098
hsa-let-7a-3                   4573.600098                  10371.000000                   5459.899902
```

1,881 miRNAs x 437 aliquots; `mirna_samples` maps the columns, `mirna_isoforms` is the 2.1 M-row long table.

```python
r = gdc.rppa("TCGA-BLCA")
r.iloc[:3, :3].round(3)
```
```
                TCGA-GV-A3QG-01A-21  TCGA-DK-A6B2-01A-21  TCGA-HQ-A2OF-01A-21
peptide_target
1433BETA                      0.054                0.053                0.090
1433EPSILON                   0.120               -0.069               -0.121
1433ZETA                      0.424                0.436                0.453
```

487 antibodies x 343 samples; `r.attrs["agid"]["1433BETA"]` is `'AGID00100'`, and `rppa_antibodies` / `rppa_samples` hold the catalog and column mapping. RPPA columns are 19-character portion ids, so `case_for_barcode` works on them too.

## Methylation

Methylation matrices are large and are built only on request, so most projects have none:

```
>>> gdc.methylation_platforms("TCGA-BLCA")
[]
>>> gdc.methylation("TCGA-BLCA")
NotMirroredError: TCGA-BLCA has no methylation tables (methylation is built only on request)
```

A TCGA-CHOL build with `methylation/betas_human_methylation_450.parquet` (486,427 probes x 45 aliquots, ten 50,000-probe row groups) shows the read paths. `probes` becomes a `probe_id in [...]` row filter: pyarrow skips row groups whose statistics exclude those ids and materialises only the matching rows. `columns` (aliquot barcodes) and `sample_type` become a column projection.

```python
mm = gdc.methylation("TCGA-CHOL", probes=["cg00000029", "cg00000165", "cg00000236"],
                     sample_type="Primary Tumor", root=".../gdcbuild_meth/TCGA-CHOL")
mm.iloc[:, :3].round(3)
```
```
            TCGA-W5-AA2X-01A-11D-A418-05  TCGA-W5-AA2Q-01A-11D-A418-05  TCGA-ZU-A8S4-01A-11D-A418-05
probe_id
cg00000029                         0.177                         0.175                         0.546
cg00000165                         0.080                         0.155                         0.271
cg00000236                         0.891                         0.876                         0.891
```
```
>>> mm.shape, mm.attrs["platform"]      # 0.009 s
((3, 36), 'human_methylation_450')
```

Memory is bounded by rows kept x columns kept x 4 bytes plus one row group of the projected columns. With no `probes`, `columns` or `sample_type` the whole platform matrix is read: 116 MB in pandas for CHOL (0.12 s locally); a 1,200-aliquot 450k project would be about 2.3 GB, so subset on large projects. `methylation_samples(project, platform)` maps the columns.

## Raw files

`manifest(project, data_type=None)` lists every open-access file with its source bucket, key, md5 and HTTPS URL. `fetch_raw` downloads one of them with a signed default boto3 client (newer GDC objects are SSE-KMS and refuse anonymous reads; any AWS principal works) and falls back to `gdc_download_url` over HTTPS on AccessDenied/403/404 or when no credentials are configured. The bytes are md5-checked against the manifest.

```python
row = gdc.manifest("TCGA-BLCA", data_type="Masked Somatic Mutation").sort_values("file_size").iloc[0]
row[["file_name", "file_size", "md5sum", "s3_bucket", "gdc_download_url"]]
```
```
file_name           8dd5c052-b294-4537-8e13-30e50c20a767.wxs.aliqu...
file_size                                                        1143
md5sum                               5d61dcadc7689ab108dd4c0528a6ccdf
s3_bucket                                                 tcga-2-open
gdc_download_url    https://api.gdc.cancer.gov/data/c8459e65-540a-...
```
```
>>> data = gdc.fetch_raw(row)                      # 0.32 s, bytes
>>> len(data), gzip.decompress(data).decode().splitlines()[0]
(1143, '#version gdc-1.0.0')
>>> gdc.fetch_raw("TCGA-BLCA", row["file_id"], dest="/tmp")   # writes /tmp/<file_name>, returns the path
```

## Error messages

Bad arguments list what is valid:

```
>>> gdc.expression("TCGA-BLCA", sample_type="Metastatic")
ValueError: no RNA-Seq aliquots with sample_type ['Metastatic']; present: {'Primary Tumor': 412, 'Solid Tissue Normal': 19}
>>> gdc.copy_number("TCGA-BLCA", workflow="nope")
ValueError: unknown gene-level copy number workflow 'nope' for TCGA-BLCA; available: ['absolute_liftover', 'ascat2', 'ascat3', 'ascatngs']
>>> gdc.manifest("TCGA-BLCA", data_type="nope")
ValueError: no files of data_type 'nope' in TCGA-BLCA; valid: Allele-specific Copy Number Segment, Biospecimen Supplement, Clinical Supplement, ...
```

## Function list

| function | returns |
|---|---|
| `projects()` | one row per mirrored project |
| `report(p)`, `tables(p)`, `readme(p)` | build report dict, table sizes, README text |
| `manifest(p, data_type=None)` | every open file with S3 location and md5 |
| `clinical(p, primary_only=True)` | cases joined to the primary diagnosis, with `os_time_days`, `os_event`, `age_at_diagnosis_years`, `stage` |
| `diagnoses`, `treatments`, `exposures`, `follow_ups`, `samples`, `aliquots` | the clinical tables as stored |
| `expression(p, genes, sample_type, kind, protein_coding, drop_par_y)` | genes x aliquots, `attrs["gene_id"]` |
| `expression_samples(p)`, `star_qc(p)` | column mapping, STAR summary rows |
| `mutations(p, genes, non_silent, columns)` | MAF rows |
| `mutation_frequency(p, non_silent=True)` | fraction of tumor barcodes mutated per gene |
| `copy_number(p, workflow, genes)` | genes x tumor aliquots, `attrs["workflow"]`, `attrs["gene_location"]` |
| `copy_number_workflows(p, kind)`, `copy_number_samples(p, workflow)` | available workflows, column mapping |
| `segments(p, kind, workflow)` | long segment tables |
| `mirna(p, kind, sample_type)`, `mirna_samples`, `mirna_isoforms` | miRNA matrices and tables |
| `rppa(p, sample_type)`, `rppa_samples`, `rppa_antibodies` | RPPA matrix (`attrs["agid"]`) and tables |
| `methylation(p, platform, probes, columns, sample_type)` | probes x aliquots betas |
| `methylation_platforms(p)`, `methylation_samples(p, platform)` | platforms, column mapping |
| `case_for_barcode(barcode)` | case id for TCGA/TARGET barcodes |
| `fetch_raw(row_or_project, file_id, dest)` | one raw GDC file, md5-verified |
