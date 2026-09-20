# scigantic-nci builders

The scripts that assemble the analysis-ready NCI companion tables from the public
Genomic Data Commons (GDC) and Imaging Data Commons (IDC) tiers, and publish them to
anonymous S3. They are the producers of the tables that the
[`scigantic-nci`](https://pypi.org/project/scigantic-nci/) package reads.

Everything here runs against public interfaces only: the GDC API, the AWS Open Data
mirrors, and idc-index. No account or agreement is needed beyond an AWS profile to
write the output bucket. A commons, a cloud resource, or an institution can run these
to produce the same companions under its own name and bucket.

Licensed MIT No Attribution (see the repository [LICENSE](../LICENSE)): use, modify and redistribute without
restriction or attribution.

## Layout

    gdc/
      dump_gdc_manifest.py          list the GDC open-access files and their AWS bucket map
      build_gdc_project_tables.py   assemble per-project tables (clinical, RNA-Seq, MAF,
                                    copy number, miRNA, RPPA, methylation); MD5-verify
                                    every source file; write README.md + BUILD_REPORT.json
      audit_gdc_project_build.py    independent re-read of the built tables by a second
                                    code path: values, row counts, column ids, case totals
      fix_cn_samples.py             repair copy-number sample labels (tumor vs normal)
      fix_pair_columns.py           re-key tumor/normal pair columns to the resolved aliquot
    idc/
      build_idc_collection_archives.py  per-collection series/study/patient indexes,
                                    clinical tables, and one redistributable sample series

## GDC

    python gdc/dump_gdc_manifest.py --out manifest
    python gdc/build_gdc_project_tables.py --manifest-dir manifest --out gdc_mirror \
        --projects TCGA-BRCA --methylation --upload
    python gdc/audit_gdc_project_build.py TCGA-BRCA --mirror gdc_mirror

Use `--projects all-bucket` for every project served from AWS. Every source file is
checked against the GDC manifest MD5 before use; a mismatch fails the project.
`--verify` re-reads S3 and checks the output against the build report. Rebuilding a
project produces byte-identical tables.

## IDC

    python idc/build_idc_collection_archives.py --out idc_mirror --collections tcga_brca --upload

Use `--collections all` for every collection. Series are read anonymously; only series
whose license permits redistribution are copied as samples.

## Dependencies

    pip install boto3 numpy pandas pyarrow idc-index

Python 3.10 or newer. `idc-index` is needed only for the IDC builder. Writing to S3
needs an AWS profile with access to the destination bucket.
