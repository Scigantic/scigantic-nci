# scigantic-nci documentation

- [gdc.md](gdc.md): `scigantic_nci.gdc`, GDC projects as pandas tables (clinical, expression, mutations, copy number, miRNA, RPPA, methylation, raw files). Every snippet was run against the mirror and its output pasted.
- [idc.md](idc.md): `scigantic_nci.idc`, Imaging Data Commons collections (series, studies, patients, slides, analysis results, clinical tables, sample DICOM, pulls from the raw buckets).
- [cli.md](cli.md): the `scigantic-nci` console script, one example per subcommand with its output.
- [table_schemas.md](table_schemas.md): the exact parquet schemas of one GDC project (TCGA-BLCA) and one IDC collection (tcga_lihc).

The README covers installation, why the mirrors exist, provenance and citation, the resolution order, and the limits.
