#!/usr/bin/env python3
"""Rewrite every copy_number/gene_level_<wf>_samples.parquet (and any other *_samples.parquet) of built projects so
sample_type/tissue_type describe the aliquot the matrix column is named after. Earlier builds described the first
aliquot of the tumor/normal pair (often the blood normal). Usage: fix_cn_samples.py <manifest-dir> <mirror-dir> [PROJECT...]"""
import json, os, sys, glob
import pandas as pd, pyarrow as pa, pyarrow.parquet as pq
sys.path.insert(0, os.path.dirname(__file__))
from build_gdc_project_tables import load_manifest, write_parquet
man_dir, mirror = sys.argv[1], sys.argv[2]
_, _, aliquot_sample = load_manifest(man_dir)
projects = sys.argv[3:] or sorted(p for p in os.listdir(mirror) if os.path.exists(os.path.join(mirror, p, "BUILD_REPORT.json")))
for p in projects:
    for f in sorted(glob.glob(os.path.join(mirror, p, "*", "*_samples.parquet"))):
        df = pd.read_parquet(f); changed = 0
        for i, row in df.iterrows():
            a = str(row["column"]).split("__")[0]
            meta = aliquot_sample.get(a)
            if not meta:
                continue
            for k in ("sample_submitter_id", "sample_type", "tissue_type", "tumor_descriptor", "case_submitter_id"):
                if meta.get(k) is not None and row.get(k) != meta[k]:
                    df.at[i, k] = meta[k]; changed += 1
        if changed:
            write_parquet(df, f)
        print(f"{p:26s} {os.path.relpath(f, os.path.join(mirror, p)):58s} {len(df):5d} rows, {changed} fields corrected")
