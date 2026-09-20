#!/usr/bin/env python3
"""
Build per-project "AI-ready" tables from the NCI Genomic Data Commons (GDC)
open-access tier and upload them to s3://scigantic-gdc-open/<PROJECT>/.

Why this exists
---------------
The GDC open buckets on AWS (tcga-2-open, gdc-target-phs000218-2-open, ...)
store every file under its own UUID prefix: 284,000 flat folders for TCGA
alone, no project or data-type hierarchy, one file per sample. A FUSE mount of
that is unusable for browsing and every analysis starts by re-assembling the
same matrices. This script assembles them once, per project, and records
exactly how.

What it builds (per project)
----------------------------
  MANIFEST.parquet             every open file of the project: ids, type,
                               workflow, size, md5, case/sample/aliquot,
                               s3 bucket + key, GDC download URL
  clinical/                    cases, diagnoses, treatments, exposures,
                               follow_ups, samples, aliquots (GDC cases API)
  rna_seq/                     STAR-Counts: genes, samples, counts_unstranded
                               (int32), tpm_unstranded (float32), star_qc
  somatic_mutations/           all masked MAFs concatenated (long)
  copy_number/                 gene-level copy number per workflow (wide),
                               segments per workflow (long)
  mirna/                       miRNA read_count + RPM (wide), isoforms (long)
  protein/                     RPPA protein_expression (wide) + antibodies
  methylation/                 beta values per platform (wide float32,
                               50k-probe row groups)   [--methylation]
  README.md, BUILD_REPORT.json

Every downloaded file is md5-verified against the GDC manifest. A mismatch
fails the project. Wide matrices assert that the feature order is identical
across files (STAR, ASCAT gene-level, miRNA, SeSAMe) and re-align by id if
not, recording it in the report.

Usage
-----
  python build_gdc_project_tables.py --manifest-dir <dir> --out <workdir> \
      --projects TCGA-BRCA,TCGA-LUAD [--methylation] [--workers 12] [--upload]
  python build_gdc_project_tables.py ... --projects all-bucket   # 60 projects served from AWS buckets
  python build_gdc_project_tables.py ... --verify --projects TCGA-BRCA  # re-read S3 and check against the report

<manifest-dir> must contain gdc_open_files.ndjson, gdc_assoc.ndjson and
project_bucket_map.json (produced by dump_gdc_manifest.py in this folder).

Read side: signed boto3 (some newer GDC objects are SSE-KMS and refuse
unsigned GETs; the key policy allows any signed principal). Projects with no
AWS bucket are read from https://api.gdc.cancer.gov/data/<file_id>.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import hashlib
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone

import boto3
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from botocore.config import Config

MIRROR_BUCKET = "scigantic-gdc-open"
MIRROR_REGION = "us-east-1"
GDC_API = "https://api.gdc.cancer.gov"

NORMAL_SAMPLE_TYPES = {
    "Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal",
    "Bone Marrow Normal", "Fibroblasts from Bone Marrow Normal",
    "EBV Immortalized Normal", "Lymphoid Normal", "Mononuclear Cells from Bone Marrow Normal",
    "Peripheral Blood Components NOS", "Blood Derived Cancer - Peripheral Blood, Post-treatment",
}

TABLE_GROUPS = ["clinical", "rna", "maf", "cnv", "mirna", "protein", "methylation"]


# --------------------------------------------------------------------------- io

class Fetcher:
    """Downloads GDC files (S3 signed, or HTTPS API), verifying md5."""

    def __init__(self, workers: int):
        self.s3 = boto3.client("s3", region_name="us-east-1",
                               config=Config(max_pool_connections=max(8, workers * 2), retries={"max_attempts": 8}))
        self.bytes_read = 0
        self.api_fallbacks = 0

    def fetch(self, row: dict) -> bytes:
        last = None
        for attempt in range(4):
            try:
                if row["s3_bucket"] and attempt < 2:
                    try:
                        body = self.s3.get_object(Bucket=row["s3_bucket"], Key=row["s3_key"])["Body"].read()
                    except self.s3.exceptions.ClientError as e:
                        if e.response.get("Error", {}).get("Code") not in ("AccessDenied", "NoSuchKey", "403", "404"):
                            raise
                        # a handful of open files are missing or locked in the AWS mirror; the GDC API still serves them
                        self.api_fallbacks += 1
                        req = urllib.request.Request(row["gdc_download_url"], headers={"User-Agent": "scigantic-gdc-mirror/1.0"})
                        with urllib.request.urlopen(req, timeout=600) as r:
                            body = r.read()
                else:
                    req = urllib.request.Request(row["gdc_download_url"], headers={"User-Agent": "scigantic-gdc-mirror/1.0"})
                    with urllib.request.urlopen(req, timeout=600) as r:
                        body = r.read()
                md5 = hashlib.md5(body).hexdigest()
                if md5 != row["md5sum"]:
                    raise ValueError(f"md5 mismatch for {row['file_id']}: got {md5} want {row['md5sum']}")
                self.bytes_read += len(body)
                return body
            except Exception as e:  # noqa: BLE001
                last = e
                time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"failed to fetch {row['file_id']} ({row['file_name']}): {last}")


def text_of(body: bytes, name: str) -> str:
    if name.endswith(".gz"):
        body = gzip.decompress(body)
    return body.decode("utf-8")


# ---------------------------------------------------------------------- manifest

def load_manifest(manifest_dir: str) -> tuple[pd.DataFrame, dict, dict]:
    files = []
    with open(os.path.join(manifest_dir, "gdc_open_files.ndjson")) as fh:
        for line in fh:
            h = json.loads(line)
            cases = h.get("cases") or []
            projects = sorted({c.get("project", {}).get("project_id", "") for c in cases})
            samples = [s for c in cases for s in (c.get("samples") or [])]
            files.append({
                "file_id": h["file_id"], "file_name": h["file_name"],
                "data_category": h.get("data_category"), "data_type": h.get("data_type"),
                "data_format": h.get("data_format"), "experimental_strategy": h.get("experimental_strategy"),
                "file_size": int(h.get("file_size") or 0), "md5sum": h.get("md5sum"),
                "project_id": projects[0] if len(projects) == 1 else "|".join(projects),
                "program": (cases[0].get("project", {}).get("program", {}) or {}).get("name") if cases else None,
                "case_id": "|".join(sorted({c.get("case_id", "") for c in cases})),
                "case_submitter_id": "|".join(sorted({c.get("submitter_id", "") for c in cases})),
                "sample_submitter_id": "|".join(sorted({s.get("submitter_id", "") for s in samples})),
                "sample_type": "|".join(sorted({s.get("sample_type", "") for s in samples})),
            })
    df = pd.DataFrame(files)

    assoc: dict[str, dict] = {}
    aliquot_sample: dict[str, dict] = {}
    with open(os.path.join(manifest_dir, "gdc_assoc.ndjson")) as fh:
        for line in fh:
            h = json.loads(line)
            ents = h.get("associated_entities") or []
            assoc[h["file_id"]] = {
                "aliquots": [e.get("entity_submitter_id") for e in ents],
                "entity_type": "|".join(sorted({e.get("entity_type", "") for e in ents})),
                "workflow_type": (h.get("analysis") or {}).get("workflow_type"),
                "platform": h.get("platform"),
            }
            for c in h.get("cases") or []:
                for s in c.get("samples") or []:
                    for p in s.get("portions") or []:
                        for a in p.get("analytes") or []:
                            for q in a.get("aliquots") or []:
                                aliquot_sample.setdefault(q["submitter_id"], {
                                    "sample_submitter_id": s.get("submitter_id"), "sample_type": s.get("sample_type"),
                                    "tissue_type": s.get("tissue_type"), "tumor_descriptor": s.get("tumor_descriptor"),
                                    "case_submitter_id": c.get("submitter_id"),
                                })
    bucket_map = json.load(open(os.path.join(manifest_dir, "project_bucket_map.json")))
    df["workflow_type"] = df.file_id.map(lambda f: assoc.get(f, {}).get("workflow_type"))
    df["platform"] = df.file_id.map(lambda f: assoc.get(f, {}).get("platform"))
    df["aliquot_submitter_id"] = df.file_id.map(lambda f: "|".join(a for a in assoc.get(f, {}).get("aliquots", []) if a))
    df["s3_bucket"] = df.project_id.map(lambda p: (bucket_map.get(p) or "").replace(" (PARTIAL)", "") or None)
    df["s3_key"] = df.file_id + "/" + df.file_name
    df.loc[df.s3_bucket.isna(), "s3_key"] = None
    df["gdc_download_url"] = GDC_API + "/data/" + df.file_id
    return df, assoc, aliquot_sample


def tumor_aliquot(aliquots: list[str], aliquot_sample: dict) -> str:
    """For paired tumor/normal files pick the tumor aliquot as the column id.

    GDC's sample tissue_type decides first: exactly one Tumor and the rest
    Normal -> the Tumor. The sample_type list is only a fallback, because it
    misses normals such as Granulocytes, Saliva, Slides and FFPE Scrolls (the
    CGCI-BLGSP and HCMI-CMDC pairs). If neither rule picks exactly one aliquot
    the id stays the sorted pipe-joined pair.
    """
    if len(aliquots) == 1:
        return aliquots[0]
    tissue = [aliquot_sample.get(a, {}).get("tissue_type") or "" for a in aliquots]
    if tissue.count("Tumor") == 1 and tissue.count("Normal") == len(aliquots) - 1:
        return aliquots[tissue.index("Tumor")]
    tumor = [a for a, t in zip(aliquots, tissue)
             if t != "Normal" and (aliquot_sample.get(a, {}).get("sample_type") or "") not in NORMAL_SAMPLE_TYPES]
    if len(tumor) == 1:
        return tumor[0]
    return "|".join(sorted(aliquots))


# ---------------------------------------------------------------------- helpers

def write_parquet(df: pd.DataFrame, path: str, row_group_size: int | None = None) -> dict:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, path, compression="zstd", row_group_size=row_group_size)
    return {"rows": table.num_rows, "columns": table.num_columns, "bytes": os.path.getsize(path)}


def unique_columns(ids: list[str], file_ids: list[str], report: dict, table: str) -> list[str]:
    seen = Counter(ids)
    out, dup = [], []
    for i, f in zip(ids, file_ids):
        if seen[i] > 1:
            out.append(f"{i}__{f[:8]}")
            dup.append(i)
        else:
            out.append(i)
    if dup:
        report.setdefault("duplicate_aliquots", {})[table] = sorted(set(dup))
    return out


def sample_table(rows: list[dict], colnames: list[str], aliquot_sample: dict) -> pd.DataFrame:
    recs = []
    for r, c in zip(rows, colnames):
        # describe the aliquot the column is named after (the tumor of a tumor/normal pair), not the first of the pair
        a = c.split("__")[0] if c else None
        if a not in aliquot_sample and r["aliquot_submitter_id"]:
            a = r["aliquot_submitter_id"].split("|")[0]
        meta = aliquot_sample.get(a, {})
        recs.append({
            "column": c, "aliquot_submitter_id": r["aliquot_submitter_id"], "file_id": r["file_id"], "file_name": r["file_name"],
            "case_submitter_id": r["case_submitter_id"], "sample_submitter_id": meta.get("sample_submitter_id") or r["sample_submitter_id"],
            "sample_type": meta.get("sample_type") or r["sample_type"], "tissue_type": meta.get("tissue_type"),
            "tumor_descriptor": meta.get("tumor_descriptor"), "workflow_type": r["workflow_type"], "platform": r["platform"],
        })
    return pd.DataFrame(recs)


def build_wide(rows: list[dict], fetcher: Fetcher, parse, feature_cols: list[str], value_specs: list[tuple[str, np.dtype]],
               workers: int, report: dict, table: str, aliquot_sample: dict, colname_fn=None):
    """Generic wide-matrix assembler. `parse(text) -> (features_df, {value_name: np.ndarray})`.

    Matrices are preallocated column-major (one column per file) so a 486k x
    1,200 methylation matrix costs 2.4 GB, not several copies of it. Verifies
    the feature order matches the first file and re-aligns by the first feature
    column when it does not (an int matrix is promoted to float32 if that leaves
    holes). Returns (features_df, {name: 2D array}, colnames, samples_df).
    """
    ref_feat = None
    ref_index = None
    mats: dict[str, np.ndarray] = {}
    realigned = 0
    ids, fids = [], []
    n = len(rows)

    def work(r):
        return r, parse(text_of(fetcher.fetch(r), r["file_name"]))

    i = 0
    with cf.ThreadPoolExecutor(workers) as ex:
        for batch_start in range(0, n, workers * 4):
            batch = rows[batch_start: batch_start + workers * 4]
            for r, (feat, vals) in ex.map(work, batch):
                if ref_feat is None:
                    ref_feat = feat.reset_index(drop=True)
                    ref_index = pd.Index(ref_feat[feature_cols[0]])
                    for name, dt in value_specs:
                        mats[name] = np.empty((len(ref_feat), n), dtype=dt, order="F")
                if len(feat) != len(ref_feat) or not (feat[feature_cols[0]].values == ref_feat[feature_cols[0]].values).all():
                    realigned += 1
                    pos = ref_index.get_indexer(feat[feature_cols[0]])
                    ok = pos >= 0
                    for name, dt in value_specs:
                        arr = np.full(len(ref_feat), np.nan, dtype=np.float64)
                        arr[pos[ok]] = np.asarray(vals[name])[ok]
                        if np.issubdtype(mats[name].dtype, np.integer) and np.isnan(arr).any():
                            mats[name] = mats[name].astype(np.float32, order="F")
                            report.setdefault("promoted_to_float", []).append(f"{table}/{name}")
                        vals[name] = arr
                for name, dt in value_specs:
                    mats[name][:, i] = vals[name]
                ids.append(colname_fn(r) if colname_fn else tumor_aliquot(r["aliquot_submitter_id"].split("|"), aliquot_sample))
                fids.append(r["file_id"])
                i += 1
            print(f"    {table}: {i}/{n} files", flush=True)
    cols = unique_columns(ids, fids, report, table)
    if realigned:
        report.setdefault("realigned_files", {})[table] = realigned
    return ref_feat, mats, cols, sample_table(rows, cols, aliquot_sample)


def wide_table(features: pd.DataFrame, matrix: np.ndarray, cols: list[str]) -> pa.Table:
    """Arrow table [feature columns..., one column per file] without a pandas copy."""
    arrays = [pa.array(features[c].to_numpy()) if features[c].dtype != object else pa.array(features[c].tolist(), pa.string()) for c in features.columns]
    names = list(features.columns)
    for j, c in enumerate(cols):
        arrays.append(pa.array(matrix[:, j]))
        names.append(c)
    return pa.Table.from_arrays(arrays, names=names)


def write_table(table: pa.Table, path: str, row_group_size: int | None = None) -> dict:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pq.write_table(table, path, compression="zstd", row_group_size=row_group_size)
    return {"rows": table.num_rows, "columns": table.num_columns, "bytes": os.path.getsize(path)}


# ----------------------------------------------------------------------- tables

def build_clinical(project: str, out: str, report: dict) -> None:
    expand = "demographic,diagnoses,diagnoses.treatments,exposures,follow_ups,samples,samples.portions.analytes.aliquots"
    filt = json.dumps({"op": "=", "content": {"field": "project.project_id", "value": project}})
    hits, frm, size = [], 0, 500
    while True:
        q = urllib.parse.urlencode({"filters": filt, "expand": expand, "size": size, "from": frm, "format": "json"})
        with urllib.request.urlopen(f"{GDC_API}/cases?{q}", timeout=600) as r:
            d = json.load(r)
        page = d["data"]["hits"]
        hits.extend(page)
        frm += size
        if frm >= d["data"]["pagination"]["total"] or not page:
            break

    def scalars(d: dict, skip=()):
        o = {}
        for k, v in d.items():
            if k in skip:
                continue
            if isinstance(v, list):
                if all(not isinstance(x, dict) for x in v):
                    o[k] = "|".join(str(x) for x in v)
            elif isinstance(v, dict):
                for k2, v2 in v.items():
                    if not isinstance(v2, (list, dict)):
                        o[f"{k}_{k2}"] = v2
            else:
                o[k] = v
        return o

    cases, dx, tx, ex, fu, samples, aliquots = [], [], [], [], [], [], []
    for h in hits:
        cid, sub = h["case_id"], h["submitter_id"]
        cases.append(scalars(h, skip=("diagnoses", "exposures", "follow_ups", "samples", "aliquot_ids", "analyte_ids", "portion_ids",
                                     "sample_ids", "slide_ids", "diagnosis_ids", "submitter_aliquot_ids", "submitter_analyte_ids",
                                     "submitter_portion_ids", "submitter_sample_ids", "submitter_slide_ids", "submitter_diagnosis_ids")))
        for d in h.get("diagnoses") or []:
            dx.append({"case_id": cid, "case_submitter_id": sub, **scalars(d, skip=("treatments",))})
            for t in d.get("treatments") or []:
                tx.append({"case_id": cid, "case_submitter_id": sub, "diagnosis_id": d.get("diagnosis_id"), **scalars(t)})
        for e in h.get("exposures") or []:
            ex.append({"case_id": cid, "case_submitter_id": sub, **scalars(e)})
        for f in h.get("follow_ups") or []:
            fu.append({"case_id": cid, "case_submitter_id": sub, **scalars(f)})
        for s in h.get("samples") or []:
            samples.append({"case_id": cid, "case_submitter_id": sub, **scalars(s, skip=("portions",))})
            for p in s.get("portions") or []:
                for a in p.get("analytes") or []:
                    for q in a.get("aliquots") or []:
                        aliquots.append({"case_submitter_id": sub, "sample_submitter_id": s.get("submitter_id"),
                                         "sample_type": s.get("sample_type"), "tissue_type": s.get("tissue_type"),
                                         "tumor_descriptor": s.get("tumor_descriptor"), "aliquot_id": q.get("aliquot_id"),
                                         "aliquot_submitter_id": q.get("submitter_id")})
    tables = {"cases": cases, "diagnoses": dx, "treatments": tx, "exposures": ex, "follow_ups": fu, "samples": samples, "aliquots": aliquots}
    for name, recs in tables.items():
        df = pd.DataFrame(recs)
        for c in list(df.columns):  # drop pure-null columns, keep the rest as parsed
            if df[c].isna().all():
                df = df.drop(columns=c)
                continue
            if df[c].dtype == object:  # GDC mixes int and str in a few fields (year_of_diagnosis in TARGET-ALL-P2); make each column one type
                vals = df[c].dropna()
                kinds = {type(v) for v in vals}
                if len(kinds) > 1 or kinds == {int} or kinds == {float}:
                    num = pd.to_numeric(vals, errors="coerce")
                    df[c] = num if num.notna().all() else df[c].map(lambda v: None if v is None or (isinstance(v, float) and pd.isna(v)) else str(v))
        if not df.empty:
            report["tables"][f"clinical/{name}"] = write_parquet(df, os.path.join(out, "clinical", f"{name}.parquet"))
    report["clinical_cases"] = len(hits)


def parse_star(text: str):
    df = pd.read_csv(io.StringIO(text), sep="\t", comment="#", dtype={"gene_id": str, "gene_name": str, "gene_type": str})
    qc = df.iloc[:4]
    body = df.iloc[4:]
    feat = body[["gene_id", "gene_name", "gene_type"]]
    return feat, {"unstranded": body["unstranded"].to_numpy(np.int64), "tpm_unstranded": body["tpm_unstranded"].to_numpy(np.float64),
                  "_qc": qc[["gene_id", "unstranded", "stranded_first", "stranded_second"]]}


def build_rna(rows, fetcher, out, report, workers, aliquot_sample):
    rows = [r for r in rows if r["workflow_type"] == "STAR - Counts"]
    if not rows:
        return
    qc_rows = []

    def parse(text):
        feat, vals = parse_star(text)
        qc_rows.append(vals.pop("_qc"))
        return feat, vals

    feat, mats, cols, samples = build_wide(rows, fetcher, parse, ["gene_id"], [("unstranded", np.int32), ("tpm_unstranded", np.float32)],
                                           workers, report, "rna_seq", aliquot_sample)
    report["tables"]["rna_seq/genes"] = write_parquet(feat, os.path.join(out, "rna_seq", "genes.parquet"))
    report["tables"]["rna_seq/samples"] = write_parquet(samples, os.path.join(out, "rna_seq", "samples.parquet"))
    report["tables"]["rna_seq/counts_unstranded"] = write_table(wide_table(feat[["gene_id"]], mats["unstranded"], cols), os.path.join(out, "rna_seq", "counts_unstranded.parquet"))
    report["tables"]["rna_seq/tpm_unstranded"] = write_table(wide_table(feat[["gene_id"]], mats["tpm_unstranded"], cols), os.path.join(out, "rna_seq", "tpm_unstranded.parquet"))
    qc = []
    for c, q in zip(cols, qc_rows):
        rec = {"column": c}
        for _, r in q.iterrows():
            rec[f"{r['gene_id']}_unstranded"] = int(r["unstranded"])
        qc.append(rec)
    report["tables"]["rna_seq/star_qc"] = write_parquet(pd.DataFrame(qc), os.path.join(out, "rna_seq", "star_qc.parquet"))


MAF_INT_COLS = ["Start_Position", "End_Position", "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count", "Entrez_Gene_Id"]


def build_maf(rows, fetcher, out, report, workers):
    if not rows:
        return
    parts = []

    def work(r):
        text = text_of(fetcher.fetch(r), r["file_name"])
        df = pd.read_csv(io.StringIO(text), sep="\t", comment="#", dtype=str, low_memory=False, keep_default_na=False, na_values=[""])
        df.insert(0, "file_id", r["file_id"])
        df.insert(1, "aliquot_submitter_id", r["aliquot_submitter_id"])
        return df

    with cf.ThreadPoolExecutor(workers) as ex:
        for i, df in enumerate(ex.map(work, rows), 1):
            parts.append(df)
            if i % 100 == 0 or i == len(rows):
                print(f"    somatic_mutations: {i}/{len(rows)} files", flush=True)
    maf = pd.concat(parts, ignore_index=True)
    for c in MAF_INT_COLS:
        if c in maf.columns:
            maf[c] = pd.to_numeric(maf[c], errors="coerce").astype("Int64")
    report["tables"]["somatic_mutations/masked_somatic_mutations"] = write_parquet(maf, os.path.join(out, "somatic_mutations", "masked_somatic_mutations.parquet"), row_group_size=200_000)
    report["maf_files"] = len(rows)


def parse_gene_cn(text: str):
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype={"gene_id": str, "gene_name": str, "chromosome": str})
    return df[["gene_id", "gene_name", "chromosome", "start", "end"]], {"copy_number": pd.to_numeric(df["copy_number"], errors="coerce").to_numpy(np.float64)}


def build_cnv(rows, fetcher, out, report, workers, aliquot_sample):
    gene_rows = [r for r in rows if r["data_type"] == "Gene Level Copy Number"]
    by_wf = defaultdict(list)
    for r in gene_rows:
        by_wf[r["workflow_type"] or "unknown"].append(r)
    for wf, wrows in by_wf.items():
        slug = wf.lower().replace(" ", "_")
        feat, mats, cols, samples = build_wide(wrows, fetcher, parse_gene_cn, ["gene_id"], [("copy_number", np.float32)], workers, report,
                                               f"copy_number/gene_level_{slug}", aliquot_sample)
        report["tables"][f"copy_number/gene_level_{slug}"] = write_table(wide_table(feat, mats["copy_number"], cols), os.path.join(out, "copy_number", f"gene_level_{slug}.parquet"))
        report["tables"][f"copy_number/gene_level_{slug}_samples"] = write_parquet(samples, os.path.join(out, "copy_number", f"gene_level_{slug}_samples.parquet"))

    seg_types = {"Copy Number Segment": "segments", "Masked Copy Number Segment": "masked_segments", "Allele-specific Copy Number Segment": "allele_specific_segments"}
    for dtype, base in seg_types.items():
        trows = [r for r in rows if r["data_type"] == dtype]
        by_wf = defaultdict(list)
        for r in trows:
            by_wf[r["workflow_type"] or "unknown"].append(r)
        for wf, wrows in by_wf.items():
            slug = wf.lower().replace(" ", "_")

            def work(r):
                df = pd.read_csv(io.StringIO(text_of(fetcher.fetch(r), r["file_name"])), sep="\t", dtype=str)
                df.insert(0, "file_id", r["file_id"])
                df.insert(1, "aliquot_submitter_id", tumor_aliquot(r["aliquot_submitter_id"].split("|"), aliquot_sample))
                return df

            parts = []
            with cf.ThreadPoolExecutor(workers) as ex:
                for i, df in enumerate(ex.map(work, wrows), 1):
                    parts.append(df)
                    if i % 200 == 0 or i == len(wrows):
                        print(f"    copy_number/{base}_{slug}: {i}/{len(wrows)} files", flush=True)
            seg = pd.concat(parts, ignore_index=True)
            for c in seg.columns:
                if c in ("Start", "End", "Num_Probes", "Copy_Number", "Major_Copy_Number", "Minor_Copy_Number"):
                    seg[c] = pd.to_numeric(seg[c], errors="coerce").astype("Int64")
                elif c == "Segment_Mean":
                    seg[c] = pd.to_numeric(seg[c], errors="coerce").astype("float32")
            report["tables"][f"copy_number/{base}_{slug}"] = write_parquet(seg, os.path.join(out, "copy_number", f"{base}_{slug}.parquet"), row_group_size=200_000)


def parse_mirna(text: str):
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype={"miRNA_ID": str, "cross-mapped": str})
    return df[["miRNA_ID"]], {"read_count": df["read_count"].to_numpy(np.int64), "reads_per_million_miRNA_mapped": df["reads_per_million_miRNA_mapped"].to_numpy(np.float64)}


def build_mirna(rows, fetcher, out, report, workers, aliquot_sample):
    mrows = [r for r in rows if r["data_type"] == "miRNA Expression Quantification"]
    if mrows:
        feat, mats, cols, samples = build_wide(mrows, fetcher, parse_mirna, ["miRNA_ID"], [("read_count", np.int32), ("reads_per_million_miRNA_mapped", np.float32)],
                                               workers, report, "mirna", aliquot_sample)
        report["tables"]["mirna/samples"] = write_parquet(samples, os.path.join(out, "mirna", "samples.parquet"))
        report["tables"]["mirna/read_count"] = write_table(wide_table(feat, mats["read_count"], cols), os.path.join(out, "mirna", "read_count.parquet"))
        report["tables"]["mirna/rpm"] = write_table(wide_table(feat, mats["reads_per_million_miRNA_mapped"], cols), os.path.join(out, "mirna", "rpm.parquet"))
    irows = [r for r in rows if r["data_type"] == "Isoform Expression Quantification"]
    if irows:
        def work(r):
            df = pd.read_csv(io.StringIO(text_of(fetcher.fetch(r), r["file_name"])), sep="\t", dtype=str)
            df.insert(0, "file_id", r["file_id"])
            df.insert(1, "aliquot_submitter_id", r["aliquot_submitter_id"])
            return df
        parts = []
        with cf.ThreadPoolExecutor(workers) as ex:
            for i, df in enumerate(ex.map(work, irows), 1):
                parts.append(df)
                if i % 200 == 0 or i == len(irows):
                    print(f"    mirna/isoforms: {i}/{len(irows)} files", flush=True)
        iso = pd.concat(parts, ignore_index=True)
        iso["read_count"] = pd.to_numeric(iso["read_count"], errors="coerce").astype("Int64")
        iso["reads_per_million_miRNA_mapped"] = pd.to_numeric(iso["reads_per_million_miRNA_mapped"], errors="coerce").astype("float32")
        report["tables"]["mirna/isoforms"] = write_parquet(iso, os.path.join(out, "mirna", "isoforms.parquet"), row_group_size=500_000)


def parse_rppa(text: str):
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype={"AGID": str, "lab_id": str, "catalog_number": str, "set_id": str, "peptide_target": str})
    return df[["AGID", "peptide_target", "lab_id", "catalog_number", "set_id"]], {"protein_expression": pd.to_numeric(df["protein_expression"], errors="coerce").to_numpy(np.float64)}


def build_protein(rows, fetcher, out, report, workers, aliquot_sample):
    prows = [r for r in rows if r["data_type"] == "Protein Expression Quantification"]
    if not prows:
        return
    feat, mats, cols, samples = build_wide(prows, fetcher, parse_rppa, ["AGID"], [("protein_expression", np.float32)], workers, report, "protein", aliquot_sample,
                                           colname_fn=lambda r: r["aliquot_submitter_id"] or r["sample_submitter_id"])
    report["tables"]["protein/antibodies"] = write_parquet(feat, os.path.join(out, "protein", "antibodies.parquet"))
    report["tables"]["protein/samples"] = write_parquet(samples, os.path.join(out, "protein", "samples.parquet"))
    report["tables"]["protein/rppa_protein_expression"] = write_table(wide_table(feat[["AGID", "peptide_target"]], mats["protein_expression"], cols), os.path.join(out, "protein", "rppa_protein_expression.parquet"))


def parse_betas(text: str):
    df = pd.read_csv(io.StringIO(text), sep="\t", header=None, names=["probe_id", "beta"], dtype={"probe_id": str}, na_values=["NA"])
    return df[["probe_id"]], {"beta": df["beta"].to_numpy(np.float64)}


def build_methylation(rows, fetcher, out, report, workers, aliquot_sample):
    mrows = [r for r in rows if r["data_type"] == "Methylation Beta Value"]
    by_plat = defaultdict(list)
    for r in mrows:
        by_plat[r["platform"] or "unknown"].append(r)
    for plat, prows in by_plat.items():
        slug = plat.lower().replace("illumina ", "").replace(" ", "_")
        feat, mats, cols, samples = build_wide(prows, fetcher, parse_betas, ["probe_id"], [("beta", np.float32)], workers, report, f"methylation/{slug}", aliquot_sample)
        report["tables"][f"methylation/betas_{slug}_samples"] = write_parquet(samples, os.path.join(out, "methylation", f"betas_{slug}_samples.parquet"))
        report["tables"][f"methylation/betas_{slug}"] = write_table(wide_table(feat, mats["beta"], cols), os.path.join(out, "methylation", f"betas_{slug}.parquet"), row_group_size=50_000)


# ------------------------------------------------------------------------ readme

def render_readme(project: str, man: pd.DataFrame, report: dict, release: str) -> str:
    by_type = man.groupby("data_type").agg(files=("file_id", "count"), GB=("file_size", lambda s: round(s.sum() / 1e9, 2))).sort_values("files", ascending=False)
    lines = [f"# {project}: GDC open-access tables", "",
             f"Built by Scigantic from the NCI Genomic Data Commons open-access tier, {release}, on {report['built_at'][:10]}.",
             f"Source objects: {report['files_read']} files, {report['bytes_read'] / 1e9:.1f} GB read and md5-verified against the GDC file manifest.",
             "", "## What is here", "",
             "| path | rows | columns |", "|---|---|---|"]
    for name, t in sorted(report["tables"].items()):
        lines.append(f"| {name}.parquet | {t['rows']} | {t['columns']} |")
    lines += ["", "## How the tables were made", "",
              "- `MANIFEST.parquet` lists every open-access file of the project with its GDC file_id, the public S3 bucket and key it lives in, and the aliquot barcode(s) it belongs to.",
              "- `clinical/` is the GDC cases API (`expand=demographic,diagnoses,treatments,exposures,follow_ups,samples`) flattened one table per entity. Lists are pipe-joined.",
              "- `rna_seq/` stacks the STAR-Counts `augmented_star_gene_counts.tsv` files column-wise: `counts_unstranded` (raw int) and `tpm_unstranded` (float32), one column per aliquot barcode, gene order as in GENCODE v36. Stranded counts and FPKM are not included; they are in the per-file TSVs listed in the manifest. The four `N_*` STAR summary rows are in `star_qc`.",
              "- `somatic_mutations/` concatenates the aliquot-level masked MAFs (`Aliquot Ensemble Somatic Variant Merging and Masking`). Every MAF column is kept; `file_id` and `aliquot_submitter_id` are prepended.",
              "- `copy_number/` gives one wide `gene_level_<workflow>` matrix per calling workflow (ASCAT2, ASCAT3, ABSOLUTE LiftOver, AscatNGS are different pipelines and are NOT interchangeable) plus long segment tables per workflow. Column ids are the tumor aliquot of each tumor/normal pair.",
              "- `mirna/` stacks BCGSC miRNA profiling: `read_count` and `rpm` wide, isoforms long.",
              "- `protein/` stacks RPPA `protein_expression` by antibody (AGID). Columns are the RPPA portion/aliquot ids.",
              "- `methylation/` (when present) stacks SeSAMe beta values per array platform (27k, 450k, EPIC, EPIC v2 are different probe sets), float32, 50,000-probe row groups so a subset of probes can be read without loading the whole matrix.",
              "", "## Sample ids", "",
              "Matrix columns are GDC aliquot submitter ids (for TCGA the 28-character barcode, e.g. `TCGA-A8-A09E-01A-11R-A10I-07`). `rna_seq/samples.parquet` and the other `*samples.parquet` tables map each column to case, sample, sample_type and file_id. A column suffixed `__xxxxxxxx` is a second file for the same aliquot; the suffix is the first 8 characters of its file_id.",
              "", "## Source files by type", "", "| data_type | files | GB |", "|---|---|---|"]
    for dt, r in by_type.iterrows():
        lines.append(f"| {dt} | {int(r['files'])} | {r['GB']} |")
    lines += ["", "## Attribution and terms", "",
              "GDC open-access data carry no access restrictions. Cite the program that generated the data (TCGA, TARGET, CPTAC, ...) and the GDC (Grossman et al., NEJM 2016, https://doi.org/10.1056/NEJMp1607591). TCGA and TARGET publication guidelines: https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing.",
              "These tables are a derived re-arrangement of GDC files; values are unchanged from the source files. The build report (`BUILD_REPORT.json`) records counts, bytes, md5 checks, duplicates and any re-alignment.",
              ""]
    return "\n".join(lines)


# ------------------------------------------------------------------------ upload

def upload_dir(local: str, project: str) -> dict:
    s3 = boto3.client("s3", region_name=MIRROR_REGION, config=Config(retries={"max_attempts": 10, "mode": "adaptive"}, max_pool_connections=32))
    paths = [os.path.join(root, f) for root, _, files in os.walk(local) for f in files]
    paths.sort(key=lambda p: p.endswith("BUILD_REPORT.json"))  # report last: its presence means "complete"

    def put(p):
        key = f"{project}/{os.path.relpath(p, local)}"
        s3.upload_file(p, MIRROR_BUCKET, key)
        if s3.head_object(Bucket=MIRROR_BUCKET, Key=key)["ContentLength"] != os.path.getsize(p):
            raise RuntimeError(f"upload size mismatch for {key}")
        return key, os.path.getsize(p)

    with cf.ThreadPoolExecutor(8) as ex:
        uploaded = dict(ex.map(put, paths[:-1]))
    k, v = put(paths[-1])
    uploaded[k] = v
    return uploaded


def verify_project(project: str) -> bool:
    s3 = boto3.client("s3", region_name=MIRROR_REGION, config=Config(retries={"max_attempts": 10, "mode": "adaptive"}, max_pool_connections=32))
    rep = json.loads(s3.get_object(Bucket=MIRROR_BUCKET, Key=f"{project}/BUILD_REPORT.json")["Body"].read())
    ok = True
    import pyarrow.fs as pafs
    fs = pafs.S3FileSystem(region=MIRROR_REGION)
    for name, t in rep["tables"].items():
        meta = pq.read_metadata(f"{MIRROR_BUCKET}/{project}/{name}.parquet", filesystem=fs)
        good = meta.num_rows == t["rows"] and meta.num_columns == t["columns"]
        ok &= good
        print(f"  {'ok ' if good else 'BAD'} {project}/{name}.parquet rows={meta.num_rows} cols={meta.num_columns} (report {t['rows']}x{t['columns']})")
    return ok


# -------------------------------------------------------------------------- main

def run_project(project: str, man_all: pd.DataFrame, aliquot_sample: dict, args, release: str) -> dict:
    t0 = time.time()
    out = os.path.join(args.out, project)
    os.makedirs(out, exist_ok=True)
    man = man_all[man_all.project_id == project].copy()
    if man.empty:
        raise SystemExit(f"no open files for {project}")
    fetcher = Fetcher(args.workers)
    report = {"project": project, "gdc_data_release": release, "built_at": datetime.now(timezone.utc).isoformat(),
              "tables": {}, "source_files_in_project": int(len(man)), "source_bucket": man.s3_bucket.iloc[0]}
    groups = args.tables.split(",")
    print(f"\n== {project}: {len(man)} open files, bucket {man.s3_bucket.iloc[0]}", flush=True)
    report["tables"]["MANIFEST"] = write_parquet(man.drop(columns=["program"]), os.path.join(out, "MANIFEST.parquet"))
    rows = man.to_dict("records")
    n_before = 0
    if "clinical" in groups:
        build_clinical(project, out, report)
    if "rna" in groups:
        build_rna([r for r in rows if r["data_type"] == "Gene Expression Quantification"], fetcher, out, report, args.workers, aliquot_sample)
    if "maf" in groups:
        build_maf([r for r in rows if r["data_type"] == "Masked Somatic Mutation"], fetcher, out, report, args.workers)
    if "cnv" in groups:
        build_cnv(rows, fetcher, out, report, args.workers, aliquot_sample)
    if "mirna" in groups:
        build_mirna(rows, fetcher, out, report, args.workers, aliquot_sample)
    if "protein" in groups:
        build_protein(rows, fetcher, out, report, args.workers, aliquot_sample)
    if "methylation" in groups:
        build_methylation(rows, fetcher, out, report, args.workers, aliquot_sample)
    report["files_read"] = int(sum(1 for r in rows if r["data_type"] in {
        "Gene Expression Quantification", "Masked Somatic Mutation", "Gene Level Copy Number", "Copy Number Segment",
        "Masked Copy Number Segment", "Allele-specific Copy Number Segment", "miRNA Expression Quantification",
        "Isoform Expression Quantification", "Protein Expression Quantification"} | ({"Methylation Beta Value"} if "methylation" in groups else set())))
    report["bytes_read"] = fetcher.bytes_read
    report["gdc_api_fallbacks"] = fetcher.api_fallbacks
    report["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(out, "README.md"), "w") as fh:
        fh.write(render_readme(project, man, report, release))
    with open(os.path.join(out, "BUILD_REPORT.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"   built {len(report['tables'])} tables in {report['elapsed_s']}s, {report['bytes_read'] / 1e9:.1f} GB read", flush=True)
    if args.upload:
        up = upload_dir(out, project)
        print(f"   uploaded {len(up)} objects, {sum(up.values()) / 1e6:.1f} MB -> s3://{MIRROR_BUCKET}/{project}/", flush=True)
    return report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--projects", required=True, help="comma list, or all-bucket, or all")
    ap.add_argument("--tables", default="clinical,rna,maf,cnv,mirna,protein")
    ap.add_argument("--methylation", action="store_true", help="also build methylation betas (large)")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--upload", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--skip-built", action="store_true", help="skip projects that already have BUILD_REPORT.json in the mirror")
    ap.add_argument("--upload-only", action="store_true", help="upload already-built local project dirs under --out, no rebuild")
    args = ap.parse_args()
    if args.methylation and "methylation" not in args.tables:
        args.tables += ",methylation"

    man_all, _assoc, aliquot_sample = load_manifest(args.manifest_dir)
    bucket_map = json.load(open(os.path.join(args.manifest_dir, "project_bucket_map.json")))
    if args.projects == "all-bucket":
        projects = sorted(p for p, b in bucket_map.items() if b and "PARTIAL" not in b)
    elif args.projects == "all":
        projects = sorted(bucket_map)
    else:
        projects = args.projects.split(",")

    if args.verify:
        ok = all(verify_project(p) for p in projects)
        return 0 if ok else 1

    s3 = boto3.client("s3", region_name=MIRROR_REGION, config=Config(retries={"max_attempts": 10, "mode": "adaptive"}, max_pool_connections=32))
    if args.upload_only:
        for p in projects:
            local = os.path.join(args.out, p)
            if not os.path.exists(os.path.join(local, "BUILD_REPORT.json")):
                print(f"== {p}: no local build, skipping")
                continue
            if args.skip_built:
                try:
                    s3.head_object(Bucket=MIRROR_BUCKET, Key=f"{p}/BUILD_REPORT.json")
                    print(f"== {p}: already in mirror, skipping")
                    continue
                except s3.exceptions.ClientError:
                    pass
            up = upload_dir(local, p)
            print(f"== {p}: uploaded {len(up)} objects, {sum(up.values()) / 1e6:.1f} MB -> s3://{MIRROR_BUCKET}/{p}/", flush=True)
        return 0
    with urllib.request.urlopen(f"{GDC_API}/status", timeout=60) as r:
        release = json.load(r)["data_release"]
    for p in projects:
        if args.skip_built:
            if os.path.exists(os.path.join(args.out, p, "BUILD_REPORT.json")):
                print(f"== {p}: already built locally, skipping")
                continue
            try:
                s3.head_object(Bucket=MIRROR_BUCKET, Key=f"{p}/BUILD_REPORT.json")
                print(f"== {p}: already in mirror, skipping")
                continue
            except s3.exceptions.ClientError:
                pass
        run_project(p, man_all, aliquot_sample, args, release)
    return 0


if __name__ == "__main__":
    sys.exit(main())
