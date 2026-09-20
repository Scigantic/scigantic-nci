#!/usr/bin/env python3
"""
Build per-collection archives for the NCI Imaging Data Commons (IDC) and
upload them to s3://scigantic-idc-open/<collection_id>/.

Why this exists
---------------
IDC's public buckets (idc-open-data, idc-open-data-two, idc-open-data-cr) hold
about a million DICOM series under one flat UUID namespace; the only way to
know what a folder is comes from the idc-index metadata tables. A per-collection
archive puts that metadata next to a browsable sample of the images, so a
collection can be explored (patients, studies, modalities, slide stains,
clinical tables) and one series pulled on demand, without BigQuery or the IDC
portal.

What it builds (per collection)
--------------------------------
  series.parquet          every DICOM series of the collection (idc-index rows:
                          patient, study, modality, body part, instance count,
                          size, license, DOI, aws bucket + series UUID)
  studies.parquet         one row per study: patient, date, modalities, sizes
  patients.parquet        one row per patient: sex, age, studies, modalities, size
  sm_series.parquet       slide-microscopy attributes (stain, fixative, objective,
                          pixel spacing, anatomy) when the collection has SM
  analysis_results.parquet  derived-data collections (segmentations, annotations)
                          that cover this collection
  clinical/<table>.parquet + clinical/dictionary.parquet   when IDC has them
  sample/<Modality>_<SeriesInstanceUID>/*.dcm   a small real series per modality
                          (whole series for radiology; thumbnail + lowest
                          pyramid levels for slide microscopy)
  samples.parquet         one row per mirrored sample with its license
  README.md, BUILD_REPORT.json

Sources: idc-index (pinned in BUILD_REPORT with the IDC data version), reads
from the IDC buckets anonymously. Sample bytes are copied verbatim, and only
from series whose license allows redistribution: CC BY-NC series are never
picked as samples (Scigantic is a commercial redistributor), so a collection
whose eligible series are all CC BY-NC gets the index tables and no sample.

Usage
-----
  python build_idc_collection_archives.py --out <workdir> --collections tcga_brca,nlst [--upload]
  python build_idc_collection_archives.py --out <workdir> --collections all [--upload] [--skip-built]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone

import concurrent.futures as cf

import boto3
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from botocore import UNSIGNED
from botocore.config import Config

MIRROR_BUCKET = "scigantic-idc-open"
MIRROR_REGION = "us-east-1"
SAMPLE_MODALITY_ORDER = ["CT", "MR", "PT", "MG", "DX", "CR", "US", "NM", "XA", "RTSTRUCT", "SEG", "SM"]
RADIOLOGY_MAX_MB = 80
SM_INSTANCE_MAX_BYTES = 12_000_000
MAX_SAMPLE_MODALITIES = 3
SAMPLE_LICENSE_POLICY = ("Sample series are copied only from series whose license allows redistribution (CC BY, the NLM terms); "
                         "CC BY-NC series are never mirrored. samples.parquet lists every mirrored sample with its license_short_name.")
NLM_COURTESY = "Courtesy of the U.S. National Library of Medicine"
SAMPLE_COLUMNS = ["modality", "SeriesInstanceUID", "PatientID", "StudyInstanceUID", "series_aws_url", "kind", "files", "bytes",
                  "instanceCount_full_series", "series_size_MB_full", "license_short_name", "source_DOI", "folder"]


def write_parquet(df: pd.DataFrame, path: str) -> dict:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    table = pa.Table.from_pandas(df.reset_index(drop=True), preserve_index=False)
    pq.write_table(table, path, compression="zstd")
    return {"rows": table.num_rows, "columns": table.num_columns, "bytes": os.path.getsize(path)}


def joined(s: pd.Series) -> str:
    return ", ".join(sorted({str(x) for x in s.dropna() if str(x)}))


class Idc:
    def __init__(self):
        from idc_index import IDCClient
        self.c = IDCClient()
        for name in ("collections_index", "clinical_index", "sm_index", "sm_instance_index", "analysis_results_index"):
            self.c.fetch_index(name)
        self.idx = self.c.index
        self.collections = self.c.collections_index
        self.clinical = self.c.clinical_index
        self.sm = self.c.sm_index
        self.smi = self.c.sm_instance_index
        self.ar = self.c.analysis_results_index
        self.version = self.c.get_idc_version()
        import importlib.metadata as im
        self.pkg = {"idc-index": im.version("idc-index"), "idc-index-data": im.version("idc-index-data")}
        self.s3 = boto3.client("s3", region_name="us-east-1", config=Config(signature_version=UNSIGNED, max_pool_connections=32))


def redistributable(series: pd.DataFrame) -> pd.DataFrame:
    """Series we may copy: anything not CC BY-NC. CC BY first so a mixed collection samples its CC BY series."""
    lic = series.license_short_name.fillna("")
    return series[~lic.str.startswith("CC BY-NC")].assign(_license_rank=(~lic.str.startswith("CC BY")).astype(int))


def pick_samples(series: pd.DataFrame, idc: Idc) -> list[dict]:
    """One small real series per modality (up to MAX_SAMPLE_MODALITIES), radiology first, never a CC BY-NC series."""
    picks = []
    series = redistributable(series)
    present = [m for m in SAMPLE_MODALITY_ORDER if m in set(series.Modality)]
    for mod in present:
        if len(picks) >= MAX_SAMPLE_MODALITIES:
            break
        sub = series[series.Modality == mod]
        sub = sub[sub._license_rank == sub._license_rank.min()]
        if mod == "SM":
            sm_rows = idc.sm[idc.sm.SeriesInstanceUID.isin(sub.SeriesInstanceUID)]
            # smallest slide by pixel matrix, so the low levels are quick to read
            if len(sm_rows):
                target = sm_rows.sort_values("max_TotalPixelMatrixColumns").iloc[0].SeriesInstanceUID
                row = sub[sub.SeriesInstanceUID == target].iloc[0]
            else:
                row = sub.sort_values("series_size_MB").iloc[0]
            inst = idc.smi[idc.smi.SeriesInstanceUID == row.SeriesInstanceUID].sort_values("instance_size")
            inst = inst[inst.instance_size <= SM_INSTANCE_MAX_BYTES]
            if inst.empty:
                continue
            picks.append({"modality": mod, "series": row, "instances": inst.head(3), "kind": "sm_levels"})
        else:
            cand = sub[(sub.series_size_MB <= RADIOLOGY_MAX_MB)]
            if mod in ("CT", "MR", "PT"):
                cand = cand[cand.instanceCount >= 10]
            if cand.empty:
                continue
            # prefer uncompressed transfer syntaxes (readable by pydicom with no decoder plugin), then the most
            # instances among small series, not a degenerate 10-slice stub
            cand = cand.assign(_compressed=cand.transfer_syntax_name.fillna("").str.contains("JPEG|RLE", regex=True))
            cand = cand.sort_values(["_compressed", "instanceCount", "series_size_MB"], ascending=[True, False, True])
            cand = cand[cand.series_size_MB <= max(cand.series_size_MB.min() * 3, 8)]
            cand = cand.sort_values(["_compressed", "instanceCount", "series_size_MB"], ascending=[True, False, True])
            picks.append({"modality": mod, "series": cand.iloc[0], "instances": None, "kind": "series"})
    return picks


def download_sample(pick: dict, out_dir: str, idc: Idc) -> dict:
    row = pick["series"]
    folder = os.path.join(out_dir, "sample", f"{pick['modality']}_{row.SeriesInstanceUID}")
    os.makedirs(folder, exist_ok=True)
    bucket, prefix = row.aws_bucket, row.crdc_series_uuid + "/"
    n, total = 0, 0
    if pick["kind"] == "series":
        keys, token = [], None
        while True:
            kw = {"Bucket": bucket, "Prefix": prefix, "MaxKeys": 1000}
            if token:
                kw["ContinuationToken"] = token
            resp = idc.s3.list_objects_v2(**kw)
            keys += [(o["Key"], o["Size"]) for o in resp.get("Contents", []) if not o["Key"].endswith("/")]  # skip the folder marker object
            token = resp.get("NextContinuationToken")
            if not token:
                break

        def get(item):
            key, size = item
            idc.s3.download_file(bucket, key, os.path.join(folder, os.path.basename(key)))
            return size

        with cf.ThreadPoolExecutor(8) as ex:
            for size in ex.map(get, keys):
                n += 1
                total += size
    else:
        for _, inst in pick["instances"].iterrows():
            key = f"{prefix}{inst.crdc_instance_uuid}.dcm"
            local = os.path.join(folder, f"{inst.crdc_instance_uuid}.dcm")
            idc.s3.download_file(bucket, key, local)
            n += 1
            total += os.path.getsize(local)
    return {"modality": pick["modality"], "SeriesInstanceUID": row.SeriesInstanceUID, "PatientID": row.PatientID,
            "StudyInstanceUID": row.StudyInstanceUID, "series_aws_url": row.series_aws_url, "kind": pick["kind"],
            "files": n, "bytes": total, "instanceCount_full_series": int(row.instanceCount), "series_size_MB_full": float(row.series_size_MB),
            "license": row.license_short_name, "source_DOI": row.source_DOI, "folder": os.path.relpath(folder, out_dir)}


def render_readme(cid: str, cmeta: pd.Series, series: pd.DataFrame, report: dict, clinical_tables: list[str], analysis: pd.DataFrame) -> str:
    by_mod = series.groupby("Modality").agg(series=("SeriesInstanceUID", "count"), patients=("PatientID", "nunique"),
                                            GB=("series_size_MB", lambda s: round(s.sum() / 1024, 2))).sort_values("series", ascending=False)
    lic = series.groupby("license_short_name").size().sort_values(ascending=False)
    dois = sorted(series.source_DOI.dropna().unique())
    lines = [f"# {cmeta.collection_name} ({cid}) on the NCI Imaging Data Commons", "",
             f"Built by Scigantic from idc-index {report['idc_index_version']} (IDC {report['idc_data_version']}) on {report['built_at'][:10]}.",
             "", f"**Cancer types:** {cmeta.cancer_types}  ", f"**Tumor locations:** {cmeta.tumor_locations}  ",
             f"**Subjects (IDC):** {cmeta.subjects}  ", f"**Program:** {cmeta.program_id}  ", f"**Supporting data (IDC):** {cmeta.supporting_data}  ",
             f"**Status:** {cmeta.status}, updated {cmeta.updated}", "", "## Description (IDC)", "", str(cmeta.description or "").strip(), "",
             "## What is in the collection", "", f"{len(series):,} DICOM series, {series.PatientID.nunique():,} patients, {series.StudyInstanceUID.nunique():,} studies, {series.series_size_MB.sum() / 1024:.1f} GB.", "",
             "| modality | series | patients | GB |", "|---|---|---|---|"]
    for mod, r in by_mod.iterrows():
        lines.append(f"| {mod} | {int(r.series):,} | {int(r.patients):,} | {r.GB} |")
    lines += ["", "## Licenses", ""] + [f"- {k}: {v:,} series" for k, v in lic.items()]
    lines += ["", "## Source DOIs", ""] + [f"- https://doi.org/{d}" for d in dois]
    lines += ["", "## Tables in this archive", "", "| path | rows | columns |", "|---|---|---|"]
    for name, t in sorted(report["tables"].items()):
        lines.append(f"| {name}.parquet | {t['rows']:,} | {t['columns']} |")
    if clinical_tables:
        lines += ["", "IDC clinical tables for this collection: " + ", ".join(clinical_tables) + ". `clinical/dictionary.parquet` gives the label and coded values of every column."]
    if len(analysis):
        lines += ["", "## Analysis results covering this collection", ""] + [f"- {r.analysis_result_title} ({r.analysis_result_id}): {r.modalities}, {r.license_short_name}, https://doi.org/{r.source_DOI}" for _, r in analysis.iterrows()]
    lines += ["", "## Sample images", ""]
    if not report["samples"]:
        lines.append("No sample series is mirrored: no series of a sampled modality is both small enough and under a license that allows "
                     f"redistribution ({', '.join(lic.index)}). Pull one from the IDC bucket (below) and respect its license"
                     + ("; CC BY-NC series are for non-commercial use only." if any(k.startswith("CC BY-NC") for k in lic.index) else "."))
    for s in report["samples"]:
        lines.append(f"- `{s['folder']}`: {s['modality']} series of patient {s['PatientID']}, {s['files']} file(s), {s['bytes'] / 1e6:.1f} MB"
                     + (f" (thumbnail and lowest pyramid levels only; the full slide is {s['series_size_MB_full']:.0f} MB in {s['instanceCount_full_series']} instances)" if s["kind"] == "sm_levels" else f" (the complete series)")
                     + f", {s['license']}")
    if report["samples"]:
        lines.append(f"\n{SAMPLE_LICENSE_POLICY}")
    lines += ["", "## Getting the rest", "",
              "Every row of `series.parquet` has `series_aws_url` (`s3://idc-open-data/<crdc_series_uuid>/*`). The buckets are public and anonymous reads work:", "",
              "```python", "import boto3, pandas as pd", "from botocore import UNSIGNED", "from botocore.config import Config",
              "s = pd.read_parquet('/mnt/archive/series.parquet')", "row = s[s.Modality == 'CT'].iloc[0]",
              "s3 = boto3.client('s3', region_name='us-east-1', config=Config(signature_version=UNSIGNED))",
              "for o in s3.list_objects_v2(Bucket=row.aws_bucket, Prefix=row.crdc_series_uuid + '/')['Contents']:",
              "    s3.download_file(row.aws_bucket, o['Key'], '/tmp/' + o['Key'].split('/')[-1])", "```", "",
              "For bulk pulls use `s5cmd` or `pip install idc-index` (`IDCClient().download_from_selection(...)`).", "",
              "## Attribution", "",
              "Cite the collection DOI(s) above and the IDC (Fedorov et al., Radiographics 2021, doi:10.1148/rg.2021210020; Fedorov et al., Nature Methods 2023, doi:10.1038/s41592-023-01893-2). Respect the per-series license in `series.parquet` (some collections mix CC BY and CC BY-NC)."]
    if any(k.startswith("National Library of Medicine") for k in lic.index):
        lines.append(f"Series under the NLM terms carry the attribution \"{NLM_COURTESY}\".")
    lines.append("")
    return "\n".join(lines)


def build_collection(cid: str, idc: Idc, out_root: str, upload: bool) -> dict:
    t0 = time.time()
    out = os.path.join(out_root, cid)
    os.makedirs(out, exist_ok=True)
    series = idc.idx[idc.idx.collection_id == cid].copy()
    if series.empty:
        raise SystemExit(f"no series for collection {cid}")
    cmeta = idc.collections[idc.collections.collection_id == cid].iloc[0]
    report = {"collection_id": cid, "collection_name": cmeta.collection_name, "idc_data_version": idc.version,
              "idc_index_version": idc.pkg["idc-index"], "idc_index_data_version": idc.pkg["idc-index-data"],
              "built_at": datetime.now(timezone.utc).isoformat(), "tables": {}, "samples": [], "sample_license_policy": SAMPLE_LICENSE_POLICY}
    print(f"\n== {cid} ({cmeta.collection_name}): {len(series):,} series, {series.PatientID.nunique():,} patients", flush=True)

    report["tables"]["series"] = write_parquet(series, os.path.join(out, "series.parquet"))
    studies = series.groupby("StudyInstanceUID").agg(PatientID=("PatientID", "first"), StudyDate=("StudyDate", "first"), StudyDescription=("StudyDescription", "first"),
                                                     modalities=("Modality", joined), n_series=("SeriesInstanceUID", "count"), instances=("instanceCount", "sum"),
                                                     size_MB=("series_size_MB", "sum"), body_parts=("BodyPartExamined", joined)).reset_index()
    report["tables"]["studies"] = write_parquet(studies, os.path.join(out, "studies.parquet"))
    patients = series.groupby("PatientID").agg(PatientSex=("PatientSex", lambda s: joined(s)), PatientAge=("PatientAge", "first"), n_studies=("StudyInstanceUID", "nunique"),
                                               n_series=("SeriesInstanceUID", "count"), modalities=("Modality", joined), size_MB=("series_size_MB", "sum")).reset_index()
    report["tables"]["patients"] = write_parquet(patients, os.path.join(out, "patients.parquet"))
    sm_rows = idc.sm[idc.sm.SeriesInstanceUID.isin(series.SeriesInstanceUID)]
    if len(sm_rows):
        report["tables"]["sm_series"] = write_parquet(sm_rows, os.path.join(out, "sm_series.parquet"))
    analysis = idc.ar[idc.ar.collections.fillna("").str.contains(rf"(?:^|,\s*){cid}(?:,|$)", regex=True)]
    if len(analysis):
        report["tables"]["analysis_results"] = write_parquet(analysis, os.path.join(out, "analysis_results.parquet"))
    clin = idc.clinical[idc.clinical.collection_id == cid]
    clinical_tables = []
    if len(clin):
        for tname in sorted(clin.short_table_name.unique()):
            try:
                df = idc.c.get_clinical_table(tname)
            except Exception as e:  # noqa: BLE001
                print(f"   clinical table {tname} failed: {e}", flush=True)
                continue
            report["tables"][f"clinical/{tname}"] = write_parquet(df, os.path.join(out, "clinical", f"{tname}.parquet"))
            clinical_tables.append(tname)
        dictionary = clin[["short_table_name", "column", "column_label", "values"]].copy()
        dictionary["values"] = dictionary["values"].astype(str)
        report["tables"]["clinical/dictionary"] = write_parquet(dictionary, os.path.join(out, "clinical", "dictionary.parquet"))
    for pick in pick_samples(series, idc):
        report["samples"].append(download_sample(pick, out, idc))
        print(f"   sample {pick['modality']}: {report['samples'][-1]['files']} files, {report['samples'][-1]['bytes'] / 1e6:.1f} MB, {report['samples'][-1]['license']}", flush=True)
    samples_df = pd.DataFrame([{**s, "license_short_name": s["license"]} for s in report["samples"]], columns=SAMPLE_COLUMNS)
    report["tables"]["samples"] = write_parquet(samples_df.astype({"files": "int64", "bytes": "int64", "instanceCount_full_series": "int64", "series_size_MB_full": "float64"})
                                                .astype({c: "string" for c in SAMPLE_COLUMNS if c not in ("files", "bytes", "instanceCount_full_series", "series_size_MB_full")}),
                                                os.path.join(out, "samples.parquet"))
    report["n_series"], report["n_patients"], report["n_studies"] = int(len(series)), int(series.PatientID.nunique()), int(series.StudyInstanceUID.nunique())
    report["total_GB"] = round(float(series.series_size_MB.sum()) / 1024, 2)
    report["modalities"] = series.Modality.value_counts().to_dict()
    report["licenses"] = series.license_short_name.value_counts().to_dict()
    report["source_DOIs"] = sorted(series.source_DOI.dropna().unique().tolist())
    report["clinical_tables"] = clinical_tables
    report["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(out, "README.md"), "w") as fh:
        fh.write(render_readme(cid, cmeta, series, report, clinical_tables, analysis))
    with open(os.path.join(out, "BUILD_REPORT.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    if upload:
        s3 = boto3.client("s3", region_name=MIRROR_REGION, config=Config(retries={"max_attempts": 10, "mode": "adaptive"}, max_pool_connections=32))
        n = 0
        for root, _, files in os.walk(out):
            for f in files:
                p = os.path.join(root, f)
                s3.upload_file(p, MIRROR_BUCKET, f"{cid}/{os.path.relpath(p, out)}")
                n += 1
        print(f"   uploaded {n} objects -> s3://{MIRROR_BUCKET}/{cid}/", flush=True)
    print(f"   built {len(report['tables'])} tables + {len(report['samples'])} sample series in {report['elapsed_s']}s", flush=True)
    return report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--collections", required=True, help="comma list or 'all'")
    ap.add_argument("--upload", action="store_true")
    ap.add_argument("--skip-built", action="store_true")
    ap.add_argument("--upload-only", action="store_true", help="upload already-built local collection dirs under --out, no rebuild")
    args = ap.parse_args()
    s3 = boto3.client("s3", region_name=MIRROR_REGION, config=Config(retries={"max_attempts": 10, "mode": "adaptive"}, max_pool_connections=32))
    if args.upload_only:
        cids = sorted(c for c in os.listdir(args.out) if os.path.exists(os.path.join(args.out, c, "BUILD_REPORT.json"))) if args.collections == "all" else args.collections.split(",")
        for cid in cids:
            if args.skip_built:
                try:
                    s3.head_object(Bucket=MIRROR_BUCKET, Key=f"{cid}/BUILD_REPORT.json")
                    print(f"== {cid}: already in mirror, skipping")
                    continue
                except s3.exceptions.ClientError:
                    pass
            out = os.path.join(args.out, cid)
            paths = [os.path.join(root, f) for root, _, files in os.walk(out) for f in files]
            paths.sort(key=lambda p: p.endswith("BUILD_REPORT.json"))  # report last: its presence means "complete"

            def put(p):
                key = f"{cid}/{os.path.relpath(p, out)}"
                s3.upload_file(p, MIRROR_BUCKET, key)
                if s3.head_object(Bucket=MIRROR_BUCKET, Key=key)["ContentLength"] != os.path.getsize(p):
                    raise RuntimeError(f"upload size mismatch for {key}")
                return os.path.getsize(p)

            with cf.ThreadPoolExecutor(16) as ex:
                sizes = list(ex.map(put, paths[:-1]))
            sizes.append(put(paths[-1]))
            n, total = len(sizes), sum(sizes)
            print(f"== {cid}: uploaded {n} objects, {total / 1e6:.1f} MB -> s3://{MIRROR_BUCKET}/{cid}/", flush=True)
        return 0
    idc = Idc()
    cids = sorted(idc.idx.collection_id.unique()) if args.collections == "all" else args.collections.split(",")
    for cid in cids:
        if args.skip_built:
            try:
                s3.head_object(Bucket=MIRROR_BUCKET, Key=f"{cid}/BUILD_REPORT.json")
                print(f"== {cid}: already in mirror, skipping")
                continue
            except s3.exceptions.ClientError:
                pass
        build_collection(cid, idc, args.out, args.upload)
    return 0


if __name__ == "__main__":
    sys.exit(main())
