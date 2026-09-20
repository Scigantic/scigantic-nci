#!/usr/bin/env python3
"""
Dump the GDC open-access file manifest that build_gdc_project_tables.py consumes.

Writes into --out:
  gdc_open_files.ndjson    every open-access file (id, name, type, size, md5, case/sample)
  gdc_assoc.ndjson         aliquot associations + workflow/platform for the matrix data types
  project_bucket_map.json  which public AWS bucket serves each project (HEAD-checked on
                           three files per project; null = only the GDC HTTPS API)

Usage:
  python dump_gdc_manifest.py --out <dir>

All three are needed because the GDC file endpoint does not return aliquot
barcodes in the same query as the bulk listing without blowing the page size,
and because the AWS bucket that mirrors a project is not recorded anywhere in
the GDC API (it comes from the Registry of Open Data listings).
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
import urllib.parse
import urllib.request

import boto3
from botocore import UNSIGNED
from botocore.config import Config

GDC_API = "https://api.gdc.cancer.gov"
MATRIX_TYPES = [
    "Gene Expression Quantification", "Methylation Beta Value", "Gene Level Copy Number",
    "Copy Number Segment", "Masked Copy Number Segment", "Allele-specific Copy Number Segment",
    "miRNA Expression Quantification", "Isoform Expression Quantification",
    "Protein Expression Quantification", "Masked Somatic Mutation",
]
# Public GDC buckets listed on registry.opendata.aws (cancer tag) as of 2026-09-13.
CANDIDATE_BUCKETS = [
    "tcga-2-open", "gdc-target-phs000218-2-open", "gdc-cptac-phs001287-2-open", "gdc-cptac-2-phs000892-2-open",
    "gdc-ccle-2-open", "gdc-cgci-phs000235-2-open", "gdc-cgci-blgsp-phs000235-2-open", "gdc-hcmi-cmdc-phs001486-2-open",
    "gdc-mmrf-commpass-phs000748-2-open", "gdc-beataml1-cohort-phs001657-2-open", "gdc-beataml1.0-crenolanib-phs001628-2-open",
    "gdc-fm-ad-phs001179-2-open", "gdc-ctsp-phs001175-2-open", "gdc-nciccr-phs001444-2-open", "gdc-wcdt-mcrpc-phs001648-2-open",
    "gdc-organoid-pancreatic-phs001611-2-open", "gdc-ohsu-cnl-phs001799-2-open", "gdc-cddp-eagle-1-phs001239-2-open",
    "gdc-exceptional-responders-er-phs001145-2-open", "gdc-mp2prt-wt-phs001965-2-open",
]


def page(endpoint: str, filters: dict, fields: str, out_path: str) -> int:
    n, frm, size = 0, 0, 10000
    with open(out_path, "w") as out:
        while True:
            q = urllib.parse.urlencode({"filters": json.dumps(filters), "fields": fields, "size": size, "from": frm, "format": "json"})
            for attempt in range(6):
                try:
                    with urllib.request.urlopen(f"{GDC_API}/{endpoint}?{q}", timeout=600) as r:
                        d = json.load(r)
                    break
                except Exception as e:  # noqa: BLE001
                    print("retry", frm, e, file=sys.stderr)
                    time.sleep(10 * (attempt + 1))
            hits = d["data"]["hits"]
            for h in hits:
                out.write(json.dumps(h) + "\n")
            n += len(hits)
            total = d["data"]["pagination"]["total"]
            print(f"{out_path}: {n}/{total}", file=sys.stderr, flush=True)
            if n >= total or not hits:
                break
            frm += size
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    open_filter = {"op": "=", "content": {"field": "access", "value": "open"}}
    page("files", open_filter,
         "file_id,file_name,data_type,data_category,data_format,file_size,experimental_strategy,md5sum,"
         "cases.case_id,cases.submitter_id,cases.project.project_id,cases.project.program.name,"
         "cases.samples.sample_type,cases.samples.submitter_id",
         os.path.join(args.out, "gdc_open_files.ndjson"))
    page("files", {"op": "and", "content": [open_filter, {"op": "in", "content": {"field": "data_type", "value": MATRIX_TYPES}}]},
         "file_id,data_type,associated_entities.entity_submitter_id,associated_entities.entity_type,associated_entities.case_id,"
         "cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.tissue_type,"
         "cases.samples.tumor_descriptor,cases.samples.portions.analytes.aliquots.submitter_id,analysis.workflow_type,platform",
         os.path.join(args.out, "gdc_assoc.ndjson"))

    # bucket audit: 3 random files per project against the candidate buckets
    by_project: dict[str, list[tuple[str, str]]] = {}
    with open(os.path.join(args.out, "gdc_open_files.ndjson")) as fh:
        for line in fh:
            h = json.loads(line)
            projs = sorted({c.get("project", {}).get("project_id", "?") for c in h.get("cases") or [{}]})
            by_project.setdefault(projs[0], []).append((h["file_id"], h["file_name"]))
    s3 = boto3.client("s3", region_name="us-east-1", config=Config(signature_version=UNSIGNED))

    def head(b, k):
        try:
            s3.head_object(Bucket=b, Key=k)
            return True
        except Exception:  # noqa: BLE001
            return False

    random.seed(7)
    result = {}
    for proj, files in sorted(by_project.items()):
        sample = random.sample(files, min(3, len(files)))
        toks = [t.lower() for t in proj.replace("_", "-").split("-")]
        cands = sorted(CANDIDATE_BUCKETS, key=lambda b: -sum(t in b for t in toks))
        found = None
        for b in cands:
            ok = [head(b, f"{fid}/{name}") for fid, name in sample]
            if all(ok):
                found = b
                break
            if any(ok):
                found = b + " (PARTIAL)"
                break
        result[proj] = found
        print(f"{proj:28s} {len(files):7d} files -> {found}", flush=True)
    with open(os.path.join(args.out, "project_bucket_map.json"), "w") as fh:
        json.dump(result, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
