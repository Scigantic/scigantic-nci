#!/usr/bin/env python3
"""Re-key tumor/normal pair ids in already-built copy_number tables with the tissue_type-first tumor_aliquot().

Builds before the fix named a paired file after BOTH aliquots ("A|B") whenever the normal's sample_type was missing
from NORMAL_SAMPLE_TYPES (Granulocytes, Saliva, Slides, FFPE Scrolls in CGCI-BLGSP and HCMI-CMDC). This rewrites,
in place, under <mirror-dir>/<PROJECT>/:
  copy_number/gene_level_<wf>.parquet            matrix column names (values untouched)
  copy_number/gene_level_<wf>_samples.parquet    `column` plus sample_submitter_id, sample_type, tissue_type,
                                                 tumor_descriptor, case_submitter_id of the resolved aliquot
  copy_number/*segments_<wf>.parquet             aliquot_submitter_id values (all other columns untouched)
  BUILD_REPORT.json                              rows/columns/bytes of rewritten tables, duplicate_aliquots
Column names go through the builder's unique_columns(), so a renamed column that collides with an existing one gets
the `__<file_id[:8]>` suffix exactly as a fresh build would; collisions are printed. Tables with nothing to change
are not rewritten. The directory may hold a subset of the project (only the files to repair plus BUILD_REPORT.json).

Usage: fix_pair_columns.py <manifest-dir> <mirror-dir> PROJECT [PROJECT...] [--dry-run]
"""
import argparse
import glob
import json
import os
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_gdc_project_tables import load_manifest, tumor_aliquot, unique_columns, write_table  # noqa: E402

META_FIELDS = ("sample_submitter_id", "sample_type", "tissue_type", "tumor_descriptor", "case_submitter_id")
SEGMENT_ROW_GROUP = 200_000  # build_cnv writes segment tables with row_group_size=200_000


def resolve(pair: str, aliquot_sample: dict) -> str:
    return tumor_aliquot(pair.split("|"), aliquot_sample) if pair and "|" in pair else pair


def set_duplicates(report: dict, table: str, scratch: dict) -> None:
    dups = report.get("duplicate_aliquots", {})
    dups.pop(table, None)
    if table in scratch.get("duplicate_aliquots", {}):
        dups[table] = scratch["duplicate_aliquots"][table]
    if dups and "duplicate_aliquots" not in report:
        # the builder adds the key while building tables, i.e. before files_read; key order is cosmetic
        items = list(report.items())
        at = next((i for i, (k, _) in enumerate(items) if k in ("realigned_files", "promoted_to_float", "files_read")), len(items))
        report.clear()
        report.update(items[:at] + [("duplicate_aliquots", dups)] + items[at:])
    elif not dups:
        report.pop("duplicate_aliquots", None)


def fix_gene_level(proj_dir: str, path: str, aliquot_sample: dict, report: dict, dry: bool) -> None:
    rel = os.path.relpath(path, proj_dir)[: -len(".parquet")]
    spath = path[: -len(".parquet")] + "_samples.parquet"
    if not os.path.exists(spath):
        raise SystemExit(f"{path}: no {os.path.basename(spath)} next to it; both are needed")
    smp_t = pq.read_table(spath)
    smp = smp_t.to_pandas()
    names = pq.read_schema(path).names
    n_feat = len(names) - len(smp)
    old = names[n_feat:]
    if old != smp["column"].tolist():
        raise SystemExit(f"{rel}: matrix columns do not match {os.path.basename(spath)} `column` in order")
    ids = [resolve(a, aliquot_sample) for a in smp["aliquot_submitter_id"]]
    scratch: dict = {}
    new = unique_columns(ids, smp["file_id"].tolist(), scratch, rel)
    changed = [i for i, (o, n) in enumerate(zip(old, new)) if o != n]
    pipes_before = sum("|" in c for c in old)
    collisions = [(old[i], new[i]) for i in changed if "__" in new[i] and "__" not in old[i]]
    print(f"  {rel}: {len(old)} columns, {pipes_before} pipe-joined -> {sum('|' in c for c in new)}, {len(changed)} renamed, "
          f"{len(collisions)} collision-suffixed")
    for o, n in collisions:
        print(f"    COLLISION {o} -> {n}")
    if not changed:
        return
    for i in changed:
        a = new[i].split("__")[0]
        meta = aliquot_sample.get(a)
        if meta is None:
            raise SystemExit(f"{rel}: resolved aliquot {a} has no sample metadata in the manifest")
        smp.at[i, "column"] = new[i]
        for k in META_FIELDS:
            smp.at[i, k] = meta.get(k)
    if len(set(new)) != len(new):
        raise SystemExit(f"{rel}: column names are not unique after the rename")
    if dry:
        return
    mat = pq.read_table(path)
    report["tables"][rel] = write_table(mat.rename_columns(names[:n_feat] + new), path)
    # cast back to the file's own schema (incl. pandas metadata) so the pandas version running this cannot change types
    report["tables"][f"{rel}_samples"] = write_table(pa.Table.from_pandas(smp, preserve_index=False).cast(smp_t.schema), spath)
    set_duplicates(report, rel, scratch)


def fix_segments(proj_dir: str, path: str, aliquot_sample: dict, report: dict, dry: bool) -> None:
    rel = os.path.relpath(path, proj_dir)[: -len(".parquet")]
    col = pq.read_table(path, columns=["aliquot_submitter_id"]).column(0).combine_chunks()
    uniq = pc.unique(col).to_pylist()
    mapping = {u: resolve(u, aliquot_sample) for u in uniq if u and "|" in u}
    rows = int(pc.sum(pc.is_in(col, value_set=pa.array(list(mapping), pa.string()))).as_py() or 0) if mapping else 0
    left = sum("|" in v for v in mapping.values())
    print(f"  {rel}: {len(col)} rows, {len(mapping)} pipe-joined ids on {rows} rows -> {left} pipe-joined ids left")
    if not mapping or dry:
        return
    t = pq.read_table(path)
    j = t.schema.get_field_index("aliquot_submitter_id")
    fixed = pa.array([mapping.get(v, v) for v in t.column(j).to_pylist()], type=t.schema.field(j).type)
    t = t.set_column(j, t.schema.field(j), fixed)
    report["tables"][rel] = write_table(t, path, row_group_size=SEGMENT_ROW_GROUP)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest_dir")
    ap.add_argument("mirror_dir")
    ap.add_argument("projects", nargs="+")
    ap.add_argument("--dry-run", action="store_true", help="print what would change, write nothing")
    args = ap.parse_args()
    _, _, aliquot_sample = load_manifest(args.manifest_dir)
    for p in args.projects:
        proj_dir = os.path.join(args.mirror_dir, p)
        rep_path = os.path.join(proj_dir, "BUILD_REPORT.json")
        report = json.load(open(rep_path))
        before = json.dumps(report, indent=1)
        print(f"== {p}")
        for path in sorted(glob.glob(os.path.join(proj_dir, "copy_number", "gene_level_*.parquet"))):
            if not path.endswith("_samples.parquet"):
                fix_gene_level(proj_dir, path, aliquot_sample, report, args.dry_run)
        for path in sorted(glob.glob(os.path.join(proj_dir, "copy_number", "*segments_*.parquet"))):
            fix_segments(proj_dir, path, aliquot_sample, report, args.dry_run)
        for path in sorted(glob.glob(os.path.join(proj_dir, "*", "*samples.parquet"))):
            if "/copy_number/" not in path and any("|" in str(c) for c in pq.read_table(path, columns=["column"]).column(0).to_pylist()):
                print(f"  WARNING {os.path.relpath(path, proj_dir)} has pipe-joined columns; not handled here")
        if not args.dry_run and json.dumps(report, indent=1) != before:
            with open(rep_path, "w") as fh:
                json.dump(report, fh, indent=1)
            print(f"  wrote {rep_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
