"""The ``scigantic-nci`` console script.

Thin argparse wrappers over :mod:`scigantic_nci.gdc` and :mod:`scigantic_nci.idc`.
Tables are written as CSV to stdout, or to ``--out`` when given; ``build-index`` writes
the bucket's catalog index as parquet. Errors from the library (``NciError``,
``ValueError``) are printed to stderr and give exit code 2.

An id is treated as an IDC collection when it contains a lowercase letter
(``tcga_lihc``, ``4d_lung``) and as a GDC project otherwise (``TCGA-LIHC``,
``CPTAC-3``); that is how the two catalogs spell their ids.
"""
from __future__ import annotations

import argparse
import sys
from typing import Any, Sequence

import pandas as pd

from . import gdc, idc
from ._store import GDC_BUCKET, IDC_BUCKET, INDEX_FILES, NciError, write_index
from ._version import __version__


def _is_idc(ident: str) -> bool:
    return any(ch.islower() for ch in ident)


def _write(df: pd.DataFrame, out: str | None, index: bool) -> None:
    if out:
        df.to_csv(out, index=index)
    else:
        df.to_csv(sys.stdout, index=index, lineterminator="\n")


def _split(value: str | None) -> list[str] | None:
    if value is None:
        return None
    items = [v.strip() for v in value.split(",") if v.strip()]
    return items or None


def _series_row(collection: str, series_uid: str | None, uuid: str | None, root: str | None) -> pd.Series:
    """The series.parquet row named by SeriesInstanceUID or crdc_series_uuid."""
    if series_uid is None and uuid is None:
        raise ValueError("pass --series-uid or --uuid")
    column, value = ("SeriesInstanceUID", series_uid) if series_uid is not None else ("crdc_series_uuid", uuid)
    df = idc.series(collection, root=root)
    hit = df[df[column] == value]
    if hit.empty:
        raise ValueError(f"{collection} has no series with {column} {value!r} ({len(df)} series)")
    return hit.iloc[0]


# ----------------------------------------------------------------- commands


def _cmd_projects(a: argparse.Namespace) -> None:
    _write(gdc.projects(), a.out, index=False)


def _cmd_collections(a: argparse.Namespace) -> None:
    _write(idc.collections(), a.out, index=False)


def _cmd_build_index(a: argparse.Namespace) -> None:
    """Write PROJECTS.parquet or COLLECTIONS.parquet, built from every BUILD_REPORT.json (never from an existing index)."""
    df = gdc._build_projects(use_index=False) if a.mirror == "gdc" else idc._build_collections(use_index=False)
    size = write_index(df, a.out)
    sys.stdout.write(f"{a.out}: {len(df)} rows, {size} bytes\n")


def _cmd_tables(a: argparse.Namespace) -> None:
    df = idc.tables(a.id, root=a.root) if _is_idc(a.id) else gdc.tables(a.id, root=a.root)
    _write(df, a.out, index=False)


def _cmd_readme(a: argparse.Namespace) -> None:
    text = idc.readme(a.id, root=a.root) if _is_idc(a.id) else gdc.readme(a.id, root=a.root)
    sys.stdout.write(text if text.endswith("\n") else text + "\n")


def _cmd_clinical(a: argparse.Namespace) -> None:
    _write(gdc.clinical(a.project, primary_only=not a.all_diagnoses, root=a.root), a.out, index=False)


def _cmd_expression(a: argparse.Namespace) -> None:
    df = gdc.expression(
        a.project,
        genes=_split(a.genes),
        sample_type=_split(a.sample_type),
        kind="counts" if a.counts else "tpm",
        root=a.root,
    )
    _write(df, a.out, index=True)


def _cmd_mutations(a: argparse.Namespace) -> None:
    df = gdc.mutations(a.project, genes=_split(a.genes), non_silent=a.non_silent, root=a.root)
    _write(df, a.out, index=False)


def _cmd_mutation_frequency(a: argparse.Namespace) -> None:
    freq = gdc.mutation_frequency(a.project, non_silent=not a.include_silent, root=a.root)
    if a.top is not None:
        freq = freq.head(a.top)
    _write(freq.to_frame(), a.out, index=True)


def _cmd_series(a: argparse.Namespace) -> None:
    df = idc.series(a.collection, modality=a.modality, max_size_mb=a.max_size_mb, root=a.root)
    _write(df, a.out, index=False)


def _cmd_pull_series(a: argparse.Namespace) -> None:
    row = _series_row(a.collection, a.series_uid, a.uuid, a.root)
    for path in idc.pull_series(row, a.dest):
        sys.stdout.write(path + "\n")


def _cmd_sample(a: argparse.Namespace) -> None:
    df = idc.samples(a.collection, root=a.root)
    if df.empty:
        raise ValueError(f"{a.collection} has no sample series")
    if a.modality is not None:
        df = df[df["modality"] == a.modality]
        if df.empty:
            have = ", ".join(idc.samples(a.collection, root=a.root)["modality"].astype(str))
            raise ValueError(f"{a.collection} has no {a.modality!r} sample; sample modalities: {have}")
    row = df.iloc[0]
    for key, value in row.items():
        sys.stdout.write(f"{key}: {value}\n")
    files = idc.sample_files(a.collection, modality=str(row["modality"]), root=a.root)
    sys.stdout.write(f"dcm_files: {len(files)}\n")
    for path in files:
        sys.stdout.write(f"  {path}\n")


def _cmd_viewer_url(a: argparse.Namespace) -> None:
    row = _series_row(a.collection, a.series_uid, a.uuid, a.root)
    sys.stdout.write(idc.viewer_url(row) + "\n")


# ------------------------------------------------------------------- parser


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scigantic-nci",
        description="GDC projects and IDC collections from the Scigantic mirrors, as CSV.",
    )
    p.add_argument("--version", action="version", version=f"scigantic-nci {__version__}")
    p.add_argument(
        "--root",
        default=None,
        help="local directory holding one project's or collection's files (the root= argument); "
        "otherwise $SCIGANTIC_NCI_ROOT, the notebook mount, or the public bucket is used",
    )
    sub = p.add_subparsers(dest="command", metavar="command")
    sub.required = True

    def out(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--out", default=None, metavar="FILE.csv", help="write CSV here instead of stdout")

    sp = sub.add_parser("projects", help="one row per mirrored GDC project")
    out(sp)
    sp.set_defaults(func=_cmd_projects)

    sp = sub.add_parser("collections", help="one row per mirrored IDC collection")
    out(sp)
    sp.set_defaults(func=_cmd_collections)

    sp = sub.add_parser(
        "build-index",
        help="write the catalog index file (PROJECTS.parquet or COLLECTIONS.parquet) for a mirror bucket",
        description="Build the projects() (gdc) or collections() (idc) frame from every project's or "
        "collection's BUILD_REPORT.json, ignoring any existing index, and write it as parquet. Upload the "
        f"result to s3://{GDC_BUCKET}/{INDEX_FILES[GDC_BUCKET]} or s3://{IDC_BUCKET}/{INDEX_FILES[IDC_BUCKET]} "
        "after every rebuild.",
    )
    sp.add_argument("mirror", choices=["gdc", "idc"])
    sp.add_argument("--out", required=True, metavar="PATH", help="parquet file to write")
    sp.set_defaults(func=_cmd_build_index)

    sp = sub.add_parser("tables", help="tables of a GDC project or IDC collection with rows, columns, bytes")
    sp.add_argument("id", metavar="PROJECT|collection")
    out(sp)
    sp.set_defaults(func=_cmd_tables)

    sp = sub.add_parser("readme", help="print the README.md of a project or collection")
    sp.add_argument("id", metavar="PROJECT|collection")
    sp.set_defaults(func=_cmd_readme)

    sp = sub.add_parser("clinical", help="GDC cases joined to the primary diagnosis, with survival columns")
    sp.add_argument("project", metavar="PROJECT")
    sp.add_argument("--all-diagnoses", action="store_true", help="keep every diagnosis row (primary_only=False)")
    out(sp)
    sp.set_defaults(func=_cmd_clinical)

    sp = sub.add_parser("expression", help="GDC STAR expression, genes x aliquots")
    sp.add_argument("project", metavar="PROJECT")
    sp.add_argument("--genes", required=True, metavar="TP53,ESR1", help="comma-separated gene symbols")
    sp.add_argument("--sample-type", default=None, metavar="TYPE", help='e.g. "Primary Tumor"; comma-separate several')
    sp.add_argument("--counts", action="store_true", help="raw STAR counts instead of TPM")
    out(sp)
    sp.set_defaults(func=_cmd_expression)

    sp = sub.add_parser("mutations", help="GDC masked somatic mutations, one row per call")
    sp.add_argument("project", metavar="PROJECT")
    sp.add_argument("--genes", default=None, metavar="TP53,RB1", help="comma-separated Hugo symbols")
    sp.add_argument("--non-silent", action="store_true", help="keep protein-changing classes only")
    out(sp)
    sp.set_defaults(func=_cmd_mutations)

    sp = sub.add_parser("mutation-frequency", help="fraction of tumor samples mutated per gene")
    sp.add_argument("project", metavar="PROJECT")
    sp.add_argument("--top", type=int, default=20, metavar="N", help="rows to print (default 20; 0 for all)")
    sp.add_argument("--include-silent", action="store_true", help="count silent calls too (non_silent=False)")
    out(sp)
    sp.set_defaults(func=_cmd_mutation_frequency)

    sp = sub.add_parser("series", help="IDC series.parquet, optionally one modality")
    sp.add_argument("collection")
    sp.add_argument("--modality", default=None, metavar="CT")
    sp.add_argument("--max-size-mb", type=float, default=None, metavar="MB", help="keep series at or below this size")
    out(sp)
    sp.set_defaults(func=_cmd_series)

    def series_id(sp: argparse.ArgumentParser) -> None:
        g = sp.add_mutually_exclusive_group(required=True)
        g.add_argument("--series-uid", default=None, metavar="UID", help="SeriesInstanceUID")
        g.add_argument("--uuid", default=None, metavar="UUID", help="crdc_series_uuid")

    sp = sub.add_parser("pull-series", help="download one series from the raw IDC bucket (anonymous)")
    sp.add_argument("collection")
    series_id(sp)
    sp.add_argument("--dest", required=True, metavar="DIR")
    sp.set_defaults(func=_cmd_pull_series)

    sp = sub.add_parser("sample", help="describe the sample series read_sample() would load")
    sp.add_argument("collection")
    sp.add_argument("--modality", default=None, metavar="CT")
    sp.set_defaults(func=_cmd_sample)

    sp = sub.add_parser("viewer-url", help="IDC web viewer link for one series")
    sp.add_argument("collection")
    series_id(sp)
    sp.set_defaults(func=_cmd_viewer_url)

    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if getattr(args, "top", None) == 0:
        args.top = None
    func: Any = args.func
    try:
        func(args)
    except (NciError, ValueError) as exc:
        sys.stderr.write(f"scigantic-nci: {exc}\n")
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
