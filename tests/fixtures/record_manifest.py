"""pytest plugin that records every mirror file the suite reads and rewrites manifest.txt.

Run against the full local mirrors (the scratchpad copies or SCIGANTIC_NCI_TEST_ROOT):

    python -m pytest -q -p tests.fixtures.record_manifest tests

The store reads parquet through pyarrow (not Python's open), so pq.read_table and
pq.read_schema are wrapped as well as the open() audit hook; os.listdir and os.scandir
(which glob and os.walk use) are wrapped too, and every regular file directly inside a
listed directory is kept so that listing-based functions (copy_number_workflows,
sample_files) see the same entries as on the full mirror. Files outside the two public
buckets (the local-only methylation build) are ignored.
"""
from __future__ import annotations

import os
import sys
from typing import Any

import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "manifest.txt")
_seen: set[str] = set()
_roots: dict[str, str] = {}


def _rec(p: Any) -> None:
    if isinstance(p, (str, os.PathLike)):
        p = os.fspath(p)
        if os.path.isfile(p):
            _seen.add(os.path.abspath(p))


def _rec_dir(p: Any) -> None:
    if isinstance(p, (str, bytes, os.PathLike)) and not isinstance(p, bytes):
        p = os.fspath(p)
        if os.path.isdir(p):
            for name in _listdir(p):
                _rec(os.path.join(p, name))


_rt, _rs = pq.read_table, pq.read_schema


def _read_table(p: Any, *a: Any, **k: Any) -> Any:
    _rec(p)
    return _rt(p, *a, **k)


def _read_schema(p: Any, *a: Any, **k: Any) -> Any:
    _rec(p)
    return _rs(p, *a, **k)


pq.read_table, pq.read_schema = _read_table, _read_schema
_listdir, _scandir = os.listdir, os.scandir


def _listdir_rec(p: Any = ".") -> Any:
    _rec_dir(p)
    return _listdir(p)


def _scandir_rec(p: Any = ".") -> Any:
    _rec_dir(p)
    return _scandir(p)


os.listdir, os.scandir = _listdir_rec, _scandir_rec
sys.addaudithook(lambda ev, args: _rec(args[0]) if ev == "open" else None)


def pytest_configure(config: Any) -> None:
    from tests import conftest

    _roots.update({conftest.GDC_MIRROR: "scigantic-gdc-open", conftest.IDC_MIRROR: "scigantic-idc-open"})


def pytest_sessionfinish(session: Any, exitstatus: int) -> None:
    rows = []
    for p in sorted(_seen):
        for root, bucket in _roots.items():
            root = os.path.abspath(root) + os.sep
            if p.startswith(root):
                rows.append(f"{bucket} {p[len(root):]} {os.path.getsize(p)}")
    with open(MANIFEST, "w") as fh:
        fh.write("# bucket key bytes: every object the test suite reads, recorded by tests/fixtures/record_manifest.py\n")
        fh.write("\n".join(sorted(rows)) + "\n")
    print(f"\nrecorded {len(rows)} objects to {MANIFEST}")
