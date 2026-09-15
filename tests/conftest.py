"""Shared pytest fixtures.

The suite reads a local copy of the two mirror buckets. The root is chosen in this order:

1. ``$SCIGANTIC_NCI_TEST_ROOT``: a directory holding ``gdc_mirror/`` and ``idc_mirror/``
   (and optionally ``gdcbuild_meth/``),
2. the full scratchpad mirrors used while building the buckets, when present on this machine,
3. ``tests/fixtures/mirror``, the subset every test reads, fetched from the public buckets by
   ``python tests/fixtures/fetch_mirror.py`` (this is what CI uses).

Set ``SCIGANTIC_NCI_TEST_FIXTURES_ONLY=1`` to force option 3 and reproduce a CI run locally.
Tests skip when no root has the data they need. Live S3 tests are marked ``network``.
"""
from __future__ import annotations

import os

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRATCH = "/private/tmp/claude-501/-Users-aaronkanzer-projects-scigantic/c826c167-1893-4b78-a1d1-5610f06f1463/scratchpad"
_FIXTURES = os.path.join(_HERE, "fixtures", "mirror")


def _pick_root() -> str:
    env = os.environ.get("SCIGANTIC_NCI_TEST_ROOT")
    if env:
        return env
    if os.environ.get("SCIGANTIC_NCI_TEST_FIXTURES_ONLY"):
        return _FIXTURES
    if os.path.isdir(os.path.join(_SCRATCH, "gdc_mirror")):
        return _SCRATCH
    return _FIXTURES


TEST_ROOT = _pick_root()
IDC_MIRROR = os.path.join(TEST_ROOT, "idc_mirror")
GDC_MIRROR = os.path.join(TEST_ROOT, "gdc_mirror")
# A TCGA-CHOL build that includes methylation/ (the public mirror was built without it).
GDC_METHYLATION_ROOT = os.path.join(TEST_ROOT, "gdcbuild_meth", "TCGA-CHOL")


def pytest_configure(config: pytest.Config) -> None:
    markers = [str(m) for m in config.getini("markers")]
    if not any(m.split(":")[0].strip() == "network" for m in markers):
        config.addinivalue_line("markers", "network: reads from public S3 (deselect with -m 'not network')")


def pytest_report_header(config: pytest.Config) -> str:
    return f"scigantic-nci test root: {TEST_ROOT} (gdc={os.path.isdir(GDC_MIRROR)}, idc={os.path.isdir(IDC_MIRROR)}, methylation={os.path.isdir(GDC_METHYLATION_ROOT)})"


@pytest.fixture
def idc_root(monkeypatch: pytest.MonkeyPatch) -> str:
    """Point the resolver at the local IDC mirror; skip when absent."""
    if not os.path.exists(os.path.join(IDC_MIRROR, "tcga_lihc", "BUILD_REPORT.json")):
        pytest.skip(f"local IDC mirror not present at {IDC_MIRROR} (run tests/fixtures/fetch_mirror.py)")
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", IDC_MIRROR)
    return IDC_MIRROR


@pytest.fixture
def gdc_root(monkeypatch: pytest.MonkeyPatch) -> str:
    """Point the resolver at the local GDC mirror; skip when absent."""
    if not os.path.exists(os.path.join(GDC_MIRROR, "TCGA-BLCA", "BUILD_REPORT.json")):
        pytest.skip(f"local GDC mirror not present at {GDC_MIRROR} (run tests/fixtures/fetch_mirror.py)")
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", GDC_MIRROR)
    return GDC_MIRROR


@pytest.fixture
def gdc_methylation_root() -> str:
    """Local TCGA-CHOL build with methylation tables, used through ``root=``; skip when absent."""
    if not os.path.exists(os.path.join(GDC_METHYLATION_ROOT, "methylation", "betas_human_methylation_450.parquet")):
        pytest.skip(f"TCGA-CHOL methylation build not present at {GDC_METHYLATION_ROOT}")
    return GDC_METHYLATION_ROOT
