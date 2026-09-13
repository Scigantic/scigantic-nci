"""Shared pytest fixtures.

Each module's tests run against a complete local mirror when one is present (the
scratchpad copies used while building the mirrors); they skip otherwise. Live S3 tests are
marked ``network``.
"""
from __future__ import annotations

import os

import pytest

_SCRATCH = "/private/tmp/claude-501/-Users-aaronkanzer-projects-scigantic/c826c167-1893-4b78-a1d1-5610f06f1463/scratchpad"
IDC_MIRROR = os.path.join(_SCRATCH, "idc_mirror")
GDC_MIRROR = os.path.join(_SCRATCH, "gdc_mirror")
# A TCGA-CHOL build that includes methylation/ (the main mirror was built without it).
GDC_METHYLATION_ROOT = os.path.join(_SCRATCH, "gdcbuild_meth", "TCGA-CHOL")


def pytest_configure(config: pytest.Config) -> None:
    markers = [str(m) for m in config.getini("markers")]
    if not any(m.split(":")[0].strip() == "network" for m in markers):
        config.addinivalue_line("markers", "network: reads from public S3 (deselect with -m 'not network')")


@pytest.fixture
def idc_root(monkeypatch: pytest.MonkeyPatch) -> str:
    """Point the resolver at the local IDC mirror (all 176 collections); skip when absent."""
    if not os.path.exists(os.path.join(IDC_MIRROR, "tcga_lihc", "BUILD_REPORT.json")):
        pytest.skip(f"local IDC mirror not present at {IDC_MIRROR}")
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", IDC_MIRROR)
    return IDC_MIRROR


@pytest.fixture
def gdc_root(monkeypatch: pytest.MonkeyPatch) -> str:
    """Point the resolver at the local GDC mirror (TCGA-ACC, -BLCA, -BRCA, -CHOL, ...); skip when absent."""
    if not os.path.exists(os.path.join(GDC_MIRROR, "TCGA-BLCA", "BUILD_REPORT.json")):
        pytest.skip(f"local GDC mirror not present at {GDC_MIRROR}")
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", GDC_MIRROR)
    return GDC_MIRROR


@pytest.fixture
def gdc_methylation_root() -> str:
    """Local TCGA-CHOL build with methylation tables, used through ``root=``; skip when absent."""
    if not os.path.exists(os.path.join(GDC_METHYLATION_ROOT, "methylation", "betas_human_methylation_450.parquet")):
        pytest.skip(f"TCGA-CHOL methylation build not present at {GDC_METHYLATION_ROOT}")
    return GDC_METHYLATION_ROOT
