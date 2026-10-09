"""Tests for the idc-index bridge in scigantic_nci.idc, against a stand-in IDCClient (no network).

The real client's calls were checked by hand against idc-index 0.13.0 (IDC data v25) on 2026-10-09.
"""
from __future__ import annotations

import importlib
from typing import Any

import pandas as pd
import pytest

from scigantic_nci import NciError, idc

class FakeClient:
    CITATION_FORMAT_BIBTEX = "bibtex-fmt"
    CITATION_FORMAT_APA = "apa-fmt"
    CITATION_FORMAT_JSON = "json-fmt"
    CITATION_FORMAT_TURTLE = "turtle-fmt"

    def __init__(self) -> None:
        self.index = pd.DataFrame(
            {
                "collection_id": ["tcga_lihc", "tcga_lihc", "tcga_lihc", "brand_new"],
                "Modality": ["CT", "CT", "MR", "CT"],
                "PatientID": ["p1", "p2", "p1", "p9"],
                "series_size_MB": [10.0, 200.0, 5.0, 1.0],
                "crdc_series_uuid": ["u1", "u2", "u3", "u4"],
                "license_short_name": ["CC BY 4.0", "CC BY-NC 4.0", "CC BY 3.0", "CC BY 4.0"],
            }
        )
        self.downloads: list[dict[str, Any]] = []

    def get_idc_version(self) -> str:
        return "v25"

    def sql_query(self, sql: str) -> pd.DataFrame:
        return pd.DataFrame({"sql": [sql]})

    def download_from_selection(self, **kw: Any) -> None:
        self.downloads.append(kw)

    def citations_from_selection(self, collection_id: str, citation_format: str) -> list[str]:
        return [f"{collection_id}:{citation_format}"]


@pytest.fixture()
def fake(monkeypatch: pytest.MonkeyPatch) -> FakeClient:
    fc = FakeClient()
    monkeypatch.setattr(idc, "client", lambda: fc)
    return fc


def test_client_without_extra_names_the_install(monkeypatch: pytest.MonkeyPatch) -> None:
    def no_idc_index(name: str) -> Any:
        raise ImportError(name)

    idc.client.cache_clear()
    monkeypatch.setattr(importlib, "import_module", no_idc_index)
    with pytest.raises(NciError, match=r"scigantic-nci\[idc\]"):
        idc.client()


def test_freshness_reports_versions_and_new_collections(fake: FakeClient, monkeypatch: pytest.MonkeyPatch) -> None:
    mirrored = pd.DataFrame(
        {"collection_id": ["tcga_lihc", "gone"], "n_series": [3, 7], "idc_data_version": ["v24", "v24"]}
    )
    monkeypatch.setattr(idc, "collections", lambda: mirrored)
    f = idc.freshness()
    assert f["mirror_version"] == "v24"
    assert f["live_version"] == "v25"
    assert f["stale"] is True
    assert f["live_series"] == 4
    assert f["mirror_series"] == 10
    assert f["collections_not_mirrored"] == ["brand_new"]
    assert f["collections_dropped_upstream"] == ["gone"]


def test_live_series_filters_like_the_mirror(fake: FakeClient) -> None:
    ct = idc.series("tcga_lihc", modality="CT", live=True)
    assert list(ct["crdc_series_uuid"]) == ["u1", "u2"]
    small = idc.series("tcga_lihc", max_size_mb=10, live=True)
    assert sorted(small["crdc_series_uuid"]) == ["u1", "u3"]
    p1 = idc.series("tcga_lihc", patient="p1", columns=["crdc_series_uuid"], live=True)
    assert list(p1.columns) == ["crdc_series_uuid"] and len(p1) == 2
    with pytest.raises(ValueError, match="modalities: CT, MR"):
        idc.series("tcga_lihc", modality="PT", live=True)
    with pytest.raises(NciError, match="not a collection"):
        idc.series("nope", live=True)
    assert len(idc.series("brand_new", live=True)) == 1


def test_query_passes_sql_through(fake: FakeClient) -> None:
    assert idc.query("select 1")["sql"][0] == "select 1"


def test_download_accepts_frames_rows_and_uuids(fake: FakeClient, tmp_path: Any) -> None:
    frame = idc.series("tcga_lihc", live=True)
    assert idc.download(frame, str(tmp_path / "a")) == ["u1", "u2", "u3"]
    assert idc.download(frame.iloc[0], str(tmp_path / "b")) == ["u1"]
    assert idc.download(["u3", "u4"], str(tmp_path / "c"), dry_run=True) == ["u3", "u4"]
    assert fake.downloads[-1]["crdc_series_uuid"] == ["u3", "u4"]
    assert fake.downloads[-1]["dry_run"] is True


def test_download_can_drop_noncommercial_series(fake: FakeClient, tmp_path: Any) -> None:
    frame = idc.series("tcga_lihc", live=True)
    assert idc.download(frame, str(tmp_path), exclude_noncommercial=True) == ["u1", "u3"]
    assert fake.downloads[-1]["crdc_series_uuid"] == ["u1", "u3"]
    with pytest.raises(NciError, match="CC BY-NC"):
        idc.download("u2", str(tmp_path), exclude_noncommercial=True)
    with pytest.raises(NciError, match="empty"):
        idc.download([], str(tmp_path))


def test_citations_formats(fake: FakeClient) -> None:
    assert idc.citations("tcga_lihc") == ["tcga_lihc:bibtex-fmt"]
    assert idc.citations("tcga_lihc", "apa") == ["tcga_lihc:apa-fmt"]
    with pytest.raises(ValueError, match="format must be"):
        idc.citations("tcga_lihc", "ris")
