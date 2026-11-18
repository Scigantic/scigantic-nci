"""idc behaviors that need only a BUILD_REPORT.json, so they run without the mirror fixture."""
from __future__ import annotations

import pytest

from scigantic_nci import idc


def _tiny_collection(tmp_path, report: dict) -> str:  # type: ignore[no-untyped-def]
    import json
    import os

    d = tmp_path / report["collection_id"]  # root= is the collection directory itself
    os.makedirs(d)
    (d / "BUILD_REPORT.json").write_text(json.dumps(report))
    return str(d)


def test_sample_files_is_empty_when_the_collection_has_no_samples(tmp_path) -> None:  # type: ignore[no-untyped-def]
    root = _tiny_collection(tmp_path, {"collection_id": "nosamples", "samples": [], "clinical_tables": []})
    assert idc.sample_files("nosamples", root=root) == []


def test_clinical_says_where_a_table_that_is_not_mirrored_can_be_read(tmp_path) -> None:  # type: ignore[no-untyped-def]
    root = _tiny_collection(
        tmp_path,
        {"collection_id": "nc", "samples": [], "clinical_tables": [], "clinical_tables_not_mirrored": {"tables": ["nc_clinical"]}},
    )
    assert idc.clinical_tables("nc", root=root) == []
    with pytest.raises(ValueError, match="get_clinical_table\\('nc_clinical'\\)"):
        idc.clinical("nc", root=root)


def test_folders_beside_the_mirror_are_not_collections() -> None:
    from scigantic_nci._store import is_data_prefix

    assert is_data_prefix("tcga_lihc") and is_data_prefix("4d_lung") and is_data_prefix("TCGA-BLCA")
    for name in ("_ask", "_thumbnails", ".aws-datasync", ""):
        assert not is_data_prefix(name)


@pytest.mark.network
def test_live_collection_listing_has_only_collections() -> None:
    ids = idc.collection_ids()
    assert ids and not [i for i in ids if i.startswith(("_", "."))]
    assert len(idc.collections()) == len(ids)
