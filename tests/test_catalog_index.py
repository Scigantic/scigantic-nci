"""The catalog index files (PROJECTS.parquet, COLLECTIONS.parquet) and ``build-index``.

Each test lays out a directory of symlinks to the local mirror's project or collection
directories (the mirrors themselves are read-only), writes an index next to them with the
console script, and checks that ``projects()`` / ``collections()`` read it and return the
same frame as the build from every BUILD_REPORT.json.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Callable, Iterator

import pandas as pd
import pyarrow.fs as pafs
import pytest

from scigantic_nci import _store, gdc, idc
from scigantic_nci._store import GDC_BUCKET, IDC_BUCKET, list_keys
from scigantic_nci.cli import main

BLCA = "TCGA-BLCA"
LIHC = "tcga_lihc"


def _clear_caches() -> None:
    gdc._project_ids.cache_clear()
    gdc._projects_frame.cache_clear()
    idc._collection_ids.cache_clear()
    idc._COLLECTIONS_CACHE.clear()


@pytest.fixture(autouse=True)
def fresh_caches() -> Iterator[None]:
    """Some tests here reuse the ``None`` cache key (no local root), so never leak frames."""
    _clear_caches()
    yield
    _clear_caches()


def _symlink_mirror(mirror: str, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    for name in sorted(os.listdir(mirror)):
        if os.path.isfile(os.path.join(mirror, name, "BUILD_REPORT.json")):
            os.symlink(os.path.join(mirror, name), dest / name)


def _layout(tmp_path: Path, mirror: str, bucket: str, flat: bool) -> tuple[str, Path]:
    """``root/<key>`` (flat) or ``root/<bucket>/<key>``; returns the root and the directory holding the ids."""
    root = tmp_path / "root"
    base = root if flat else root / bucket
    _symlink_mirror(mirror, base)
    return str(root), base


def _spy(monkeypatch: pytest.MonkeyPatch, module: Any, name: str) -> list[str]:
    """Record the ids a module's row builder is called for."""
    built: list[str] = []
    real: Callable[[str], dict[str, Any]] = getattr(module, name)

    def wrapper(key: str) -> dict[str, Any]:
        built.append(key)
        return real(key)

    monkeypatch.setattr(module, name, wrapper)
    return built


@pytest.fixture
def gdc_expected(gdc_root: str, monkeypatch: pytest.MonkeyPatch) -> pd.DataFrame:
    """projects() built from every BUILD_REPORT.json of the local GDC mirror."""
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", gdc_root)
    return gdc._build_projects(use_index=False)


@pytest.fixture
def idc_expected(idc_root: str, monkeypatch: pytest.MonkeyPatch) -> pd.DataFrame:
    """collections() built from every BUILD_REPORT.json of the local IDC mirror."""
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", idc_root)
    return idc._build_collections(use_index=False)


# --------------------------------------------------------------------- build-index


@pytest.mark.parametrize("flat", [True, False])
def test_projects_reads_local_index(
    gdc_root: str, gdc_expected: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], flat: bool
) -> None:
    root, base = _layout(tmp_path, gdc_root, GDC_BUCKET, flat)
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", root)
    out = base / "PROJECTS.parquet"
    assert main(["build-index", "gdc", "--out", str(out)]) == 0
    assert capsys.readouterr().out == f"{out}: 57 rows, {out.stat().st_size} bytes\n"
    built = _spy(monkeypatch, gdc, "_project_row")
    p = gdc.projects()
    assert built == []  # every row came from the index
    pd.testing.assert_frame_equal(p, gdc_expected)
    assert int(p.set_index("project").loc[BLCA, "cases"]) == 412
    # the index file is not a project
    assert "PROJECTS.parquet" not in gdc._known_projects()


def test_collections_reads_local_index(
    idc_root: str, idc_expected: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    root, base = _layout(tmp_path, idc_root, IDC_BUCKET, flat=False)
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", root)
    out = base / "COLLECTIONS.parquet"
    assert main(["build-index", "idc", "--out", str(out)]) == 0
    assert capsys.readouterr().out == f"{out}: 176 rows, {out.stat().st_size} bytes\n"
    built = _spy(monkeypatch, idc, "_report_row")
    c = idc.collections()
    assert built == []
    pd.testing.assert_frame_equal(c, idc_expected)
    assert len(idc.collection_ids()) == 176 and "COLLECTIONS.parquet" not in idc.collection_ids()


def test_build_index_ignores_existing_index(
    gdc_root: str, idc_root: str, gdc_expected: pd.DataFrame, idc_expected: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    for mirror, bucket, kind, expected, key in (
        (gdc_root, GDC_BUCKET, "gdc", gdc_expected, "project"),
        (idc_root, IDC_BUCKET, "idc", idc_expected, "collection_id"),
    ):
        root, base = _layout(tmp_path / kind, mirror, bucket, flat=True)
        monkeypatch.setenv("SCIGANTIC_NCI_ROOT", root)
        stale = expected.copy()
        stale["built_at"] = "stale"
        stale.to_parquet(base / _store.INDEX_FILES[bucket], index=False)
        out = tmp_path / f"{kind}.parquet"
        assert main(["build-index", kind, "--out", str(out)]) == 0
        assert capsys.readouterr().out == f"{out}: {len(expected)} rows, {out.stat().st_size} bytes\n"
        written = pd.read_parquet(out)
        assert "stale" not in set(written["built_at"])
        pd.testing.assert_frame_equal(written, expected)
        assert written[key].is_unique
        assert not os.path.exists(str(out) + ".part")


def test_write_index_creates_missing_directories(tmp_path: Path) -> None:
    out = tmp_path / "not" / "yet" / "PROJECTS.parquet"
    df = pd.DataFrame({"project": ["TCGA-X"], "cases": [3]})
    size = _store.write_index(df, str(out))
    assert out.stat().st_size == size and not (out.parent / "PROJECTS.parquet.part").exists()
    pd.testing.assert_frame_equal(pd.read_parquet(out), df)


def test_build_index_usage_errors(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        main(["build-index", "gdc"])  # --out is required
    with pytest.raises(SystemExit):
        main(["build-index", "tcga", "--out", "x.parquet"])
    assert "invalid choice" in capsys.readouterr().err


# ------------------------------------------------------------ partial and bad indexes


def test_projects_merges_ids_missing_from_index(
    gdc_root: str, gdc_expected: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, base = _layout(tmp_path, gdc_root, GDC_BUCKET, flat=True)
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", root)
    index = gdc_expected[gdc_expected["project"] != BLCA].copy()
    index.loc[index["project"] == "TCGA-ACC", "built_at"] = "from-index"
    index["extra_column"] = 1  # a newer index may carry more columns; they are ignored
    index.to_parquet(base / "PROJECTS.parquet", index=False)
    built = _spy(monkeypatch, gdc, "_project_row")
    p = gdc.projects()
    assert built == [BLCA]
    assert list(p.columns) == list(gdc_expected.columns)
    assert p.set_index("project").loc["TCGA-ACC", "built_at"] == "from-index"
    want = gdc_expected.copy()
    want.loc[want["project"] == "TCGA-ACC", "built_at"] = "from-index"
    pd.testing.assert_frame_equal(p, want)


def test_collections_merges_ids_missing_from_index(
    idc_root: str, idc_expected: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, base = _layout(tmp_path, idc_root, IDC_BUCKET, flat=False)
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", root)
    idc_expected[~idc_expected["collection_id"].isin([LIHC, "nlst"])].to_parquet(base / "COLLECTIONS.parquet", index=False)
    built = _spy(monkeypatch, idc, "_report_row")
    c = idc.collections()
    assert sorted(built) == ["nlst", LIHC]
    pd.testing.assert_frame_equal(c, idc_expected)


def _drop_column(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=["built_at"])


def _rename_column(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={df.columns[4]: "renamed"})


def _duplicate_id(df: pd.DataFrame) -> pd.DataFrame:
    return pd.concat([df, df.iloc[:1]], ignore_index=True)


@pytest.mark.parametrize("corrupt", [_drop_column, _rename_column, _duplicate_id, None], ids=["dropped", "renamed", "duplicate-id", "not-parquet"])
def test_bad_index_falls_back_entirely(
    gdc_root: str,
    idc_root: str,
    gdc_expected: pd.DataFrame,
    idc_expected: pd.DataFrame,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    corrupt: Callable[[pd.DataFrame], pd.DataFrame] | None,
) -> None:
    for mirror, bucket, expected, module, builder, frame in (
        (gdc_root, GDC_BUCKET, gdc_expected, gdc, "_project_row", gdc.projects),
        (idc_root, IDC_BUCKET, idc_expected, idc, "_report_row", idc.collections),
    ):
        root, base = _layout(tmp_path / bucket, mirror, bucket, flat=False)
        monkeypatch.setenv("SCIGANTIC_NCI_ROOT", root)
        path = base / _store.INDEX_FILES[bucket]
        if corrupt is None:
            path.write_bytes(b"PAR1 this is not a parquet file")
        else:
            stale = expected.copy()
            stale["built_at"] = "stale"
            corrupt(stale).to_parquet(path, index=False)
        built = _spy(monkeypatch, module, builder)
        got = frame()
        assert len(built) == len(expected)
        pd.testing.assert_frame_equal(got, expected)


# ------------------------------------------------------------- the bucket-root index


def test_bucket_index_through_s3_filesystem(
    gdc_root: str, gdc_expected: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The S3 path, with pyarrow's S3 filesystem swapped for a local tree laid out like the bucket."""
    bucket_dir = tmp_path / "s3" / GDC_BUCKET
    _symlink_mirror(gdc_root, bucket_dir)
    fs = pafs.SubTreeFileSystem(str(tmp_path / "s3"), pafs.LocalFileSystem())
    monkeypatch.setattr(_store, "_s3fs", lambda: fs)
    monkeypatch.delenv("SCIGANTIC_NCI_ROOT", raising=False)
    monkeypatch.setenv("SCIGANTIC_MOUNT_PATH", str(tmp_path / "no-mount"))

    # no index yet: every row is built from the "bucket"
    built = _spy(monkeypatch, gdc, "_project_row")
    pd.testing.assert_frame_equal(gdc.projects(), gdc_expected)
    assert len(built) == 57
    _clear_caches()

    # a root file is not a key: list_keys (what 0.1.1 enumerates ids with) still sees 57 projects
    index = gdc_expected[gdc_expected["project"] != "TCGA-CHOL"].copy()
    index.loc[index["project"] == "TCGA-ACC", "built_at"] = "from-index"
    index.to_parquet(bucket_dir / "PROJECTS.parquet", index=False)
    assert list_keys(GDC_BUCKET) == gdc_expected["project"].tolist()

    # with the index: only the project missing from it is built
    built.clear()
    p = gdc.projects()
    assert built == ["TCGA-CHOL"]
    assert p.set_index("project").loc["TCGA-ACC", "built_at"] == "from-index"
    _clear_caches()

    # a project on the notebook mount is read from the mount, never from the bucket's index
    monkeypatch.setenv("SCIGANTIC_MOUNT_PATH", os.path.join(gdc_root, "TCGA-ACC"))
    built.clear()
    p = gdc.projects()
    assert sorted(built) == ["TCGA-ACC", "TCGA-CHOL"]
    pd.testing.assert_frame_equal(p, gdc_expected)


def test_local_mirror_never_reads_bucket_index(
    gdc_root: str, gdc_expected: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Projects read from $SCIGANTIC_NCI_ROOT without a local index are built from their files."""
    calls: list[str | None] = []
    real = _store.read_index

    def spy(bucket: str, columns: Any, key_column: str, local_dir: str | None = None) -> dict[str, dict[str, Any]]:
        calls.append(local_dir)
        return real(bucket, columns, key_column, local_dir)

    monkeypatch.setattr(_store, "read_index", spy)
    root, _ = _layout(tmp_path, gdc_root, GDC_BUCKET, flat=True)
    monkeypatch.setenv("SCIGANTIC_NCI_ROOT", root)
    pd.testing.assert_frame_equal(gdc.projects(), gdc_expected)
    assert calls == []
