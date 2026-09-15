"""Where the bytes come from.

Every GDC project and IDC collection lives at one prefix of a public bucket:

    s3://scigantic-gdc-open/<PROJECT>/       e.g. TCGA-BRCA/
    s3://scigantic-idc-open/<collection_id>/ e.g. tcga_brca/

A `Store` resolves one such prefix to a place it can read parquet and files
from, in this order:

1. an explicit ``root`` argument (a local directory holding that prefix's files),
2. ``$SCIGANTIC_NCI_ROOT/<bucket>/<key>`` if that directory exists (local mirrors, tests),
3. the Scigantic notebook mount at ``$SCIGANTIC_MOUNT_PATH`` (default ``/mnt/archive``)
   when its ``BUILD_REPORT.json`` names the same project / collection,
4. the public bucket, read anonymously through pyarrow's S3 filesystem.

Anonymous S3 works because both buckets are public-read; no credentials are
ever required, and none are used even if present (so a pod role or a laptop
profile cannot change what you get).

Each bucket also has one catalog index file at its root, ``PROJECTS.parquet`` in
the GDC bucket and ``COLLECTIONS.parquet`` in the IDC bucket, holding exactly the
``gdc.projects()`` / ``idc.collections()`` frame. It is a file, not a prefix, so
``list_keys`` (directories only) never mistakes it for an id. `catalog_rows` reads
it in place of one BUILD_REPORT.json per id; ids always come from listing.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Callable, Iterator, Sequence

import pandas as pd
import pyarrow.parquet as pq

GDC_BUCKET = "scigantic-gdc-open"
IDC_BUCKET = "scigantic-idc-open"
REGION = "us-east-1"
#: Catalog index file at the root of each bucket (see `catalog_rows`).
INDEX_FILES = {GDC_BUCKET: "PROJECTS.parquet", IDC_BUCKET: "COLLECTIONS.parquet"}


class NciError(Exception):
    """Base class for scigantic-nci errors."""


class NotMirroredError(NciError):
    """The project / collection has no prefix in the mirror bucket."""


@lru_cache(maxsize=1)
def _s3fs() -> Any:
    import pyarrow.fs as pafs

    return pafs.S3FileSystem(anonymous=True, region=REGION)


@dataclass(frozen=True)
class Store:
    bucket: str
    key: str  # prefix without trailing slash, e.g. "TCGA-BRCA" or "tcga_brca"
    local: str | None  # local directory when reading from disk, else None

    @property
    def location(self) -> str:
        return self.local if self.local else f"s3://{self.bucket}/{self.key}/"

    def _path(self, rel: str) -> str:
        return os.path.join(self.local, rel) if self.local else f"{self.bucket}/{self.key}/{rel}"

    def exists(self, rel: str) -> bool:
        if self.local:
            return os.path.exists(os.path.join(self.local, rel))
        import pyarrow.fs as pafs

        info = _s3fs().get_file_info(self._path(rel))
        return bool(info.type != pafs.FileType.NotFound)

    def read_parquet(self, rel: str, columns: list[str] | None = None, filters: Any = None) -> pd.DataFrame:
        kw: dict[str, Any] = {}
        if columns is not None:
            kw["columns"] = columns
        if filters is not None:
            kw["filters"] = filters
        if self.local:
            return pd.DataFrame(pq.read_table(os.path.join(self.local, rel), **kw).to_pandas())
        return pd.DataFrame(pq.read_table(self._path(rel), filesystem=_s3fs(), **kw).to_pandas())

    def parquet_schema(self, rel: str) -> list[str]:
        if self.local:
            return list(pq.read_schema(os.path.join(self.local, rel)).names)
        return list(pq.read_schema(self._path(rel), filesystem=_s3fs()).names)

    def read_bytes(self, rel: str) -> bytes:
        if self.local:
            with open(os.path.join(self.local, rel), "rb") as fh:
                return fh.read()
        with _s3fs().open_input_stream(self._path(rel)) as f:
            return bytes(f.read())

    def read_json(self, rel: str) -> Any:
        return json.loads(self.read_bytes(rel).decode("utf-8"))

    def listdir(self, rel: str = "") -> list[str]:
        """Names directly under ``rel`` (files and directories), no recursion."""
        if self.local:
            p = os.path.join(self.local, rel) if rel else self.local
            return sorted(os.listdir(p)) if os.path.isdir(p) else []
        import pyarrow.fs as pafs

        base = self._path(rel).rstrip("/")
        infos = _s3fs().get_file_info(pafs.FileSelector(base, recursive=False))
        return sorted(os.path.basename(i.path.rstrip("/")) for i in infos)

    def walk_files(self, rel: str = "") -> Iterator[str]:
        """Every file path under ``rel``, relative to the store root."""
        if self.local:
            base = os.path.join(self.local, rel) if rel else self.local
            for dp, _, files in os.walk(base):
                for f in sorted(files):
                    yield os.path.relpath(os.path.join(dp, f), self.local)
            return
        import pyarrow.fs as pafs

        base = self._path(rel).rstrip("/")
        for i in _s3fs().get_file_info(pafs.FileSelector(base, recursive=True)):
            if i.type == pafs.FileType.File:
                yield os.path.relpath(i.path, f"{self.bucket}/{self.key}")


def _env_base(bucket: str, key: str) -> str | None:
    """The ``$SCIGANTIC_NCI_ROOT/<bucket>`` or ``$SCIGANTIC_NCI_ROOT`` directory holding ``<key>/BUILD_REPORT.json``."""
    env_root = os.environ.get("SCIGANTIC_NCI_ROOT")
    if env_root:
        for base in (os.path.join(env_root, bucket), env_root):
            if os.path.exists(os.path.join(base, key, "BUILD_REPORT.json")):
                return base
    return None


def _mount() -> str:
    return os.environ.get("SCIGANTIC_MOUNT_PATH", "/mnt/archive")


def _mount_key(mount: str, bucket: str) -> str | None:
    """The project or collection id named by the mount's BUILD_REPORT.json, if any."""
    rep = os.path.join(mount, "BUILD_REPORT.json")
    if not os.path.exists(rep):
        return None
    try:
        with open(rep) as fh:
            d = json.load(fh)
    except (OSError, ValueError):
        return None
    field = "project" if bucket == GDC_BUCKET else "collection_id"
    value = d.get(field) if isinstance(d, dict) else None
    return value if isinstance(value, str) else None


def _mount_matches(mount: str, bucket: str, key: str) -> bool:
    return _mount_key(mount, bucket) == key


def resolve(bucket: str, key: str, root: str | None = None) -> Store:
    """Find where ``<bucket>/<key>/`` can be read from (see module docstring)."""
    if root:
        if not os.path.exists(os.path.join(root, "BUILD_REPORT.json")):
            raise NotMirroredError(f"{root} has no BUILD_REPORT.json; not a scigantic-nci archive")
        return Store(bucket, key, root)
    base = _env_base(bucket, key)
    if base is not None:
        return Store(bucket, key, os.path.join(base, key))
    mount = _mount()
    if _mount_matches(mount, bucket, key):
        return Store(bucket, key, mount)
    store = Store(bucket, key, None)
    if not store.exists("BUILD_REPORT.json"):
        raise NotMirroredError(f"s3://{bucket}/{key}/ is not in the mirror (no BUILD_REPORT.json)")
    return store


def list_keys(bucket: str) -> list[str]:
    """Top-level prefixes of a mirror bucket (project ids or collection ids)."""
    import pyarrow.fs as pafs

    infos = _s3fs().get_file_info(pafs.FileSelector(bucket, recursive=False))
    return sorted(os.path.basename(i.path.rstrip("/")) for i in infos if i.type == pafs.FileType.Directory)


def read_index(bucket: str, columns: Sequence[str], key_column: str, local_dir: str | None = None) -> dict[str, dict[str, Any]]:
    """Rows of a bucket's catalog index keyed by id, or ``{}`` when it cannot be used.

    Reads ``<local_dir>/<index file>`` when ``local_dir`` is given, else the file at the
    bucket root, anonymously. Absent, unreadable, missing any of ``columns`` or with a
    null or repeated id: ``{}``, so the caller builds every row itself. Extra columns are
    ignored.
    """
    name = INDEX_FILES[bucket]
    try:
        if local_dir is not None:
            df = pq.read_table(os.path.join(local_dir, name)).to_pandas()
        else:
            with _s3fs().open_input_file(f"{bucket}/{name}") as fh:
                df = pq.read_table(fh).to_pandas()
    except Exception:  # the index is an optimisation; any failure means build the rows
        return {}
    if not set(columns) <= set(df.columns) or df[key_column].isna().any() or not df[key_column].is_unique:
        return {}
    return {str(r[key_column]): r for r in df[list(columns)].to_dict("records")}


def write_index(df: pd.DataFrame, path: str) -> int:
    """Write a catalog frame as an index file; returns its size in bytes.

    Plain parquet with no pandas metadata and ``string`` rather than ``large_string``
    columns, so every supported pandas and pyarrow reads it back with the dtypes it would
    give the frame built from the BUILD_REPORT.json files. Written to ``<path>.part`` and
    renamed, so a failed build never leaves a truncated index.
    """
    import pyarrow as pa

    table = pa.Table.from_pandas(df, preserve_index=False)
    schema = pa.schema([pa.field(f.name, pa.string()) if pa.types.is_large_string(f.type) else f for f in table.schema])
    table = table.cast(schema).replace_schema_metadata(None)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".part"
    pq.write_table(table, tmp)
    os.replace(tmp, path)
    return os.path.getsize(path)


def _local_index_dir(bucket: str) -> str | None:
    """``$SCIGANTIC_NCI_ROOT/<bucket>`` or ``$SCIGANTIC_NCI_ROOT`` when it holds the bucket's index file."""
    env_root = os.environ.get("SCIGANTIC_NCI_ROOT")
    if env_root:
        for base in (os.path.join(env_root, bucket), env_root):
            if os.path.isfile(os.path.join(base, INDEX_FILES[bucket])):
                return base
    return None


def catalog_rows(
    bucket: str,
    ids: Sequence[str],
    build_row: Callable[[str], dict[str, Any]],
    columns: Sequence[str],
    key_column: str,
    use_index: bool = True,
    workers: int = 16,
) -> list[dict[str, Any]]:
    """One catalog row per id, in the order of ``ids``.

    With ``use_index`` one index is consulted first: a local copy of the index file under
    ``$SCIGANTIC_NCI_ROOT`` when there is one (used for every id not on the notebook
    mount), otherwise the file at the bucket root, used only for the ids that resolve to
    the bucket (ids read from a local directory or the mount never trigger a network read
    and are never shadowed by the bucket's copy). Every id the index does not cover is
    built by ``build_row`` from its BUILD_REPORT.json, ``workers`` at a time.
    ``use_index=False`` builds every row, which is how the index itself is written.
    """
    rows: dict[str, dict[str, Any]] = {}
    if use_index and ids:
        mounted = _mount_key(_mount(), bucket)
        local_dir = _local_index_dir(bucket)
        candidates: list[str] = []
        for k in ids:
            if _env_base(bucket, k) is not None:  # read from a local directory
                if local_dir is not None:
                    candidates.append(k)
            elif k != mounted:  # read from the bucket
                candidates.append(k)
        if candidates:
            index = read_index(bucket, columns, key_column, local_dir)
            rows.update((k, index[k]) for k in candidates if k in index)
    missing = [k for k in ids if k not in rows]
    if missing:
        with ThreadPoolExecutor(max_workers=max(1, min(workers, len(missing)))) as pool:
            rows.update(zip(missing, pool.map(build_row, missing)))
    return [rows[k] for k in ids]
