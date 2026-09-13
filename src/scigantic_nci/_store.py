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
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Iterator

import pandas as pd
import pyarrow.parquet as pq

GDC_BUCKET = "scigantic-gdc-open"
IDC_BUCKET = "scigantic-idc-open"
REGION = "us-east-1"


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


def _mount_matches(mount: str, bucket: str, key: str) -> bool:
    rep = os.path.join(mount, "BUILD_REPORT.json")
    if not os.path.exists(rep):
        return False
    try:
        with open(rep) as fh:
            d = json.load(fh)
    except (OSError, ValueError):
        return False
    field = "project" if bucket == GDC_BUCKET else "collection_id"
    return bool(d.get(field) == key)


def resolve(bucket: str, key: str, root: str | None = None) -> Store:
    """Find where ``<bucket>/<key>/`` can be read from (see module docstring)."""
    if root:
        if not os.path.exists(os.path.join(root, "BUILD_REPORT.json")):
            raise NotMirroredError(f"{root} has no BUILD_REPORT.json; not a scigantic-nci archive")
        return Store(bucket, key, root)
    env_root = os.environ.get("SCIGANTIC_NCI_ROOT")
    if env_root:
        cand = os.path.join(env_root, bucket, key)
        if os.path.exists(os.path.join(cand, "BUILD_REPORT.json")):
            return Store(bucket, key, cand)
        cand = os.path.join(env_root, key)
        if os.path.exists(os.path.join(cand, "BUILD_REPORT.json")):
            return Store(bucket, key, cand)
    mount = os.environ.get("SCIGANTIC_MOUNT_PATH", "/mnt/archive")
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
