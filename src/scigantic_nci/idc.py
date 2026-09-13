"""Imaging Data Commons (IDC) collections mirrored at s3://scigantic-idc-open/<collection_id>/.

Every collection has series.parquet, studies.parquet and patients.parquet (one row per
DICOM series, study and patient), a README.md and a BUILD_REPORT.json. Collections with
slide microscopy also have sm_series.parquet; 55 have analysis_results.parquet; 59 have
clinical/<table>.parquet with a clinical/dictionary.parquet. Under sample/ each collection
keeps one or more small DICOM series (a whole radiology series, or the thumbnail plus the
two lowest pyramid levels of one slide) so a notebook can look at pixels without pulling
anything from the raw IDC buckets.

Typical use::

    import scigantic_nci as nci

    nci.idc.collections()                                   # 176 rows, one per collection
    s = nci.idc.series("tcga_lihc", modality="CT")          # 777 rows
    vol = nci.idc.read_sample("tcga_lihc", "CT")            # Volume, array (36, 512, 512) float32
    files = nci.idc.pull_series(s.iloc[0], "/tmp/ct")       # anonymous read from idc-open-data

Everything reads through ``_store.resolve`` so a mounted archive, a local mirror
(``$SCIGANTIC_NCI_ROOT``) or the public bucket all behave the same. Only ``pull_series``
and ``series_size`` talk to the raw IDC buckets, anonymously through boto3.
"""
from __future__ import annotations

import importlib
import io
import math
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Callable

import numpy as np
import pandas as pd

from ._store import IDC_BUCKET, NciError, NotMirroredError, Store, list_keys, resolve

RAW_BUCKET_DEFAULT = "idc-open-data"
VIEWER_BASE = "https://viewer.imaging.datacommons.cancer.gov/viewer"

_SM_LIST_COLUMNS = [
    "embeddingMedium_CodeMeaning",
    "embeddingMedium_code_designator_value_str",
    "tissueFixative_CodeMeaning",
    "tissueFixative_code_designator_value_str",
    "staining_usingSubstance_CodeMeaning",
    "staining_usingSubstance_code_designator_value_str",
]

_ANALYSIS_COLUMNS = [
    "analysis_result_id",
    "analysis_result_title",
    "source_DOI",
    "source_url",
    "subjects",
    "collections",
    "modalities",
    "updated",
    "license_url",
    "license_long_name",
    "license_short_name",
    "description",
    "citation",
]


# ---------------------------------------------------------------------------
# collection ids and stores


def _local_root() -> str | None:
    return os.environ.get("SCIGANTIC_NCI_ROOT") or None


@lru_cache(maxsize=4)
def _collection_ids(local_root: str | None) -> tuple[str, ...]:
    """Collection ids known to this process: the local mirror when one is configured, else S3."""
    if local_root:
        for base in (os.path.join(local_root, IDC_BUCKET), local_root):
            if not os.path.isdir(base):
                continue
            ids = sorted(
                d for d in os.listdir(base) if os.path.exists(os.path.join(base, d, "BUILD_REPORT.json"))
            )
            if ids:
                return tuple(ids)
    return tuple(list_keys(IDC_BUCKET))


def collection_ids() -> list[str]:
    """Sorted ids of every mirrored collection (cached per process)."""
    return list(_collection_ids(_local_root()))


def _store(collection: str, root: str | None) -> Store:
    try:
        return resolve(IDC_BUCKET, collection, root)
    except NotMirroredError as e:
        if root:
            raise
        ids = collection_ids()
        raise NotMirroredError(
            f"{collection!r} is not a mirrored IDC collection ({e}). "
            f"Valid ids ({len(ids)}): {', '.join(ids)}"
        ) from None


# ---------------------------------------------------------------------------
# catalog


def report(collection: str, root: str | None = None) -> dict[str, Any]:
    """The collection's BUILD_REPORT.json as a dict (tables, samples, modalities, licenses, counts)."""
    d = _store(collection, root).read_json("BUILD_REPORT.json")
    return dict(d)


def readme(collection: str, root: str | None = None) -> str:
    """The collection's README.md (describes every table and its provenance)."""
    return _store(collection, root).read_bytes("README.md").decode("utf-8")


def tables(collection: str, root: str | None = None) -> pd.DataFrame:
    """Every table of the collection with rows, columns and bytes (from the build report).

    Same shape as ``gdc.tables``: columns table, rows, columns, bytes, sorted by table.
    """
    t = report(collection, root).get("tables", {})
    df = pd.DataFrame(
        [{"table": k, "rows": v.get("rows"), "columns": v.get("columns"), "bytes": v.get("bytes")} for k, v in t.items()],
        columns=["table", "rows", "columns", "bytes"],
    )
    return df.sort_values("table").reset_index(drop=True)


def _report_row(cid: str) -> dict[str, Any]:
    d = report(cid)
    return {
        "collection_id": d["collection_id"],
        "collection_name": d.get("collection_name"),
        "n_patients": d.get("n_patients"),
        "n_studies": d.get("n_studies"),
        "n_series": d.get("n_series"),
        "total_GB": d.get("total_GB"),
        "modalities": ", ".join(d.get("modalities", {})),
        "licenses": ", ".join(d.get("licenses", {})),
        "clinical_tables": len(d.get("clinical_tables", [])),
        "idc_data_version": d.get("idc_data_version"),
        "built_at": d.get("built_at"),
    }


_COLLECTIONS_CACHE: dict[str, pd.DataFrame] = {}


def collections() -> pd.DataFrame:
    """One row per mirrored collection, from each BUILD_REPORT.json (cached per process).

    Columns: collection_id, collection_name, n_patients, n_studies, n_series, total_GB,
    modalities (comma-joined, most series first), licenses, clinical_tables (count),
    idc_data_version, built_at. Reads 176 small JSON files with 16 threads.
    """
    key = _local_root() or ""
    cached = _COLLECTIONS_CACHE.get(key)
    if cached is None:
        ids = collection_ids()
        with ThreadPoolExecutor(max_workers=16) as pool:
            rows = list(pool.map(_report_row, ids))
        cached = pd.DataFrame(rows).sort_values("collection_id").reset_index(drop=True)
        _COLLECTIONS_CACHE[key] = cached
    return cached.copy()


# ---------------------------------------------------------------------------
# per-collection tables


def _modalities(rep: dict[str, Any]) -> list[str]:
    return list(rep.get("modalities", {}))


def _check_modality(rep: dict[str, Any], modality: str | None) -> None:
    if modality is None:
        return
    valid = _modalities(rep)
    if modality not in valid:
        raise ValueError(
            f"{rep['collection_id']} has no {modality!r} series; modalities: {', '.join(valid)}"
        )


def series(
    collection: str,
    modality: str | None = None,
    patient: str | list[str] | None = None,
    max_size_mb: float | None = None,
    columns: list[str] | None = None,
    root: str | None = None,
) -> pd.DataFrame:
    """One row per DICOM series (series.parquet), optionally filtered.

    ``modality`` must be one of the collection's modalities (ValueError lists them);
    ``patient`` is one PatientID or a list; ``max_size_mb`` keeps series at or below that
    size (series_size_MB). ``columns`` reads a subset of the 31 columns. Filters are pushed
    down to pyarrow so a filtered read of a 22,032-row collection does not build the full frame.
    """
    st = _store(collection, root)
    filters: list[tuple[str, str, Any]] = []
    if modality is not None:
        _check_modality(st.read_json("BUILD_REPORT.json"), modality)
        filters.append(("Modality", "==", modality))
    if patient is not None:
        pats = [patient] if isinstance(patient, str) else list(patient)
        filters.append(("PatientID", "in", pats))
    if max_size_mb is not None:
        filters.append(("series_size_MB", "<=", float(max_size_mb)))
    df = st.read_parquet("series.parquet", columns=columns, filters=filters or None)
    return df.reset_index(drop=True)


def studies(collection: str, root: str | None = None) -> pd.DataFrame:
    """One row per study: StudyInstanceUID, PatientID, StudyDate, StudyDescription, modalities,
    n_series, instances, size_MB, body_parts."""
    return _store(collection, root).read_parquet("studies.parquet")


def patients(collection: str, root: str | None = None) -> pd.DataFrame:
    """One row per patient: PatientID, PatientSex, PatientAge, n_studies, n_series, modalities, size_MB."""
    return _store(collection, root).read_parquet("patients.parquet")


def _join_list(v: Any) -> Any:
    if v is None:
        return None
    if isinstance(v, (list, tuple, np.ndarray)):
        return ", ".join(str(x) for x in v)
    if isinstance(v, float) and math.isnan(v):
        return None
    return v


def slides(collection: str, root: str | None = None) -> pd.DataFrame:
    """Slide microscopy series: sm_series.parquet joined to series.parquet on SeriesInstanceUID.

    The list-valued code columns (embedding medium, fixative, stain) are flattened to
    comma-joined strings, and a ``stain`` column holds the flattened
    staining_usingSubstance_CodeMeaning. Raises ValueError when the collection has no slides.
    """
    st = _store(collection, root)
    if not st.exists("sm_series.parquet"):
        rep = st.read_json("BUILD_REPORT.json")
        raise ValueError(
            f"{collection} has no slide microscopy (sm_series.parquet); modalities: "
            f"{', '.join(_modalities(rep))}"
        )
    sm = st.read_parquet("sm_series.parquet")
    for col in _SM_LIST_COLUMNS:
        if col in sm.columns:
            sm[col] = sm[col].map(_join_list)
    sm["stain"] = sm["staining_usingSubstance_CodeMeaning"]
    ser = st.read_parquet("series.parquet")
    ser = ser.drop(columns=[c for c in ser.columns if c in sm.columns and c != "SeriesInstanceUID"])
    out = sm.merge(ser, on="SeriesInstanceUID", how="left")
    return out.reset_index(drop=True)


def analysis_results(collection: str, root: str | None = None) -> pd.DataFrame:
    """Analysis results (segmentations, annotations) attached to the collection; empty frame when none."""
    st = _store(collection, root)
    if not st.exists("analysis_results.parquet"):
        return pd.DataFrame(columns=_ANALYSIS_COLUMNS)
    return st.read_parquet("analysis_results.parquet")


def clinical_tables(collection: str, root: str | None = None) -> list[str]:
    """Names of the collection's clinical tables (59 collections have at least one)."""
    rep = _store(collection, root).read_json("BUILD_REPORT.json")
    return [str(t) for t in rep.get("clinical_tables", [])]


def clinical(collection: str, table: str | None = None, root: str | None = None) -> pd.DataFrame:
    """One clinical table (clinical/<table>.parquet); default is the first listed in BUILD_REPORT."""
    st = _store(collection, root)
    tables = [str(t) for t in st.read_json("BUILD_REPORT.json").get("clinical_tables", [])]
    if not tables:
        raise ValueError(f"{collection} has no clinical tables")
    if table is None:
        table = tables[0]
    if table not in tables:
        raise ValueError(f"{collection} has no clinical table {table!r}; tables: {', '.join(tables)}")
    return st.read_parquet(f"clinical/{table}.parquet")


def clinical_dictionary(collection: str, root: str | None = None) -> pd.DataFrame:
    """clinical/dictionary.parquet: short_table_name, column, column_label, values (one row per column)."""
    st = _store(collection, root)
    if not st.exists("clinical/dictionary.parquet"):
        raise ValueError(f"{collection} has no clinical tables")
    return st.read_parquet("clinical/dictionary.parquet")


# ---------------------------------------------------------------------------
# samples


def samples(collection: str, root: str | None = None) -> pd.DataFrame:
    """The sample series kept under sample/ (from BUILD_REPORT samples), one row each.

    Columns: modality, kind (series or sm_levels), files, bytes, instanceCount_full_series,
    series_size_MB_full, PatientID, StudyInstanceUID, SeriesInstanceUID, series_aws_url,
    license, source_DOI, folder. Empty for the one collection without samples.
    """
    rep = _store(collection, root).read_json("BUILD_REPORT.json")
    cols = [
        "modality",
        "kind",
        "files",
        "bytes",
        "instanceCount_full_series",
        "series_size_MB_full",
        "PatientID",
        "StudyInstanceUID",
        "SeriesInstanceUID",
        "series_aws_url",
        "license",
        "source_DOI",
        "folder",
    ]
    rows = [{c: s.get(c) for c in cols} for s in rep.get("samples", [])]
    return pd.DataFrame(rows, columns=cols)


def _pick_sample(rep: dict[str, Any], modality: str | None) -> dict[str, Any]:
    samp = list(rep.get("samples", []))
    if not samp:
        raise ValueError(f"{rep['collection_id']} has no sample series")
    if modality is None:
        return dict(samp[0])
    for s in samp:
        if s.get("modality") == modality:
            return dict(s)
    have = ", ".join(str(s.get("modality")) for s in samp)
    raise ValueError(f"{rep['collection_id']} has no {modality!r} sample; sample modalities: {have}")


def sample_files(collection: str, modality: str | None = None, root: str | None = None) -> list[str]:
    """Relative paths (under the collection prefix) of the sample .dcm files, sorted.

    With ``modality`` only that sample's folder is listed.
    """
    st = _store(collection, root)
    if modality is None:
        base = "sample"
    else:
        base = _pick_sample(st.read_json("BUILD_REPORT.json"), modality)["folder"]
    return sorted(p for p in st.walk_files(base) if p.endswith(".dcm"))


def _pydicom() -> Any:
    try:
        import pydicom
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise ImportError(
            "pydicom is required to read DICOM samples; install it with "
            'pip install "scigantic-nci[dicom]"'
        ) from e
    return pydicom


def _dcmread(data: bytes) -> Any:
    return _pydicom().dcmread(io.BytesIO(data))


def _pixels(ds: Any, index: int | None = None) -> Any:
    """Decoded pixel array of a dataset, or one frame of a multi-frame instance.

    Turns the pydicom "no decoder plugin" RuntimeError into an NciError with an install hint.
    """
    try:
        if index is not None:
            try:
                px = importlib.import_module("pydicom.pixels")
            except ImportError:  # pydicom < 3 decodes every frame
                return np.asarray(ds.pixel_array[index])
            return np.asarray(px.pixel_array(ds, index=index))
        return np.asarray(ds.pixel_array)
    except (RuntimeError, NotImplementedError) as e:
        ts = str(getattr(getattr(ds, "file_meta", None), "TransferSyntaxUID", ""))
        raise NciError(
            f"cannot decode pixel data (transfer syntax {ts}): {e}. Compressed DICOM needs a "
            'decoder plugin: pip install "pylibjpeg[all]" (or pillow for JPEG baseline).'
        ) from e


def _um_per_px(ds: Any) -> float | None:
    """Micrometres per pixel from SharedFunctionalGroupsSequence/PixelMeasures, else None."""
    try:
        ps = ds.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing
        return float(ps[0]) * 1000.0
    except (AttributeError, IndexError, TypeError, ValueError):
        return None


def _pixel_spacing(ds: Any) -> tuple[float | None, float | None]:
    ps = ds.get("PixelSpacing")
    if ps is None:
        try:
            ps = ds.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing
        except (AttributeError, IndexError, TypeError):
            return (None, None)
    try:
        return (float(ps[0]), float(ps[1]))
    except (IndexError, TypeError, ValueError):
        return (None, None)


@dataclass
class Volume:
    """A radiology sample series decoded to an array.

    ``array`` is float32 with shape (n_instances, rows, cols[, channels]); instances are
    sorted by ImagePositionPatient[2] then InstanceNumber, CT values are rescaled to
    Hounsfield units, and a multi-frame instance contributes its middle frame. When the
    instances do not share one shape ``array`` is None, ``ragged`` is True and ``images``
    holds the individual 2D arrays in the same order. ``spacing`` is (slice, row, column)
    in mm with None where unknown; ``datasets`` are the sorted pydicom datasets.
    """

    modality: str
    description: str
    array: np.ndarray[Any, Any] | None = field(repr=False)
    spacing: tuple[float | None, float | None, float | None]
    datasets: list[Any] = field(repr=False)
    images: list[np.ndarray[Any, Any]] | None = field(default=None, repr=False)
    ragged: bool = False

    @property
    def shape(self) -> tuple[int, ...] | None:
        return None if self.array is None else tuple(self.array.shape)

    def __repr__(self) -> str:
        if self.ragged and self.images is not None:
            shapes = sorted({im.shape for im in self.images})
            body = f"ragged, {len(self.images)} images with shapes {shapes}"
        else:
            body = f"shape={self.shape}"
        return (
            f"Volume(modality={self.modality!r}, description={self.description!r}, "
            f"{body}, spacing={self.spacing}, n_instances={len(self.datasets)})"
        )


@dataclass
class SlideLevel:
    """One instance of a DICOM WSI pyramid: its size, ImageType and where it lives."""

    index: int
    width: int
    height: int
    image_type: str
    rows: int
    cols: int
    n_frames: int
    um_per_px: float | None
    path: str


@dataclass
class SlideLevels:
    """The pyramid levels of a sample slide (sm_levels sample), lowest resolution first.

    ``levels`` lists each instance's width (TotalPixelMatrixColumns), ImageType (THUMBNAIL or
    VOLUME), tile geometry and micrometres per pixel. ``level(i)`` assembles that level's
    tiles into one array (see ``assemble_level``).
    """

    modality: str
    description: str
    levels: list[SlideLevel]
    datasets: list[Any] = field(repr=False)
    _reader: Callable[[int], bytes] = field(repr=False)

    def level(self, i: int) -> np.ndarray[Any, Any]:
        """Assemble level ``i`` (index into ``levels``) into a (rows, cols[, channels]) array."""
        if not 0 <= i < len(self.levels):
            raise IndexError(f"level {i} out of range; {len(self.levels)} levels")
        return assemble_level(self._reader(i))

    @property
    def widths(self) -> list[int]:
        return [lv.width for lv in self.levels]


def assemble_level(path: str | bytes | Any) -> np.ndarray[Any, Any]:
    """Assemble the tiled frames of one WSI instance into a full level.

    ``path`` is a local .dcm path, the file's bytes, or an already-read pydicom Dataset.
    Tiles are row-major with ``ncol = ceil(TotalPixelMatrixColumns / Columns)``; the canvas has
    shape ``(TotalPixelMatrixRows, TotalPixelMatrixColumns) + frames.shape[3:]`` so RGB slides
    keep their channel axis and fluorescence slides (SamplesPerPixel 1) have none. A level with
    NumberOfFrames == 1 is one tile and is returned as is.
    """
    if isinstance(path, (bytes, bytearray)):
        ds = _dcmread(bytes(path))
    elif isinstance(path, str):
        ds = _pydicom().dcmread(path)
    else:
        ds = path
    frames = _pixels(ds)
    n_frames = int(ds.get("NumberOfFrames", 1) or 1)
    rows, cols = int(ds.Rows), int(ds.Columns)
    if n_frames == 1:
        tile: np.ndarray[Any, Any] = np.asarray(frames)
        if tile.ndim >= 3 and tile.shape[0] == 1 and tile.shape[1:3] == (rows, cols):
            tile = tile[0]
        return tile
    if frames.ndim < 3:
        frames = frames[np.newaxis]
    total_rows = int(ds.get("TotalPixelMatrixRows", rows))
    total_cols = int(ds.get("TotalPixelMatrixColumns", cols))
    ncol = math.ceil(total_cols / cols)
    nrow = math.ceil(total_rows / rows)
    canvas = np.zeros((nrow * rows, ncol * cols) + frames.shape[3:], dtype=frames.dtype)
    for k in range(min(n_frames, frames.shape[0])):
        r, c = divmod(k, ncol)
        canvas[r * rows : (r + 1) * rows, c * cols : (c + 1) * cols] = frames[k]
    return canvas[:total_rows, :total_cols]


def _sort_key(ds: Any) -> tuple[float, int]:
    z = 0.0
    ipp = ds.get("ImagePositionPatient")
    if ipp is not None:
        try:
            z = float(ipp[2])
        except (IndexError, TypeError, ValueError):
            z = 0.0
    inst = ds.get("InstanceNumber")
    try:
        n = int(inst) if inst is not None else 0
    except (TypeError, ValueError):
        n = 0
    return (z, n)


def _read_volume(sample: dict[str, Any], datasets: list[Any]) -> Volume:
    datasets = sorted(datasets, key=_sort_key)
    modality = str(sample.get("modality") or datasets[0].get("Modality", ""))
    description = str(datasets[0].get("SeriesDescription", "") or "")
    images: list[np.ndarray[Any, Any]] = []
    for ds in datasets:
        n_frames = int(ds.get("NumberOfFrames", 1) or 1)
        img = _pixels(ds, index=n_frames // 2) if n_frames > 1 else _pixels(ds)
        arr = np.asarray(img, dtype=np.float32)
        if modality == "CT":
            slope = float(ds.get("RescaleSlope", 1) or 1)
            intercept = float(ds.get("RescaleIntercept", 0) or 0)
            arr = arr * slope + intercept
        images.append(arr)
    dy, dx = _pixel_spacing(datasets[0])
    dz: float | None = None
    zs = [_sort_key(ds)[0] for ds in datasets if ds.get("ImagePositionPatient") is not None]
    if len(zs) > 1:
        diffs = np.diff(np.unique(np.asarray(zs)))
        if diffs.size:
            dz = float(np.median(diffs))
    if dz is None:
        for tag in ("SpacingBetweenSlices", "SliceThickness"):
            v = datasets[0].get(tag)
            if v is not None:
                try:
                    dz = float(v)
                    break
                except (TypeError, ValueError):
                    pass
    shapes = {im.shape for im in images}
    if len(shapes) == 1:
        return Volume(modality, description, np.stack(images), (dz, dy, dx), datasets)
    return Volume(modality, description, None, (dz, dy, dx), datasets, images=images, ragged=True)


def _read_slide(sample: dict[str, Any], st: Store, rels: list[str], datasets: list[Any]) -> SlideLevels:
    order = sorted(range(len(datasets)), key=lambda i: int(datasets[i].get("TotalPixelMatrixColumns", 0) or 0))
    datasets = [datasets[i] for i in order]
    rels = [rels[i] for i in order]
    levels = []
    for i, (ds, rel) in enumerate(zip(datasets, rels)):
        it = ds.get("ImageType")
        image_type = str(it[2]) if it is not None and len(it) > 2 else ""
        levels.append(
            SlideLevel(
                index=i,
                width=int(ds.get("TotalPixelMatrixColumns", ds.Columns)),
                height=int(ds.get("TotalPixelMatrixRows", ds.Rows)),
                image_type=image_type,
                rows=int(ds.Rows),
                cols=int(ds.Columns),
                n_frames=int(ds.get("NumberOfFrames", 1) or 1),
                um_per_px=_um_per_px(ds),
                path=os.path.join(st.local, rel) if st.local else f"s3://{st.bucket}/{st.key}/{rel}",
            )
        )
    description = str(datasets[0].get("SeriesDescription", "") or "")

    def reader(i: int) -> bytes:
        return st.read_bytes(rels[i])

    return SlideLevels(str(sample.get("modality", "SM")), description, levels, datasets, reader)


def read_sample(collection: str, modality: str | None = None, root: str | None = None) -> Volume | SlideLevels:
    """Load one sample series with pydicom (needs the ``dicom`` extra).

    Radiology samples (kind ``series``) return a ``Volume``; slide samples (kind ``sm_levels``)
    return a ``SlideLevels`` whose ``level(i)`` assembles a pyramid level. ``modality`` picks
    which sample when a collection has several (ValueError lists them); default is the first.
    """
    st = _store(collection, root)
    rep = st.read_json("BUILD_REPORT.json")
    sample = _pick_sample(rep, modality)
    rels = sorted(p for p in st.walk_files(sample["folder"]) if p.endswith(".dcm"))
    if not rels:
        raise NciError(f"{collection}: no .dcm files under {sample['folder']}")
    datasets = [_dcmread(st.read_bytes(r)) for r in rels]
    if sample.get("kind") == "sm_levels" or sample.get("modality") == "SM":
        return _read_slide(sample, st, rels, datasets)
    return _read_volume(sample, datasets)


# ---------------------------------------------------------------------------
# raw IDC objects


def _raw_client() -> Any:
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    return boto3.client("s3", region_name="us-east-1", config=Config(signature_version=UNSIGNED))


def _bucket_and_uuid(row_or_uuid: Any, bucket: str | None) -> tuple[str, str]:
    if isinstance(row_or_uuid, str):
        return (bucket or RAW_BUCKET_DEFAULT, row_or_uuid)
    try:
        uuid = str(row_or_uuid["crdc_series_uuid"])
    except (KeyError, TypeError, IndexError):
        raise ValueError(
            "expected a series.parquet row (with crdc_series_uuid and aws_bucket) or a crdc_series_uuid string"
        ) from None
    if bucket is None:
        b = row_or_uuid.get("aws_bucket") if hasattr(row_or_uuid, "get") else None
        bucket = str(b) if b else RAW_BUCKET_DEFAULT
    return (bucket, uuid)


def _list_raw(client: Any, bucket: str, uuid: str) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    token: str | None = None
    while True:
        kw: dict[str, Any] = {"Bucket": bucket, "Prefix": f"{uuid}/"}
        if token:
            kw["ContinuationToken"] = token
        resp = client.list_objects_v2(**kw)
        out.extend(o for o in resp.get("Contents", []) if not str(o["Key"]).endswith("/"))
        if not resp.get("IsTruncated"):
            return out
        token = resp.get("NextContinuationToken")


def pull_series(row_or_uuid: Any, dest: str, bucket: str | None = None) -> list[str]:
    """Download one series from the raw IDC bucket into ``dest`` (anonymous, 8 threads).

    ``row_or_uuid`` is a series.parquet row (aws_bucket and crdc_series_uuid are used) or a
    crdc_series_uuid string (bucket defaults to idc-open-data). Objects live at
    ``s3://<bucket>/<uuid>/<instance>.dcm``; the folder-marker key ending in '/' is skipped.
    Returns the sorted local file paths.
    """
    bucket, uuid = _bucket_and_uuid(row_or_uuid, bucket)
    client = _raw_client()
    objs = _list_raw(client, bucket, uuid)
    if not objs:
        raise NciError(f"no objects under s3://{bucket}/{uuid}/")
    os.makedirs(dest, exist_ok=True)

    def fetch(key: str) -> str:
        local = os.path.join(dest, os.path.basename(key))
        client.download_file(bucket, key, local)
        return local

    with ThreadPoolExecutor(max_workers=8) as pool:
        files = list(pool.map(fetch, [str(o["Key"]) for o in objs]))
    return sorted(files)


def series_size(row_or_uuid: Any, bucket: str | None = None) -> float:
    """Size of one series in MB: series_size_MB from a series.parquet row, or summed from S3 for a uuid."""
    if not isinstance(row_or_uuid, str):
        try:
            return float(row_or_uuid["series_size_MB"])
        except (KeyError, TypeError, IndexError, ValueError):
            pass
    b, uuid = _bucket_and_uuid(row_or_uuid, bucket)
    objs = _list_raw(_raw_client(), b, uuid)
    return float(sum(int(o["Size"]) for o in objs)) / 1e6


def viewer_url(row: Any) -> str:
    """Link to the IDC web viewer for a series.parquet row (StudyInstanceUID, SeriesInstanceUID)."""
    try:
        study = str(row["StudyInstanceUID"])
        ser = str(row["SeriesInstanceUID"])
    except (KeyError, TypeError, IndexError):
        raise ValueError("expected a row with StudyInstanceUID and SeriesInstanceUID") from None
    return f"{VIEWER_BASE}/{study}?seriesInstanceUID={ser}"


__all__ = [
    "SlideLevel",
    "SlideLevels",
    "Volume",
    "analysis_results",
    "assemble_level",
    "clinical",
    "clinical_dictionary",
    "clinical_tables",
    "collection_ids",
    "collections",
    "patients",
    "pull_series",
    "read_sample",
    "readme",
    "report",
    "sample_files",
    "samples",
    "series",
    "series_size",
    "slides",
    "studies",
    "tables",
    "viewer_url",
]
