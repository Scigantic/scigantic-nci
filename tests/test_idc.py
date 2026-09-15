"""Tests for scigantic_nci.idc against the local IDC mirror (plus two live-S3 tests).

Every number asserted here was measured on the mirror built 2026-09-13 (IDC data v24).
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from scigantic_nci import NciError, NotMirroredError, idc

pytestmark = pytest.mark.usefixtures("idc_root")


# ---------------------------------------------------------------------------
# catalog


def test_collections_catalog() -> None:
    c = idc.collections()
    assert len(c) == 176
    assert list(c.columns) == [
        "collection_id",
        "collection_name",
        "n_patients",
        "n_studies",
        "n_series",
        "total_GB",
        "modalities",
        "licenses",
        "clinical_tables",
        "idc_data_version",
        "built_at",
    ]
    lihc = c.set_index("collection_id").loc["tcga_lihc"]
    assert lihc["collection_name"] == "TCGA-LIHC"
    assert lihc["n_series"] == 3457
    assert lihc["n_patients"] == 377
    assert lihc["n_studies"] == 614
    assert lihc["modalities"] == "MR, SEG, SM, CT, PT"
    assert lihc["clinical_tables"] == 2
    assert lihc["idc_data_version"] == "v24"
    assert int(c["n_series"].sum()) == 1_032_911
    assert int(c["n_patients"].sum()) == 85_442
    assert c.sort_values("n_series").iloc[-1]["collection_id"] == "nlst"
    assert (c["clinical_tables"] > 0).sum() == 59
    # cached per process: same content, but a copy
    c2 = idc.collections()
    assert c2 is not c and c2.equals(c)


def test_collection_ids_and_unknown_collection() -> None:
    ids = idc.collection_ids()
    assert len(ids) == 176 and ids == sorted(ids)
    with pytest.raises(NotMirroredError) as ei:
        idc.series("not_a_collection")
    msg = str(ei.value)
    assert "Valid ids (176)" in msg and "tcga_lihc" in msg and "acrin_6698" in msg


def test_report_and_readme() -> None:
    rep = idc.report("tcga_lihc")
    assert rep["collection_id"] == "tcga_lihc"
    assert rep["modalities"] == {"MR": 910, "SEG": 899, "SM": 870, "CT": 777, "PT": 1}
    assert rep["tables"]["series"]["rows"] == 3457
    txt = idc.readme("tcga_lihc")
    assert "series.parquet" in txt


# ---------------------------------------------------------------------------
# tcga_lihc tables


def test_series_full_and_filters() -> None:
    s = idc.series("tcga_lihc")
    assert s.shape == (3457, 31)
    assert s["PatientID"].nunique() == 377
    assert s["Modality"].value_counts().to_dict() == {"MR": 910, "SEG": 899, "SM": 870, "CT": 777, "PT": 1}
    ct = idc.series("tcga_lihc", modality="CT")
    assert len(ct) == 777 and set(ct["Modality"]) == {"CT"}
    small = idc.series("tcga_lihc", modality="CT", max_size_mb=1)
    assert len(small) == 96 and small["series_size_MB"].max() <= 1
    one = idc.series("tcga_lihc", patient="TCGA-DD-A3A6")
    assert one["Modality"].value_counts().to_dict() == {"SM": 3, "CT": 2, "SEG": 2}
    two = idc.series("tcga_lihc", patient=["TCGA-DD-A3A6", "TCGA-DD-A4NG"])
    assert two["PatientID"].nunique() == 2
    sub = idc.series("tcga_lihc", modality="PT", columns=["SeriesInstanceUID", "instanceCount"])
    assert list(sub.columns) == ["SeriesInstanceUID", "instanceCount"] and len(sub) == 1
    with pytest.raises(ValueError, match="modalities: MR, SEG, SM, CT, PT"):
        idc.series("tcga_lihc", modality="XX")


def test_studies_patients() -> None:
    st = idc.studies("tcga_lihc")
    assert st.shape == (614, 9)
    assert st["PatientID"].nunique() == 377
    p = idc.patients("tcga_lihc")
    assert p.shape == (377, 7)
    assert p["n_series"].sum() == 3457


def test_slides_join_and_stain_flattening() -> None:
    sl = idc.slides("tcga_lihc")
    assert len(sl) == 870
    assert set(sl["Modality"]) == {"SM"}
    assert sl["stain"].tolist() == ["hematoxylin stain, water soluble eosin stain"] * 870
    assert sl["staining_usingSubstance_CodeMeaning"].iloc[0] == "hematoxylin stain, water soluble eosin stain"
    assert isinstance(sl["embeddingMedium_CodeMeaning"].iloc[0], str)
    assert sl["max_TotalPixelMatrixColumns"].max() == 209439
    assert sl["series_size_MB"].notna().all()
    with pytest.raises(ValueError, match="no slide microscopy"):
        idc.slides("acrin_6698")


def test_analysis_results() -> None:
    ar = idc.analysis_results("tcga_lihc")
    assert ar["analysis_result_id"].tolist() == ["bamf_aimi_annotations", "tcga_sbu_til_maps"]
    assert ar["subjects"].tolist() == [4226, 7600]
    empty = idc.analysis_results("tcga_chol")
    assert empty.empty and len(empty.columns) == 13


def test_clinical_lihc() -> None:
    assert idc.clinical_tables("tcga_lihc") == [
        "bamf_aimi_annotations_liver_ct_qa_results",
        "bamf_aimi_annotations_liver_mr_qa_results",
    ]
    assert idc.clinical("tcga_lihc").shape == (98, 15)
    assert idc.clinical("tcga_lihc", "bamf_aimi_annotations_liver_mr_qa_results").shape == (72, 15)
    d = idc.clinical_dictionary("tcga_lihc")
    assert d.shape == (30, 4) and list(d.columns) == ["short_table_name", "column", "column_label", "values"]
    with pytest.raises(ValueError, match="tables: "):
        idc.clinical("tcga_lihc", "nope")
    with pytest.raises(ValueError, match="no clinical tables"):
        idc.clinical("tcga_chol")
    with pytest.raises(ValueError, match="no clinical tables"):
        idc.clinical_dictionary("tcga_chol")


# ---------------------------------------------------------------------------
# acrin_6698


def test_acrin_6698_tables() -> None:
    assert idc.patients("acrin_6698").shape == (385, 7)
    s = idc.series("acrin_6698", columns=["Modality", "PatientID"])
    assert len(s) == 18747
    assert s["Modality"].value_counts().to_dict() == {"MR": 16534, "SEG": 2213}
    assert idc.clinical_tables("acrin_6698") == ["acrin_6698_clinical"]
    clin = idc.clinical("acrin_6698")
    assert clin.shape == (385, 32)
    assert clin["tcia_patient_id"].nunique() == 385
    assert idc.clinical_dictionary("acrin_6698").shape == (32, 4)


# ---------------------------------------------------------------------------
# samples


def test_samples_and_sample_files() -> None:
    sm = idc.samples("tcga_lihc")
    assert sm["modality"].tolist() == ["CT", "MR", "PT"]
    assert sm["kind"].tolist() == ["series"] * 3
    assert sm["files"].tolist() == [36, 111, 215]
    ct_files = idc.sample_files("tcga_lihc", "CT")
    assert len(ct_files) == 36 and all(f.endswith(".dcm") for f in ct_files)
    assert ct_files[0].startswith("sample/CT_1.3.6.1.4.1.14519.5.2.1.3344.4008.1590978269182606855431778873.8/")
    assert len(idc.sample_files("tcga_lihc")) == 36 + 111 + 215
    with pytest.raises(ValueError, match="sample modalities: CT, MR, PT"):
        idc.sample_files("tcga_lihc", "SM")
    # the one collection without any sample series
    assert idc.samples("b_mode_and_ceus_liver").empty
    with pytest.raises(ValueError, match="no sample series"):
        idc.read_sample("b_mode_and_ceus_liver")


def test_samples_license_short_name() -> None:
    lihc = idc.samples("tcga_lihc")
    cols = list(lihc.columns)
    assert cols[cols.index("license") + 1] == "license_short_name"
    assert lihc["license_short_name"].tolist() == ["CC BY 3.0"] * 3
    assert lihc["license_short_name"].equals(lihc["license"].rename("license_short_name"))
    # a collection that mixes CC BY 4.0 and CC BY-NC 3.0 series: the mirror never copies a
    # CC BY-NC series, so its one sample is a CC BY 4.0 SEG and its CT/RTSTRUCT are not sampled
    nsclc = idc.samples("nsclc_radiomics")
    assert nsclc["modality"].tolist() == ["SEG"]
    assert nsclc["license_short_name"].tolist() == ["CC BY 4.0"]
    assert set(idc.report("nsclc_radiomics")["licenses"]) == {"CC BY 4.0", "CC BY-NC 3.0"}
    ser = idc.series("nsclc_radiomics", columns=["SeriesInstanceUID", "license_short_name"])
    lic = ser.set_index("SeriesInstanceUID")["license_short_name"]
    assert nsclc["SeriesInstanceUID"].map(lic).tolist() == nsclc["license_short_name"].tolist()
    assert "license_short_name" in idc.samples("b_mode_and_ceus_liver").columns
    # an all-CC BY-NC collection keeps its index tables and has no sample at all
    assert idc.samples("phantom_fda").empty
    assert set(idc.report("phantom_fda")["licenses"]) == {"CC BY-NC 3.0"}


def test_read_sample_ct_volume() -> None:
    pytest.importorskip("pydicom")
    v = idc.read_sample("tcga_lihc", "CT")
    assert isinstance(v, idc.Volume)
    assert v.modality == "CT" and v.description == "AXIAL"
    assert v.array is not None and v.array.shape == (36, 512, 512) and v.array.dtype == np.float32
    assert not v.ragged and v.shape == (36, 512, 512)
    # Hounsfield: air is -1000 after rescale (this scanner stored RescaleIntercept 0, so raw == HU)
    assert float(v.array.min()) == -1000.0 and float(v.array.max()) == 1323.0
    assert v.spacing == (5.0, 0.585938, 0.585938)
    zs = [float(ds.ImagePositionPatient[2]) for ds in v.datasets]
    assert zs == sorted(zs) and len(v.datasets) == 36
    assert "shape=(36, 512, 512)" in repr(v)


def test_read_sample_mr_default_is_first_sample() -> None:
    pytest.importorskip("pydicom")
    v = idc.read_sample("tcga_lihc")
    assert isinstance(v, idc.Volume) and v.modality == "CT"
    mr = idc.read_sample("tcga_lihc", "MR")
    assert isinstance(mr, idc.Volume)
    assert mr.shape == (111, 144, 192) and mr.spacing == (6.5, 1.9791666269302, 1.9791666269302)


def test_tcga_chol_slide_levels() -> None:
    pytest.importorskip("pydicom")
    assert idc.report("tcga_chol")["modalities"] == {"SM": 110}
    assert idc.samples("tcga_chol")["kind"].tolist() == ["sm_levels"]
    sl = idc.read_sample("tcga_chol")
    assert isinstance(sl, idc.SlideLevels)
    assert sl.description == "Frozen HE TP TSA"
    assert [lv.image_type for lv in sl.levels] == ["THUMBNAIL", "VOLUME", "VOLUME"]
    assert sl.widths == [708, 1992, 7968]
    assert [lv.height for lv in sl.levels] == [768, 2158, 8633]
    assert [lv.n_frames for lv in sl.levels] == [1, 81, 1224]
    assert [round(lv.um_per_px or 0, 2) for lv in sl.levels] == [11.37, 4.04, 1.01]
    assert all(os.path.exists(lv.path) for lv in sl.levels)
    thumb = sl.level(0)
    assert thumb.shape == (768, 708, 3) and thumb.dtype == np.uint8
    with pytest.raises(IndexError):
        sl.level(3)
    try:
        lvl1 = sl.level(1)
    except NciError as e:  # JPEG baseline tiles need a decoder plugin
        assert "pylibjpeg" in str(e)
        pytest.skip("no JPEG decoder plugin installed")
    assert lvl1.shape == (2158, 1992, 3)


def test_htan_ohsu_fluorescence_no_channel_axis() -> None:
    pytest.importorskip("pydicom")
    sl = idc.read_sample("htan_ohsu")
    assert isinstance(sl, idc.SlideLevels)
    assert sl.description == "mIHC"
    assert len(sl.levels) == 3
    assert all(lv.image_type == "VOLUME" and lv.n_frames == 1 for lv in sl.levels)
    assert all((lv.height, lv.width) == (103, 143) for lv in sl.levels)
    assert all(ds.SamplesPerPixel == 1 for ds in sl.datasets)
    img = sl.level(0)
    assert img.shape == (103, 143) and img.ndim == 2
    assert round(sl.levels[0].um_per_px or 0, 3) == 8.067


def test_multiframe_instance() -> None:
    # remind's US sample is one instance holding a 39-frame cine loop; read_sample takes the middle frame
    pydicom = pytest.importorskip("pydicom")
    files = idc.sample_files("remind", "US")
    assert len(files) == 1
    ds = pydicom.dcmread(os.path.join(_mirror(), "remind", files[0]))
    assert int(ds.NumberOfFrames) == 39 and (ds.Rows, ds.Columns) == (599, 725)
    v = idc.read_sample("remind", "US")
    assert isinstance(v, idc.Volume) and v.shape == (1, 599, 725)
    assert np.array_equal(v.array[0], ds.pixel_array[19].astype(v.array.dtype))


def _mirror() -> str:
    return os.environ["SCIGANTIC_NCI_ROOT"]


class _FakeTiled:
    """Minimal stand-in for a tiled WSI Dataset so the assembly is tested without a decoder."""

    def __init__(self, rows: int, cols: int, total_rows: int, total_cols: int, channels: int) -> None:
        self.Rows, self.Columns = rows, cols
        self.TotalPixelMatrixRows, self.TotalPixelMatrixColumns = total_rows, total_cols
        ncol = -(-total_cols // cols)
        nrow = -(-total_rows // rows)
        self.NumberOfFrames = nrow * ncol
        shape = (self.NumberOfFrames, rows, cols) + ((channels,) if channels else ())
        # frame k is filled with the value k so every tile's placement is checkable
        self.pixel_array = np.arange(self.NumberOfFrames, dtype=np.uint8).reshape((-1,) + (1,) * (len(shape) - 1)) * np.ones(shape, np.uint8)

    def get(self, key: str, default: object = None) -> object:
        return getattr(self, key, default)


def test_assemble_level_row_major_and_cropping() -> None:
    pytest.importorskip("pydicom")
    ds = _FakeTiled(rows=4, cols=4, total_rows=10, total_cols=11, channels=3)
    assert ds.NumberOfFrames == 9  # 3 rows x 3 cols of tiles
    lvl = idc.assemble_level(ds)
    assert lvl.shape == (10, 11, 3)
    assert lvl[0, 0, 0] == 0 and lvl[0, 4, 0] == 1 and lvl[0, 8, 0] == 2
    assert lvl[4, 0, 0] == 3 and lvl[9, 10, 2] == 8
    mono = idc.assemble_level(_FakeTiled(rows=4, cols=4, total_rows=5, total_cols=4, channels=0))
    assert mono.shape == (5, 4) and mono[4, 0] == 1
    single = _FakeTiled(rows=6, cols=7, total_rows=6, total_cols=7, channels=3)
    assert idc.assemble_level(single).shape == (6, 7, 3)


# ---------------------------------------------------------------------------
# raw-bucket helpers (offline parts)


def test_viewer_url_and_series_size_from_row() -> None:
    s = idc.series("tcga_lihc", modality="CT", max_size_mb=0.5)
    row = s.sort_values("series_size_MB").iloc[0]
    assert row["crdc_series_uuid"] == "6a3683fe-035c-4932-ac6c-14c0c9005da6"
    assert idc.viewer_url(row) == (
        "https://viewer.imaging.datacommons.cancer.gov/viewer/"
        "1.3.6.1.4.1.14519.5.2.1.3344.4008.164131928981776564121389657757"
        "?seriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335"
    )
    assert idc.series_size(row) == 0.409608
    assert idc.viewer_url(row.to_dict()) == idc.viewer_url(row)
    with pytest.raises(ValueError):
        idc.viewer_url("not a row")
    with pytest.raises(ValueError):
        idc.pull_series(pd.Series({"x": 1}), "/nonexistent")


# ---------------------------------------------------------------------------
# live S3


@pytest.mark.network
def test_live_acrin_6698_series_column_subset(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SCIGANTIC_NCI_ROOT", raising=False)
    monkeypatch.setenv("SCIGANTIC_MOUNT_PATH", "/nonexistent-mount")
    s = idc.series("acrin_6698", modality="SEG", columns=["SeriesInstanceUID", "PatientID", "series_size_MB"])
    assert list(s.columns) == ["SeriesInstanceUID", "PatientID", "series_size_MB"]
    assert len(s) == 2213
    assert s["PatientID"].nunique() == 385


@pytest.mark.network
def test_live_pull_smallest_ct_series(tmp_path: object) -> None:
    pydicom = pytest.importorskip("pydicom")
    s = idc.series("tcga_lihc", modality="CT")
    row = s.sort_values("series_size_MB").iloc[0]
    dest = os.path.join(str(tmp_path), "ct")
    files = idc.pull_series(row, dest)
    assert len(files) == int(row["instanceCount"]) == 1
    assert all(f.endswith(".dcm") and os.path.dirname(f) == dest for f in files)
    assert abs(sum(os.path.getsize(f) for f in files) / 1e6 - float(row["series_size_MB"])) < 0.01
    ds = pydicom.dcmread(files[0])
    assert ds.Modality == "CT" and ds.SeriesInstanceUID == row["SeriesInstanceUID"]
    assert ds.pixel_array.shape == (int(ds.Rows), int(ds.Columns))
    # same series by uuid string, bucket defaulted to idc-open-data
    assert idc.series_size(row["crdc_series_uuid"]) == pytest.approx(float(row["series_size_MB"]), abs=1e-6)
