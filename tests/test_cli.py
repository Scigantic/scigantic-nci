"""Tests for the scigantic-nci console script, run through ``main()`` against the local mirrors.

Numbers are the ones measured in test_gdc.py and test_idc.py (TCGA-BLCA, tcga_lihc).
"""
from __future__ import annotations

import io
import os

import pandas as pd
import pytest

from scigantic_nci.cli import main

BLCA = "TCGA-BLCA"
LIHC = "tcga_lihc"
CT_UID = "1.3.6.1.4.1.14519.5.2.1.3344.4008.615841149820206875475663260335"


def _csv(capsys: pytest.CaptureFixture[str]) -> pd.DataFrame:
    out = capsys.readouterr().out
    return pd.read_csv(io.StringIO(out))


# ------------------------------------------------------------------------ GDC


def test_projects_csv(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["projects"]) == 0
    df = _csv(capsys)
    assert {"project", "cases", "n_tables", "gdc_data_release"} <= set(df.columns)
    assert int(df.set_index("project").loc[BLCA, "cases"]) == 412


def test_tables_gdc(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["tables", BLCA]) == 0
    df = _csv(capsys)
    assert list(df.columns) == ["table", "rows", "columns", "bytes"]
    assert len(df) == 35
    assert int(df.set_index("table").loc["MANIFEST", "rows"]) == 10603


def test_tables_with_explicit_root(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    root = os.path.join(gdc_root, BLCA)
    assert main(["--root", root, "tables", BLCA]) == 0
    assert len(_csv(capsys)) == 35


def test_readme_gdc(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["readme", BLCA]) == 0
    out = capsys.readouterr().out
    assert out.startswith("# TCGA-BLCA") and out.endswith("\n")


def test_clinical_to_file(gdc_root: str, tmp_path: os.PathLike[str], capsys: pytest.CaptureFixture[str]) -> None:
    out = os.path.join(str(tmp_path), "clin.csv")
    assert main(["clinical", BLCA, "--out", out]) == 0
    assert capsys.readouterr().out == ""
    df = pd.read_csv(out)
    assert len(df) == 412 and df["case_id"].is_unique
    assert {"os_time_days", "os_event", "age_at_diagnosis_years", "stage"} <= set(df.columns)
    assert int(df["os_event"].sum()) == 182


def test_expression_genes_and_sample_type(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["expression", BLCA, "--genes", "TP53,ESR1", "--sample-type", "Solid Tissue Normal"]) == 0
    df = _csv(capsys).set_index("gene_name")
    assert sorted(df.index) == ["ESR1", "TP53"]
    assert df.shape == (2, 19)
    assert df.columns.str.startswith("TCGA-").all()


def test_expression_counts_are_integers(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["expression", BLCA, "--genes", "TP53", "--counts", "--sample-type", "Solid Tissue Normal"]) == 0
    df = _csv(capsys).set_index("gene_name")
    assert df.shape == (1, 19)
    assert (df.dtypes == "int64").all()


def test_mutations_non_silent(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["mutations", BLCA, "--genes", "TP53,RB1", "--non-silent"]) == 0
    df = _csv(capsys)
    assert set(df["Hugo_Symbol"]) == {"TP53", "RB1"}
    assert "Silent" not in set(df["Variant_Classification"])
    assert {"Tumor_Sample_Barcode", "HGVSp_Short", "t_depth"} <= set(df.columns)
    assert len(df.columns) == 17


def test_mutation_frequency_top(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["mutation-frequency", BLCA, "--top", "3"]) == 0
    df = _csv(capsys)
    assert list(df.columns) == ["Hugo_Symbol", "fraction_mutated"]
    assert df["Hugo_Symbol"].tolist() == ["TP53", "TTN", "KMT2D"]
    assert round(float(df["fraction_mutated"].iloc[0]), 3) == 0.488


def test_unknown_project_exit_code_2(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["clinical", "TCGA-NOPE"]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "TCGA-NOPE" in captured.err and BLCA in captured.err


def test_bad_sample_type_exit_code_2(gdc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["expression", BLCA, "--genes", "TP53", "--sample-type", "Metastatic"]) == 2
    assert "Metastatic" in capsys.readouterr().err


# ------------------------------------------------------------------------ IDC


def test_collections_csv(idc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["collections"]) == 0
    df = _csv(capsys)
    assert len(df) == 176
    assert int(df.set_index("collection_id").loc[LIHC, "n_series"]) == 3457


def test_tables_idc(idc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["tables", LIHC]) == 0
    df = _csv(capsys)
    assert list(df.columns) == ["table", "rows", "columns", "bytes"]
    assert int(df.set_index("table").loc["series", "rows"]) == 3457
    assert "clinical/dictionary" in set(df["table"])


def test_readme_idc(idc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["readme", LIHC]) == 0
    assert capsys.readouterr().out.startswith("# TCGA-LIHC (tcga_lihc)")


def test_series_modality(idc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["series", LIHC, "--modality", "CT"]) == 0
    df = _csv(capsys)
    assert len(df) == 777 and (df["Modality"] == "CT").all()
    assert main(["series", LIHC, "--modality", "XX"]) == 2
    assert "modalities: MR, SEG, SM, CT, PT" in capsys.readouterr().err


def test_sample_describes_ct(idc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["sample", LIHC, "--modality", "CT"]) == 0
    out = capsys.readouterr().out
    assert "modality: CT\n" in out and "kind: series\n" in out
    assert "PatientID: TCGA-DD-A3A6\n" in out
    assert "dcm_files: 36\n" in out
    assert out.count(".dcm\n") == 36
    assert main(["sample", LIHC, "--modality", "SM"]) == 2
    assert "sample modalities: CT, MR, PT" in capsys.readouterr().err


def test_viewer_url_by_series_uid_and_uuid(idc_root: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["viewer-url", LIHC, "--series-uid", CT_UID]) == 0
    url = capsys.readouterr().out.strip()
    assert url == (
        "https://viewer.imaging.datacommons.cancer.gov/viewer/"
        "1.3.6.1.4.1.14519.5.2.1.3344.4008.164131928981776564121389657757"
        f"?seriesInstanceUID={CT_UID}"
    )
    assert main(["viewer-url", LIHC, "--uuid", "6a3683fe-035c-4932-ac6c-14c0c9005da6"]) == 0
    assert capsys.readouterr().out.strip() == url
    assert main(["viewer-url", LIHC, "--series-uid", "1.2.3"]) == 2
    assert "no series with SeriesInstanceUID '1.2.3'" in capsys.readouterr().err


@pytest.mark.network
def test_pull_series_smallest_ct(idc_root: str, tmp_path: os.PathLike[str], capsys: pytest.CaptureFixture[str]) -> None:
    dest = os.path.join(str(tmp_path), "ct")
    assert main(["pull-series", LIHC, "--series-uid", CT_UID, "--dest", dest]) == 0
    paths = capsys.readouterr().out.splitlines()
    assert len(paths) == 1 and paths[0].endswith(".dcm") and os.path.dirname(paths[0]) == dest
    assert abs(os.path.getsize(paths[0]) / 1e6 - 0.409608) < 0.01
