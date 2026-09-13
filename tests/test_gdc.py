"""Tests for scigantic_nci.gdc against the local GDC mirror (see conftest.gdc_root).

Every number here was measured on the Data Release 46.0 build of the mirror.
"""
from __future__ import annotations

import gzip
import os

import numpy as np
import pandas as pd
import pytest

from scigantic_nci import gdc
from scigantic_nci._store import GDC_BUCKET, NciError, NotMirroredError, Store

BLCA = "TCGA-BLCA"


# ------------------------------------------------------------------ project level


def test_projects_lists_mirror(gdc_root: str) -> None:
    p = gdc.projects()
    assert {"project", "gdc_data_release", "built_at", "n_tables", "cases", "files_read", "bytes_read"} <= set(p.columns)
    assert p["project"].is_unique
    row = p.set_index("project").loc[BLCA]
    assert row["cases"] == 412
    assert row["n_tables"] == 35
    assert row["source_files_in_project"] == 10603
    assert str(row["gdc_data_release"]).startswith("Data Release 46.0")
    # cached per process: a second call is a copy of the same frame
    p2 = gdc.projects()
    assert p2 is not p and p2.equals(p)


def test_report_tables_readme(gdc_root: str) -> None:
    rep = gdc.report(BLCA)
    assert rep["project"] == BLCA
    assert rep["tables"]["clinical/cases"]["rows"] == 412
    t = gdc.tables(BLCA)
    assert list(t.columns) == ["table", "rows", "columns", "bytes"]
    assert len(t) == 35
    assert t.set_index("table").loc["MANIFEST", "rows"] == 10603
    assert gdc.readme(BLCA).startswith("# TCGA-BLCA")


def test_manifest(gdc_root: str) -> None:
    m = gdc.manifest(BLCA)
    assert len(m) == 10603
    assert {"file_id", "md5sum", "s3_bucket", "s3_key", "gdc_download_url", "data_type"} <= set(m.columns)
    maf = gdc.manifest(BLCA, data_type="Masked Somatic Mutation")
    assert len(maf) == 415
    assert (maf["data_type"] == "Masked Somatic Mutation").all()
    with pytest.raises(ValueError, match="Masked Somatic Mutation"):
        gdc.manifest(BLCA, data_type="nope")


def test_unknown_project_lists_valid_ids(gdc_root: str) -> None:
    with pytest.raises(NotMirroredError, match=BLCA):
        gdc.clinical("TCGA-NOPE")


# ---------------------------------------------------------------------- clinical


def test_clinical_one_row_per_case(gdc_root: str) -> None:
    c = gdc.clinical(BLCA)
    assert len(c) == 412 and c["case_id"].is_unique
    assert c.attrs["stage_column"] == "ajcc_pathologic_stage"
    assert c.attrs["project"] == BLCA
    # clashing diagnosis columns get the _dx suffix
    assert {"submitter_id_dx", "state_dx", "updated_datetime_dx"} <= set(c.columns)
    # helper columns
    assert c["os_event"].dtype == bool and int(c["os_event"].sum()) == 182
    assert int((c["demographic_vital_status"] == "Dead").sum()) == 182
    assert int(c["os_time_days"].isna().sum()) == 50
    dead = c[c["os_event"]]
    assert dead["os_time_days"].equals(dead["demographic_days_to_death"])
    alive = c[~c["os_event"]]
    assert alive["os_time_days"].equals(alive["days_to_last_follow_up"])
    years = c["age_at_diagnosis_years"].dropna()
    assert 34 < years.min() < 35 and years.max() == pytest.approx(90.0, abs=0.01)
    assert c["stage"].value_counts().to_dict() == {"Stage III": 141, "Stage IV": 136, "Stage II": 131, "Stage I": 1}


def test_clinical_all_diagnoses(gdc_root: str) -> None:
    c = gdc.clinical(BLCA, primary_only=False)
    assert len(c) == 1014 and c["case_id"].nunique() == 412


@pytest.mark.parametrize(
    "project,stage_column",
    [
        ("TCGA-ACC", "ensat_pathologic_stage"),
        ("TCGA-OV", "figo_stage"),
        ("TCGA-THYM", "masaoka_stage"),
        ("TCGA-LGG", "tumor_grade"),
        ("TCGA-GBM", None),
    ],
)
def test_clinical_stage_column_choice(gdc_root: str, project: str, stage_column: str | None) -> None:
    c = gdc.clinical(project)
    assert c.attrs["stage_column"] == stage_column
    if stage_column is None:
        assert c["stage"].isna().all()
    else:
        assert c["stage"].notna().sum() == c[stage_column].notna().sum() > 0


def test_clinical_thin_wrappers(gdc_root: str) -> None:
    assert len(gdc.diagnoses(BLCA)) == 1014
    assert len(gdc.treatments(BLCA)) == 2306
    assert len(gdc.exposures(BLCA)) == 677
    assert len(gdc.follow_ups(BLCA)) == 1967
    assert len(gdc.samples(BLCA)) == 1240
    assert len(gdc.aliquots(BLCA)) == 4836
    with pytest.raises(NotMirroredError, match="clinical/exposures"):
        gdc.exposures("TCGA-CHOL")


# -------------------------------------------------------------------- expression


def test_expression_samples(gdc_root: str) -> None:
    s = gdc.expression_samples(BLCA)
    assert len(s) == 431
    assert s["sample_type"].value_counts().to_dict() == {"Primary Tumor": 412, "Solid Tissue Normal": 19}
    assert len(gdc.star_qc(BLCA)) == 431


def test_expression_gene_subset(gdc_root: str) -> None:
    e = gdc.expression(BLCA, genes=["TP53", "ESR1", "NOPE"], sample_type="Primary Tumor")
    assert e.shape == (2, 412)
    assert set(e.index) == {"TP53", "ESR1"} and e.index.name == "gene_name"
    assert e.attrs["missing_genes"] == ["NOPE"]
    assert e.attrs["gene_id"]["TP53"] == "ENSG00000141510.18"
    assert e.attrs["gene_type"]["TP53"] == "protein_coding"
    assert e.dtypes.iloc[0] == np.float32
    assert e.attrs["kind"] == "tpm"
    # attrs survive slicing and concat (a plain Series in attrs would break repr/concat)
    assert pd.concat([e.iloc[:, :2], e.iloc[:, 2:4]], axis=1).attrs["gene_id"]["ESR1"] == "ENSG00000091831.24"
    with pytest.raises(ValueError, match="none of"):
        gdc.expression(BLCA, genes=["NOPE"])


def test_expression_full_matrix_and_duplicates(gdc_root: str) -> None:
    e = gdc.expression(BLCA, sample_type="Solid Tissue Normal")
    # 60,660 GENCODE v36 rows minus 44 _PAR_Y ids minus 1,189 further duplicated names
    assert e.shape == (59427, 19)
    assert e.index.is_unique
    dropped = e.attrs["dropped_duplicate_names"]
    assert len(dropped) == 67 and "Y_RNA" in dropped and "SNORD116" in dropped
    pc = gdc.expression(BLCA, sample_type="Solid Tissue Normal", protein_coding=True, kind="counts")
    assert pc.shape == (19938, 19)
    assert pc.dtypes.iloc[0] == np.int32
    with pytest.raises(ValueError, match="kind"):
        gdc.expression(BLCA, kind="fpkm")
    with pytest.raises(ValueError, match="Metastatic"):
        gdc.expression(BLCA, sample_type="Metastatic")


# --------------------------------------------------------------------- mutations


def test_mutations(gdc_root: str) -> None:
    m = gdc.mutations(BLCA)
    assert m.shape == (117053, 17)
    assert list(m.columns) == list(gdc.DEFAULT_MAF_COLUMNS)
    assert m["Tumor_Sample_Barcode"].nunique() == 414
    assert gdc.mutations(BLCA, columns="all").shape == (117053, 142)
    tp53 = gdc.mutations(BLCA, genes=["TP53"])
    assert len(tp53) == 233 and (tp53["Hugo_Symbol"] == "TP53").all()
    ns = gdc.mutations(BLCA, genes="TP53", non_silent=True, columns=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"])
    assert len(ns) == 226 and ns["Tumor_Sample_Barcode"].nunique() == 202
    assert set(ns["Variant_Classification"]) <= gdc.NON_SILENT
    assert list(ns.columns) == ["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"]
    with pytest.raises(ValueError, match="unknown MAF columns"):
        gdc.mutations(BLCA, columns=["Hugo_Symbol", "Nope"])


def test_mutation_frequency(gdc_root: str) -> None:
    f = gdc.mutation_frequency(BLCA)
    assert f.attrs["n_samples"] == 414
    assert f.index[0] == "TP53"
    assert f["TP53"] == pytest.approx(202 / 414)
    assert abs(f["TP53"] - 0.49) < 0.01
    assert f["TTN"] == pytest.approx(0.425, abs=0.001)
    everything = gdc.mutation_frequency(BLCA, non_silent=False)
    assert everything.index[0] == "TTN" and everything["TTN"] == pytest.approx(0.522, abs=0.001)


# ------------------------------------------------------------------- copy number


def test_copy_number_workflows(gdc_root: str) -> None:
    assert gdc.copy_number_workflows(BLCA) == ["absolute_liftover", "ascat2", "ascat3", "ascatngs"]
    assert gdc.copy_number_workflows(BLCA, "segments") == ["dnacopy", "gatk4_cnv"]
    assert gdc.copy_number_workflows(BLCA, "masked") == ["dnacopy"]
    assert gdc.copy_number_workflows(BLCA, "allele_specific") == ["ascat2", "ascat3", "ascatngs"]
    with pytest.raises(ValueError, match="kind"):
        gdc.copy_number_workflows(BLCA, "nope")


def test_copy_number_matrix(gdc_root: str) -> None:
    cn = gdc.copy_number(BLCA, genes=["ERBB2", "CDKN2A", "TP53"])
    assert cn.shape == (3, 392) and cn.attrs["workflow"] == "ascat3"
    assert cn.attrs["gene_location"]["TP53"] == ("chr17", 7661779, 7687538)
    assert cn.attrs["gene_id"]["ERBB2"].startswith("ENSG00000141736")
    full = gdc.copy_number(BLCA, workflow="ascat3")
    assert full.shape == (59390, 392) and full.index.is_unique
    assert gdc.copy_number(BLCA, workflow="ascatngs").shape[1] == 13
    with pytest.raises(ValueError, match="available"):
        gdc.copy_number(BLCA, workflow="nope")
    with pytest.raises(ValueError, match="none of"):
        gdc.copy_number(BLCA, genes=["NOPE"])


def test_copy_number_samples_sample_type_is_not_the_column(gdc_root: str) -> None:
    s = gdc.copy_number_samples(BLCA)
    assert len(s) == 392 and s.attrs["workflow"] == "ascat3"
    # every matrix column is the tumor aliquot (-01A) although sample_type often says normal
    assert s["column"].str[13:15].eq("01").all()
    assert s["sample_type"].value_counts().to_dict() == {"Blood Derived Normal": 195, "Primary Tumor": 186, "Solid Tissue Normal": 11}


def test_segments(gdc_root: str) -> None:
    seg = gdc.segments(BLCA)
    assert seg.shape == (588674, 8) and seg.attrs["workflow"] == "dnacopy"
    assert gdc.segments(BLCA, workflow="gatk4_cnv").shape == (2386633, 8)
    assert gdc.segments(BLCA, "masked").shape == (266389, 8)
    asc = gdc.segments(BLCA, "allele_specific")
    assert asc.shape == (49849, 9) and asc.attrs["workflow"] == "ascat3"
    assert {"Copy_Number", "Major_Copy_Number", "Minor_Copy_Number"} <= set(asc.columns)
    with pytest.raises(ValueError, match="kind"):
        gdc.segments(BLCA, "nope")


# ---------------------------------------------------------------- miRNA and RPPA


def test_mirna(gdc_root: str) -> None:
    rpm = gdc.mirna(BLCA)
    assert rpm.shape == (1881, 437) and rpm.index.name == "miRNA_ID"
    counts = gdc.mirna(BLCA, kind="counts", sample_type="Solid Tissue Normal")
    assert counts.shape == (1881, 19)
    assert len(gdc.mirna_samples(BLCA)) == 437
    with pytest.raises(ValueError, match="kind"):
        gdc.mirna(BLCA, kind="tpm")


def test_rppa(gdc_root: str) -> None:
    r = gdc.rppa(BLCA)
    assert r.shape == (487, 343) and r.index.name == "peptide_target"
    assert r.attrs["agid"]["1433BETA"] == "AGID00100"
    assert len(gdc.rppa_antibodies(BLCA)) == 487
    assert len(gdc.rppa_samples(BLCA)) == 343


# ------------------------------------------------------------------- methylation


def test_methylation_absent(gdc_root: str) -> None:
    assert gdc.methylation_platforms(BLCA) == []
    with pytest.raises(NotMirroredError, match="methylation"):
        gdc.methylation(BLCA)


def test_methylation_probe_and_column_subsets(gdc_methylation_root: str) -> None:
    root = gdc_methylation_root
    assert gdc.methylation_platforms("TCGA-CHOL", root=root) == ["human_methylation_450"]
    s = gdc.methylation_samples("TCGA-CHOL", root=root)
    assert len(s) == 45 and s["sample_type"].value_counts().to_dict() == {"Primary Tumor": 36, "Solid Tissue Normal": 9}
    two = gdc.methylation("TCGA-CHOL", probes=["cg00000029", "cg00000165", "nope"], root=root)
    assert two.shape == (2, 45) and two.index.name == "probe_id" and two.attrs["platform"] == "human_methylation_450"
    assert two.dtypes.iloc[0] == np.float32
    assert 0.0 <= float(two.loc["cg00000029"].max()) <= 1.0
    normals = gdc.methylation("TCGA-CHOL", sample_type="Solid Tissue Normal", root=root)
    assert normals.shape == (486427, 9)
    one = gdc.methylation("TCGA-CHOL", columns=[s["column"].iloc[0]], probes="cg00000029", root=root)
    assert one.shape == (1, 1)
    with pytest.raises(ValueError, match="unknown methylation columns"):
        gdc.methylation("TCGA-CHOL", columns=["nope"], root=root)
    with pytest.raises(ValueError, match="platform"):
        gdc.methylation("TCGA-CHOL", platform="epic", root=root)
    with pytest.raises(ValueError, match="probe_ids"):
        gdc.methylation("TCGA-CHOL", probes=["nope"], root=root)


# ----------------------------------------------------------------------- helpers


def test_case_for_barcode() -> None:
    assert gdc.case_for_barcode("TCGA-DK-A6B2-01A-11R-A30C-07") == "TCGA-DK-A6B2"
    assert gdc.case_for_barcode("TCGA-DK-A6B2") == "TCGA-DK-A6B2"
    assert gdc.case_for_barcode("TARGET-30-PAAPFA-01A-01R") == "TARGET-30-PAAPFA"
    assert gdc.case_for_barcode("C3L-00001-01") == "C3L-00001-01"


def test_mutation_barcodes_map_to_cases(gdc_root: str) -> None:
    m = gdc.mutations(BLCA, genes=["TP53"], columns=["Tumor_Sample_Barcode"])
    cases = set(gdc.clinical(BLCA)["submitter_id"])
    assert set(m["Tumor_Sample_Barcode"].map(gdc.case_for_barcode)) <= cases


# ----------------------------------------------------------------------- network


@pytest.mark.network
def test_fetch_raw_smallest_maf(gdc_root: str, tmp_path: os.PathLike[str]) -> None:
    row = gdc.manifest(BLCA, data_type="Masked Somatic Mutation").sort_values("file_size").iloc[0]
    assert row["file_size"] == 1143
    data = gdc.fetch_raw(row)
    assert isinstance(data, bytes) and len(data) == 1143
    assert gzip.decompress(data).startswith(b"#version gdc-1.0.0")
    path = gdc.fetch_raw(BLCA, row["file_id"], dest=str(tmp_path))
    assert isinstance(path, str) and os.path.basename(path) == row["file_name"] and os.path.getsize(path) == 1143
    bad = row.copy()
    bad["md5sum"] = "0" * 32
    with pytest.raises(NciError, match="md5 mismatch"):
        gdc.fetch_raw(bad)
    with pytest.raises(ValueError, match="not in the TCGA-BLCA manifest"):
        gdc.fetch_raw(BLCA, "nope")


@pytest.mark.network
def test_live_s3_chol(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SCIGANTIC_NCI_ROOT", raising=False)
    if not Store(GDC_BUCKET, "TCGA-CHOL", None).exists("BUILD_REPORT.json"):
        pytest.skip("s3://scigantic-gdc-open/TCGA-CHOL/ is not uploaded yet")
    assert gdc.report("TCGA-CHOL")["project"] == "TCGA-CHOL"
    assert len(gdc.clinical("TCGA-CHOL")) == 51
    s = gdc.expression_samples("TCGA-CHOL")
    assert s["sample_type"].value_counts().to_dict() == {"Primary Tumor": 35, "Solid Tissue Normal": 9}
