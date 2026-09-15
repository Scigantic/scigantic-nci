"""GDC projects: clinical, RNA-Seq, mutations, copy number, miRNA, RPPA, methylation.

Every function takes the GDC project id first (``"TCGA-BRCA"``) and an optional
``root`` last. ``root`` is a local directory holding that project's tables; when it
is omitted the resolver in ``_store`` picks ``$SCIGANTIC_NCI_ROOT``, the notebook
mount, or the public bucket ``s3://scigantic-gdc-open/<PROJECT>/``.

Wide matrices (expression, copy number, miRNA, RPPA, methylation) come back as
features x aliquots DataFrames. Column ids are aliquot barcodes; the sibling
``*_samples`` functions map each column to ``case_submitter_id`` and ``sample_type``.
Functions that take ``genes`` accept gene symbols and match on ``gene_name``.

Nothing here prints, and column subsets are pushed down to parquet so asking for
five genes or two probes never materialises a whole matrix in pandas.
"""
from __future__ import annotations

import hashlib
import json
import os
from functools import lru_cache
from typing import Any, Iterable, Sequence

import pandas as pd

from ._store import GDC_BUCKET, NciError, NotMirroredError, Store, catalog_rows, list_keys, resolve

__all__ = [
    "NON_SILENT",
    "DEFAULT_MAF_COLUMNS",
    "projects",
    "report",
    "tables",
    "readme",
    "manifest",
    "clinical",
    "diagnoses",
    "treatments",
    "exposures",
    "follow_ups",
    "samples",
    "aliquots",
    "expression",
    "expression_samples",
    "star_qc",
    "mutations",
    "mutation_frequency",
    "copy_number",
    "copy_number_samples",
    "copy_number_workflows",
    "segments",
    "mirna",
    "mirna_samples",
    "mirna_isoforms",
    "rppa",
    "rppa_samples",
    "rppa_antibodies",
    "methylation",
    "methylation_samples",
    "methylation_platforms",
    "case_for_barcode",
    "fetch_raw",
]

#: MAF Variant_Classification values that change the protein (the maftools
#: "nonSyn" set). Everything else (Silent, Intron, UTRs, Flanks, IGR, RNA,
#: Splice_Region) counts as silent for ``non_silent=True``.
NON_SILENT: frozenset[str] = frozenset(
    {
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Nonstop_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Splice_Site",
        "Translation_Start_Site",
    }
)

#: Columns ``mutations()`` reads when ``columns`` is None (17 of the 142).
DEFAULT_MAF_COLUMNS: tuple[str, ...] = (
    "file_id",
    "aliquot_submitter_id",
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "HGVSc",
    "HGVSp_Short",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
)

_EXPRESSION_KINDS = {"tpm": "rna_seq/tpm_unstranded.parquet", "counts": "rna_seq/counts_unstranded.parquet"}
_MIRNA_KINDS = {"rpm": "mirna/rpm.parquet", "counts": "mirna/read_count.parquet"}
_SEGMENT_KINDS = {"segments": "segments", "masked": "masked_segments", "allele_specific": "allele_specific_segments"}
_WORKFLOW_PREFERENCE = ("ascat3", "ascat2", "absolute_liftover", "ascatngs", "dnacopy", "gatk4_cnv")


class _Frozen(dict[str, Any]):
    """A dict for ``DataFrame.attrs`` that pandas can carry through every operation for free.

    pandas deep-copies ``attrs`` on most operations and compares them on concat; a
    plain 60,000-entry dict would add about 15 ms to every slice. Deep copy and copy
    of this class return the same object, so treat it as read-only.
    """

    def __deepcopy__(self, memo: dict[int, Any]) -> _Frozen:
        return self

    def __copy__(self) -> _Frozen:
        return self


# --------------------------------------------------------------------------- ids


def _local_roots() -> list[str]:
    """Directories that may hold ``<PROJECT>/BUILD_REPORT.json`` locally."""
    out: list[str] = []
    env_root = os.environ.get("SCIGANTIC_NCI_ROOT")
    if env_root:
        out.extend(p for p in (os.path.join(env_root, GDC_BUCKET), env_root) if os.path.isdir(p))
    return out


@lru_cache(maxsize=8)
def _project_ids(env_root: str | None) -> tuple[str, ...]:
    ids: set[str] = set()
    for base in _local_roots():
        for name in sorted(os.listdir(base)):
            rep = os.path.join(base, name, "BUILD_REPORT.json")
            if os.path.isfile(rep):
                try:
                    with open(rep) as fh:
                        d = json.load(fh)
                except (OSError, ValueError):
                    continue
                if isinstance(d, dict) and d.get("project") == name:
                    ids.add(name)
    mount = os.environ.get("SCIGANTIC_MOUNT_PATH", "/mnt/archive")
    rep = os.path.join(mount, "BUILD_REPORT.json")
    if os.path.isfile(rep):
        try:
            with open(rep) as fh:
                d = json.load(fh)
            if isinstance(d, dict) and isinstance(d.get("project"), str):
                ids.add(d["project"])
        except (OSError, ValueError):
            pass
    try:
        ids.update(list_keys(GDC_BUCKET))
    except Exception:  # offline; local ids are still valid
        if not ids:
            raise
    return tuple(sorted(ids))


def _known_projects() -> tuple[str, ...]:
    return _project_ids(os.environ.get("SCIGANTIC_NCI_ROOT"))


def _store(project: str, root: str | None) -> Store:
    if root is None:
        known = _known_projects()
        if project not in known:
            raise NotMirroredError(
                f"{project!r} is not a mirrored GDC project. "
                f"Mirrored ({len(known)}): {', '.join(known) if known else 'none yet'}"
            )
    return resolve(GDC_BUCKET, project, root)


def _read(store: Store, rel: str, columns: list[str] | None = None, filters: Any = None) -> pd.DataFrame:
    """``store.read_parquet`` with a NotMirroredError naming what is there when ``rel`` is absent."""
    try:
        return store.read_parquet(rel, columns=columns, filters=filters)
    except (FileNotFoundError, OSError) as exc:
        if store.exists(rel):
            raise
        group = os.path.dirname(rel)
        have = [f for f in store.listdir(group) if f.endswith(".parquet")] if group else store.listdir()
        raise NotMirroredError(
            f"{store.key} has no {rel} (in {group or 'root'}: {', '.join(have) if have else 'nothing'})"
        ) from exc


def _as_list(x: str | Iterable[str] | None) -> list[str] | None:
    if x is None:
        return None
    return [x] if isinstance(x, str) else list(x)


# ------------------------------------------------------------------ project level


_PROJECT_COLUMNS: tuple[str, ...] = (
    "project",
    "gdc_data_release",
    "built_at",
    "n_tables",
    "cases",
    "source_files_in_project",
    "files_read",
    "bytes_read",
)


def _project_row(pid: str) -> dict[str, Any]:
    store = resolve(GDC_BUCKET, pid)
    rep = store.read_json("BUILD_REPORT.json")
    n_cases: int | None = None
    if store.exists("clinical/cases.parquet"):
        n_cases = len(store.read_parquet("clinical/cases.parquet", columns=["case_id"]))
    return {
        "project": rep.get("project", pid),
        "gdc_data_release": rep.get("gdc_data_release"),
        "built_at": rep.get("built_at"),
        "n_tables": len(rep.get("tables", {})),
        "cases": n_cases,
        "source_files_in_project": rep.get("source_files_in_project"),
        "files_read": rep.get("files_read"),
        "bytes_read": rep.get("bytes_read"),
    }


def _build_projects(use_index: bool = True) -> pd.DataFrame:
    """The projects() frame, uncached. ``use_index=False`` ignores PROJECTS.parquet (to write it)."""
    rows = catalog_rows(GDC_BUCKET, _known_projects(), _project_row, _PROJECT_COLUMNS, "project", use_index=use_index)
    return pd.DataFrame(rows, columns=list(_PROJECT_COLUMNS))


@lru_cache(maxsize=8)
def _projects_frame(env_root: str | None) -> pd.DataFrame:
    return _build_projects()


def projects() -> pd.DataFrame:
    """One row per mirrored project, from each BUILD_REPORT.json plus the case count.

    Cached for the life of the process (per ``$SCIGANTIC_NCI_ROOT`` value).
    Columns: project, gdc_data_release, built_at, n_tables, cases,
    source_files_in_project, files_read, bytes_read.

    The ids come from listing the bucket (and any local copy). Rows are read from the
    bucket's ``PROJECTS.parquet`` index when it has them and built from the project's own
    files otherwise (16 at a time). A project rebuilt after the index was written keeps its
    previous row (``built_at`` and counts) until the index is regenerated.
    """
    return _projects_frame(os.environ.get("SCIGANTIC_NCI_ROOT")).copy()


def report(project: str, root: str | None = None) -> dict[str, Any]:
    """The project's BUILD_REPORT.json as a dict (release, build time, per-table rows/columns/bytes)."""
    d = _store(project, root).read_json("BUILD_REPORT.json")
    return dict(d)


def tables(project: str, root: str | None = None) -> pd.DataFrame:
    """Every table of the project with rows, columns and bytes (from the build report)."""
    t = report(project, root).get("tables", {})
    df = pd.DataFrame(
        [{"table": k, "rows": v.get("rows"), "columns": v.get("columns"), "bytes": v.get("bytes")} for k, v in t.items()],
        columns=["table", "rows", "columns", "bytes"],
    )
    return df.sort_values("table").reset_index(drop=True)


def readme(project: str, root: str | None = None) -> str:
    """The project's README.md text (what each table is and how it was built)."""
    return _store(project, root).read_bytes("README.md").decode("utf-8")


def manifest(project: str, data_type: str | None = None, root: str | None = None) -> pd.DataFrame:
    """Every open-access file of the project (MANIFEST.parquet), optionally one ``data_type``.

    Columns include file_id, file_name, data_type, file_size, md5sum, case_submitter_id,
    sample_type, aliquot_submitter_id, s3_bucket, s3_key and gdc_download_url.
    Rows of a tumor/normal pair carry both ids pipe-joined in aliquot_submitter_id.
    """
    store = _store(project, root)
    if data_type is None:
        return _read(store, "MANIFEST.parquet")
    df = _read(store, "MANIFEST.parquet", filters=[("data_type", "==", data_type)])
    if df.empty:
        valid = sorted(_read(store, "MANIFEST.parquet", columns=["data_type"])["data_type"].dropna().unique())
        raise ValueError(f"no files of data_type {data_type!r} in {project}; valid: {', '.join(valid)}")
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------- clinical


def _clinical_table(project: str, name: str, root: str | None) -> pd.DataFrame:
    return _read(_store(project, root), f"clinical/{name}.parquet")


def diagnoses(project: str, root: str | None = None) -> pd.DataFrame:
    """clinical/diagnoses: more rows than cases; ``diagnosis_is_primary_disease`` flags the primary."""
    return _clinical_table(project, "diagnoses", root)


def treatments(project: str, root: str | None = None) -> pd.DataFrame:
    """clinical/treatments, one row per treatment record (keyed by case_id and diagnosis_id)."""
    return _clinical_table(project, "treatments", root)


def exposures(project: str, root: str | None = None) -> pd.DataFrame:
    """clinical/exposures (tobacco, alcohol, ...). Absent in some projects (raises NotMirroredError)."""
    return _clinical_table(project, "exposures", root)


def follow_ups(project: str, root: str | None = None) -> pd.DataFrame:
    """clinical/follow_ups, one row per follow-up visit."""
    return _clinical_table(project, "follow_ups", root)


def samples(project: str, root: str | None = None) -> pd.DataFrame:
    """clinical/samples: one row per biospecimen sample (sample_type, tissue_type, preservation_method...)."""
    return _clinical_table(project, "samples", root)


def aliquots(project: str, root: str | None = None) -> pd.DataFrame:
    """clinical/aliquots: aliquot_submitter_id -> sample_submitter_id, case_submitter_id, sample_type.

    This is the join table for non-TCGA/TARGET barcodes (see ``case_for_barcode``).
    """
    return _clinical_table(project, "aliquots", root)


def _is_true(s: pd.Series) -> pd.Series:
    return s.map(lambda v: v is True or str(v).lower() == "true").astype(bool)


def clinical(project: str, primary_only: bool = True, root: str | None = None) -> pd.DataFrame:
    """Cases joined to their diagnosis, one row per case, with survival helper columns.

    With ``primary_only`` (default) each case is joined to the single diagnosis row
    flagged ``diagnosis_is_primary_disease``; a case without such a row keeps NaN
    diagnosis columns but a non-null ``case_submitter_id`` (the case's submitter_id).
    With ``primary_only=False`` every diagnosis row is kept, so a case with
    recurrence/metastasis rows appears more than once.

    Column names that exist in both tables (state, submitter_id, updated_datetime)
    get the suffix ``_dx`` on the diagnosis side. Helper columns added at the end:

    - ``os_time_days``: demographic_days_to_death when vital status is Dead,
      else days_to_last_follow_up from the joined diagnosis row.
    - ``os_event``: True when vital status is Dead.
    - ``age_at_diagnosis_years``: age_at_diagnosis (days) / 365.25.
    - ``stage``: a copy of the best-covered column whose name ends in ``_stage``
      (ajcc_pathologic_stage, figo_stage, ensat_pathologic_stage, ...), falling
      back to tumor_grade when no stage column has values. The source column
      name is in ``df.attrs["stage_column"]`` (None if nothing qualified).
    """
    store = _store(project, root)
    cases = _read(store, "clinical/cases.parquet")
    dx = _read(store, "clinical/diagnoses.parquet")
    if primary_only:
        if "diagnosis_is_primary_disease" in dx.columns:
            dx = dx[_is_true(dx["diagnosis_is_primary_disease"])]
        dx = dx.drop_duplicates("case_id", keep="first")
    df = cases.merge(dx, on="case_id", how="left", suffixes=("", "_dx"))
    # case_submitter_id comes from the diagnoses side of the left join, so a case without a
    # (primary) diagnosis row would carry NaN there; the case's own submitter_id is the same
    # barcode, so take it from the cases table for every row.
    if "submitter_id" in cases.columns:
        df["case_submitter_id"] = df["submitter_id"].to_numpy()

    vital = df.get("demographic_vital_status", pd.Series(index=df.index, dtype=object))
    dead = vital.astype(str).str.lower().eq("dead")
    dtd = pd.to_numeric(df.get("demographic_days_to_death", pd.Series(index=df.index, dtype=float)), errors="coerce")
    dlf = pd.to_numeric(df.get("days_to_last_follow_up", pd.Series(index=df.index, dtype=float)), errors="coerce")
    df["os_time_days"] = dtd.where(dead, dlf)
    df["os_event"] = dead.astype(bool)
    age_days = pd.to_numeric(df.get("age_at_diagnosis", pd.Series(index=df.index, dtype=float)), errors="coerce")
    df["age_at_diagnosis_years"] = age_days / 365.25

    stage_cols = [c for c in df.columns if c.endswith("_stage") and c != "ajcc_staging_system_edition"]
    coverage = {c: int(df[c].notna().sum()) for c in stage_cols}
    stage_column: str | None = None
    if coverage and max(coverage.values()) > 0:
        stage_column = max(stage_cols, key=lambda c: coverage[c])
    elif "tumor_grade" in df.columns and df["tumor_grade"].notna().any():
        stage_column = "tumor_grade"
    df["stage"] = df[stage_column] if stage_column else pd.Series(index=df.index, dtype=object)
    df.attrs["stage_column"] = stage_column
    df.attrs["stage_coverage"] = coverage
    df.attrs["project"] = project
    return df


# ------------------------------------------------------------------- expression


def _genes_table(store: Store) -> pd.DataFrame:
    return _read(store, "rna_seq/genes.parquet")


def _select_columns(smp: pd.DataFrame, sample_type: str | Iterable[str] | None, what: str) -> list[str]:
    """Matrix columns (aliquot barcodes) for the requested sample_type(s)."""
    wanted = _as_list(sample_type)
    if wanted is None:
        return list(smp["column"])
    hit = smp[smp["sample_type"].isin(wanted)]
    if hit.empty:
        counts = smp["sample_type"].value_counts().to_dict()
        raise ValueError(f"no {what} aliquots with sample_type {wanted}; present: {counts}")
    return list(hit["column"])


def _gene_matrix(
    df: pd.DataFrame,
    genes: list[str] | None,
    drop_par_y: bool,
    protein_coding: bool,
) -> pd.DataFrame:
    """Turn a gene_id/gene_name/gene_type + aliquot-column frame into a gene_name-indexed matrix."""
    if drop_par_y:
        df = df[~df["gene_id"].astype(str).str.endswith("_PAR_Y")]
    if protein_coding:
        if "gene_type" not in df.columns:
            raise ValueError("this table has no gene_type column; protein_coding filtering needs rna_seq/genes")
        df = df[df["gene_type"] == "protein_coding"]
    missing: list[str] = []
    if genes is not None:
        present = set(df["gene_name"])
        missing = [g for g in genes if g not in present]
        if len(missing) == len(genes):
            raise ValueError(
                f"none of {genes} match a gene_name (symbols are GENCODE v36 names, e.g. TP53, ESR1, ERBB2)"
            )
        df = df[df["gene_name"].isin(genes)]
    dup_mask = df["gene_name"].duplicated(keep="first")
    dropped = sorted(df.loc[dup_mask, "gene_name"].unique().tolist())
    df = df[~dup_mask]
    meta_cols = [c for c in ("gene_id", "gene_name", "gene_type", "chromosome", "start", "end") if c in df.columns]
    out = df.drop(columns=meta_cols).set_index(pd.Index(df["gene_name"], name="gene_name"))
    names = df["gene_name"].tolist()
    out.attrs["gene_id"] = _Frozen(zip(names, df["gene_id"].tolist()))
    if "gene_type" in df.columns:
        out.attrs["gene_type"] = _Frozen(zip(names, df["gene_type"].tolist()))
    if {"chromosome", "start", "end"} <= set(df.columns):
        out.attrs["gene_location"] = _Frozen(zip(names, zip(df["chromosome"].tolist(), df["start"].tolist(), df["end"].tolist())))
    out.attrs["dropped_duplicate_names"] = dropped
    out.attrs["missing_genes"] = missing
    return out


def expression_samples(project: str, root: str | None = None) -> pd.DataFrame:
    """rna_seq/samples: matrix ``column`` -> aliquot_submitter_id, file_id, case_submitter_id, sample_type, ..."""
    return _read(_store(project, root), "rna_seq/samples.parquet")


def star_qc(project: str, root: str | None = None) -> pd.DataFrame:
    """rna_seq/star_qc: the four STAR summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous) per column."""
    return _read(_store(project, root), "rna_seq/star_qc.parquet")


def expression(
    project: str,
    genes: str | Iterable[str] | None = None,
    sample_type: str | Iterable[str] | None = None,
    kind: str = "tpm",
    protein_coding: bool = False,
    drop_par_y: bool = True,
    root: str | None = None,
) -> pd.DataFrame:
    """STAR gene expression as a genes x aliquots matrix indexed by gene_name.

    ``kind`` is ``"tpm"`` (tpm_unstranded, float32) or ``"counts"`` (counts_unstranded,
    raw int32). ``sample_type`` keeps only aliquots of that type (``"Primary Tumor"``,
    ``"Solid Tissue Normal"``, ``"Metastatic"``; a list is allowed); the columns are
    chosen from rna_seq/samples first so only those aliquot columns are read.
    ``genes`` is a list of symbols; they are matched on gene_name after the read,
    and symbols that do not exist end up in ``df.attrs["missing_genes"]`` (all missing
    raises ValueError). ``protein_coding`` keeps gene_type == protein_coding.

    ``drop_par_y`` (default True) removes the 44 ``_PAR_Y`` pseudoautosomal duplicate
    ids. Duplicated gene names that remain (67 names, mostly snoRNA/snRNA/Y_RNA
    families that GENCODE names identically) keep their first occurrence; the
    dropped names are listed in ``df.attrs["dropped_duplicate_names"]``. The Ensembl
    id and gene_type of each row are in ``df.attrs["gene_id"]`` and
    ``df.attrs["gene_type"]`` (read-only dicts keyed by gene_name).

    A full TCGA-BRCA TPM matrix (60,616 x 1,231, float32) is about 300 MB in pandas.
    """
    if kind not in _EXPRESSION_KINDS:
        raise ValueError(f"kind must be one of {sorted(_EXPRESSION_KINDS)}, got {kind!r}")
    store = _store(project, root)
    smp = expression_samples(project, root)
    cols = _select_columns(smp, sample_type, "RNA-Seq")
    gene_list = _as_list(genes)
    genes_df = _genes_table(store)
    filters: Any = None
    if gene_list is not None:
        ids = genes_df.loc[genes_df["gene_name"].isin(gene_list), "gene_id"].tolist()
        filters = [("gene_id", "in", ids)] if ids else None
    mat = _read(store, _EXPRESSION_KINDS[kind], columns=["gene_id", *cols], filters=filters)
    df = genes_df.merge(mat, on="gene_id", how="inner")
    out = _gene_matrix(df, gene_list, drop_par_y, protein_coding)
    out.attrs.update({"project": project, "kind": kind, "sample_type": _as_list(sample_type)})
    return out


# -------------------------------------------------------------------- mutations


_MAF = "somatic_mutations/masked_somatic_mutations.parquet"


def mutations(
    project: str,
    genes: str | Iterable[str] | None = None,
    non_silent: bool = False,
    columns: Sequence[str] | str | None = None,
    root: str | None = None,
) -> pd.DataFrame:
    """Masked somatic mutations (all aliquot MAFs stacked), one row per variant call.

    ``genes`` restricts to those Hugo_Symbol values with a parquet row filter
    (pushed down, so a handful of genes is read in milliseconds); a symbol with no
    calls simply contributes no rows. ``non_silent`` keeps the ``NON_SILENT``
    Variant_Classification set. ``columns`` is None for ``DEFAULT_MAF_COLUMNS``
    (17 columns), ``"all"`` for every column (142 in TCGA), or an explicit list.

    ``Tumor_Sample_Barcode`` is the tumor aliquot; ``aliquot_submitter_id`` holds the
    normal and tumor aliquots pipe-joined as in the manifest. Use ``case_for_barcode``
    (TCGA/TARGET) or ``aliquots()`` to reach the case.
    """
    store = _store(project, root)
    schema = None
    cols: list[str] | None
    if columns is None:
        schema = store.parquet_schema(_MAF)
        cols = [c for c in DEFAULT_MAF_COLUMNS if c in schema]
    elif isinstance(columns, str):
        if columns != "all":
            raise ValueError("columns must be None, 'all' or a list of MAF column names")
        cols = None
    else:
        schema = store.parquet_schema(_MAF)
        bad = [c for c in columns if c not in schema]
        if bad:
            raise ValueError(f"unknown MAF columns {bad}; available: {', '.join(schema)}")
        cols = list(columns)
    filters: list[tuple[str, str, list[str]]] = []
    gene_list = _as_list(genes)
    if gene_list is not None:
        filters.append(("Hugo_Symbol", "in", gene_list))
    if non_silent:
        filters.append(("Variant_Classification", "in", sorted(NON_SILENT)))
    df = _read(store, _MAF, columns=cols, filters=filters or None).reset_index(drop=True)
    df.attrs.update({"project": project, "non_silent": non_silent, "genes": gene_list})
    return df


def mutation_frequency(project: str, non_silent: bool = True, root: str | None = None) -> pd.Series:
    """Fraction of sequenced tumor samples (unique Tumor_Sample_Barcode in the MAF) carrying
    at least one call per gene, sorted descending. The denominator counts every barcode
    in the MAF, before the ``non_silent`` filter; it is in ``s.attrs["n_samples"]``.
    """
    store = _store(project, root)
    df = _read(store, _MAF, columns=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"])
    n = int(df["Tumor_Sample_Barcode"].nunique())
    if non_silent:
        df = df[df["Variant_Classification"].isin(NON_SILENT)]
    freq = df.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique().sort_values(ascending=False) / max(n, 1)
    freq.name = "fraction_mutated"
    freq.attrs.update({"project": project, "n_samples": n, "non_silent": non_silent})
    return freq


# ------------------------------------------------------------------ copy number


def _cn_workflows(store: Store, kind: str) -> list[str]:
    if kind == "gene_level":
        prefix = "gene_level_"
        names = [
            f[len(prefix) : -len(".parquet")]
            for f in store.listdir("copy_number")
            if f.startswith(prefix) and f.endswith(".parquet") and not f.endswith("_samples.parquet")
        ]
    else:
        prefix = _SEGMENT_KINDS[kind] + "_"
        names = [f[len(prefix) : -len(".parquet")] for f in store.listdir("copy_number") if f.startswith(prefix) and f.endswith(".parquet")]
    return sorted(names)


def _pick_workflow(available: list[str], workflow: str | None, what: str, project: str) -> str:
    if not available:
        raise NotMirroredError(f"{project} has no {what} tables in copy_number/")
    if workflow is None:
        for w in _WORKFLOW_PREFERENCE:
            if w in available:
                return w
        return available[0]
    if workflow not in available:
        raise ValueError(f"unknown {what} workflow {workflow!r} for {project}; available: {available}")
    return workflow


def copy_number_workflows(project: str, kind: str = "gene_level", root: str | None = None) -> list[str]:
    """Workflows with a table of this ``kind``: ``"gene_level"`` (ascat2, ascat3, absolute_liftover,
    ascatngs), ``"segments"`` (dnacopy, gatk4_cnv), ``"masked"`` (dnacopy) or ``"allele_specific"``
    (ascat2, ascat3, ascatngs). Empty list when the project has none.
    """
    if kind != "gene_level" and kind not in _SEGMENT_KINDS:
        raise ValueError(f"kind must be 'gene_level' or one of {sorted(_SEGMENT_KINDS)}")
    store = _store(project, root)
    if "copy_number" not in store.listdir():
        return []
    return _cn_workflows(store, kind)


def copy_number_samples(project: str, workflow: str | None = None, root: str | None = None) -> pd.DataFrame:
    """copy_number/gene_level_<workflow>_samples: matrix column -> case_submitter_id, file_id, ...

    ``sample_type``, ``sample_submitter_id`` and ``case_submitter_id`` describe the matrix
    column's own aliquot, which is the tumor member of the tumor/normal pair the file was
    called on (in TCGA-BLCA all 392 ASCAT3 columns are primary tumor barcodes and all 392
    rows say Primary Tumor). For the paired workflows ``aliquot_submitter_id`` holds both
    aliquots of the pair pipe-joined, in no fixed order. Tables built before 2026-09-13
    reported the pair's first aliquot as ``sample_type`` instead. Every matrix column is a
    single tumor aliquot.
    """
    store = _store(project, root)
    wf = _pick_workflow(_cn_workflows(store, "gene_level"), workflow, "gene-level copy number", project)
    df = _read(store, f"copy_number/gene_level_{wf}_samples.parquet")
    df.attrs["workflow"] = wf
    return df


def copy_number(
    project: str,
    workflow: str | None = None,
    genes: str | Iterable[str] | None = None,
    drop_par_y: bool = True,
    root: str | None = None,
) -> pd.DataFrame:
    """Gene-level copy number (float total copy number, NaN where not called) as genes x tumor aliquots.

    ``workflow`` defaults to ascat3 when present, else the first available; the choice
    is in ``df.attrs["workflow"]``. The workflows are different pipelines and are not
    interchangeable. ``genes`` is pushed down as a gene_name row filter. The index is
    gene_name; ``df.attrs["gene_id"]`` maps it to the Ensembl id and
    ``df.attrs["gene_location"]`` to (chromosome, start, end). As in
    ``expression()``, ``_PAR_Y`` ids are dropped and duplicated names keep the first row.
    """
    store = _store(project, root)
    wf = _pick_workflow(_cn_workflows(store, "gene_level"), workflow, "gene-level copy number", project)
    gene_list = _as_list(genes)
    filters: Any = [("gene_name", "in", gene_list)] if gene_list is not None else None
    df = _read(store, f"copy_number/gene_level_{wf}.parquet", filters=filters)
    if gene_list is not None and df.empty:
        raise ValueError(f"none of {gene_list} match a gene_name in copy_number/gene_level_{wf}")
    out = _gene_matrix(df, gene_list, drop_par_y, False)
    out.attrs.update({"project": project, "workflow": wf})
    return out


def segments(project: str, kind: str = "segments", workflow: str | None = None, root: str | None = None) -> pd.DataFrame:
    """Long segment tables: ``kind`` is ``"segments"`` (Segment_Mean log2 ratios; dnacopy or
    gatk4_cnv), ``"masked"`` (masked_segments_dnacopy, germline CNV removed) or
    ``"allele_specific"`` (ASCAT Copy_Number / Major / Minor per segment). ``workflow``
    defaults to the preferred one present (ascat3, dnacopy, ...) and is recorded in attrs.
    """
    if kind not in _SEGMENT_KINDS:
        raise ValueError(f"kind must be one of {sorted(_SEGMENT_KINDS)}, got {kind!r}")
    store = _store(project, root)
    wf = _pick_workflow(_cn_workflows(store, kind), workflow, f"{kind} segment", project)
    df = _read(store, f"copy_number/{_SEGMENT_KINDS[kind]}_{wf}.parquet")
    df.attrs.update({"project": project, "kind": kind, "workflow": wf})
    return df


# -------------------------------------------------------------- miRNA and RPPA


def mirna_samples(project: str, root: str | None = None) -> pd.DataFrame:
    """mirna/samples: matrix column -> case_submitter_id, sample_type, file_id, ..."""
    return _read(_store(project, root), "mirna/samples.parquet")


def mirna(project: str, kind: str = "rpm", sample_type: str | Iterable[str] | None = None, root: str | None = None) -> pd.DataFrame:
    """BCGSC miRNA profiling as miRNA_ID x aliquots: ``"rpm"`` (reads per million miRNA
    mapped) or ``"counts"`` (read_count). ``sample_type`` selects columns through
    mirna/samples so only those columns are read.
    """
    if kind not in _MIRNA_KINDS:
        raise ValueError(f"kind must be one of {sorted(_MIRNA_KINDS)}, got {kind!r}")
    store = _store(project, root)
    cols = _select_columns(mirna_samples(project, root), sample_type, "miRNA")
    df = _read(store, _MIRNA_KINDS[kind], columns=["miRNA_ID", *cols]).set_index("miRNA_ID")
    df.attrs.update({"project": project, "kind": kind})
    return df


def mirna_isoforms(project: str, root: str | None = None) -> pd.DataFrame:
    """mirna/isoforms: long table (file_id, aliquot_submitter_id, miRNA_ID, isoform_coords,
    read_count, reads_per_million_miRNA_mapped, cross-mapped, miRNA_region). Millions of rows.
    """
    return _read(_store(project, root), "mirna/isoforms.parquet")


def rppa_samples(project: str, root: str | None = None) -> pd.DataFrame:
    """protein/samples: RPPA matrix column (a 19-character portion id) -> case_submitter_id, sample_type."""
    return _read(_store(project, root), "protein/samples.parquet")


def rppa_antibodies(project: str, root: str | None = None) -> pd.DataFrame:
    """protein/antibodies: AGID, peptide_target, lab_id, catalog_number, set_id."""
    return _read(_store(project, root), "protein/antibodies.parquet")


def rppa(project: str, sample_type: str | Iterable[str] | None = None, root: str | None = None) -> pd.DataFrame:
    """RPPA protein expression as peptide_target x samples (float). The antibody id of
    each row is in ``df.attrs["agid"]``. Columns are RPPA portion ids (``TCGA-GV-A3QG-01A-21``),
    mapped to cases by ``rppa_samples``.
    """
    store = _store(project, root)
    cols = _select_columns(rppa_samples(project, root), sample_type, "RPPA")
    df = _read(store, "protein/rppa_protein_expression.parquet", columns=["AGID", "peptide_target", *cols])
    out = df.drop(columns=["AGID", "peptide_target"]).set_index(pd.Index(df["peptide_target"], name="peptide_target"))
    out.attrs["agid"] = _Frozen(zip(df["peptide_target"].tolist(), df["AGID"].tolist()))
    out.attrs["project"] = project
    return out


# ------------------------------------------------------------------ methylation


def methylation_platforms(project: str, root: str | None = None) -> list[str]:
    """Array platforms with a betas table, e.g. ``["human_methylation_450"]``. Empty when the
    project was built without methylation (most projects: the matrices are large).
    """
    store = _store(project, root)
    if "methylation" not in store.listdir():
        return []
    prefix = "betas_"
    return sorted(
        f[len(prefix) : -len(".parquet")]
        for f in store.listdir("methylation")
        if f.startswith(prefix) and f.endswith(".parquet") and not f.endswith("_samples.parquet")
    )


def _pick_platform(project: str, platform: str | None, root: str | None) -> str:
    plats = methylation_platforms(project, root)
    if not plats:
        raise NotMirroredError(f"{project} has no methylation tables (methylation is built only on request)")
    if platform is None:
        if len(plats) == 1:
            return plats[0]
        rep = report(project, root).get("tables", {})
        return max(plats, key=lambda p: int(rep.get(f"methylation/betas_{p}", {}).get("columns", 0)))
    if platform not in plats:
        raise ValueError(f"unknown methylation platform {platform!r} for {project}; available: {plats}")
    return platform


def methylation_samples(project: str, platform: str | None = None, root: str | None = None) -> pd.DataFrame:
    """methylation/betas_<platform>_samples: matrix column -> case_submitter_id, sample_type, ..."""
    plat = _pick_platform(project, platform, root)
    df = _read(_store(project, root), f"methylation/betas_{plat}_samples.parquet")
    df.attrs["platform"] = plat
    return df


def methylation(
    project: str,
    platform: str | None = None,
    probes: str | Iterable[str] | None = None,
    columns: str | Iterable[str] | None = None,
    sample_type: str | Iterable[str] | None = None,
    root: str | None = None,
) -> pd.DataFrame:
    """SeSAMe beta values as probe_id x aliquots (float32), for one array platform.

    ``platform`` defaults to the only one present, or the one with the most columns
    (recorded in ``df.attrs["platform"]``). ``probes`` (cg ids) becomes a parquet row
    filter: the file is written in 50,000-probe row groups, so pyarrow reads only
    the row groups whose statistics can contain those ids and materialises only the
    matching rows. ``columns`` is a list of aliquot barcodes to read, and
    ``sample_type`` picks columns through the samples table; both become a parquet
    column projection. Memory is therefore bounded by rows kept x columns kept x 4
    bytes plus one row group of the projected columns (50,000 x columns x 4 bytes).

    With neither ``probes`` nor ``columns`` nor ``sample_type`` the whole matrix is
    read: 486,427 probes x 45 aliquots is 116 MB in pandas (TCGA-CHOL 450k); a 1,200-aliquot
    450k project is about 2.3 GB, so pass probes or columns for large projects.
    """
    plat = _pick_platform(project, platform, root)
    store = _store(project, root)
    rel = f"methylation/betas_{plat}.parquet"
    smp = methylation_samples(project, plat, root)
    col_list = _as_list(columns)
    if col_list is not None:
        known = set(store.parquet_schema(rel))
        bad = [c for c in col_list if c not in known]
        if bad:
            raise ValueError(f"unknown methylation columns {bad[:5]}{'...' if len(bad) > 5 else ''}; see methylation_samples()")
    if sample_type is not None:
        by_type = _select_columns(smp, sample_type, "methylation")
        col_list = [c for c in col_list if c in set(by_type)] if col_list is not None else by_type
    probe_list = _as_list(probes)
    filters: Any = [("probe_id", "in", probe_list)] if probe_list is not None else None
    read_cols = None if col_list is None else ["probe_id", *col_list]
    df = _read(store, rel, columns=read_cols, filters=filters)
    if probe_list is not None and df.empty:
        raise ValueError(f"none of {probe_list[:5]} are probe_ids on {plat}")
    out = df.set_index("probe_id")
    out.attrs.update({"project": project, "platform": plat})
    return out


# ---------------------------------------------------------------------- helpers


def case_for_barcode(barcode: str) -> str:
    """Case (patient) id from an aliquot or sample barcode.

    TCGA barcodes: the first 12 characters (``TCGA-DK-A6B2-01A-11R-A30C-07`` -> ``TCGA-DK-A6B2``).
    TARGET barcodes: the first 16 (``TARGET-30-PAAPFA-01A-01R`` -> ``TARGET-30-PAAPFA``).
    Anything else is returned unchanged: CPTAC, HCMI, MMRF and other programs use
    ids without a positional structure, so join them to cases through
    ``aliquots()`` (aliquot_submitter_id -> case_submitter_id) or the ``*_samples`` tables.
    """
    if barcode.startswith("TCGA-"):
        return barcode[:12]
    if barcode.startswith("TARGET-"):
        return barcode[:16]
    return barcode


def _manifest_row(project: str, file_id: str, root: str | None) -> pd.Series:
    df = manifest(project, root=root)
    hit = df[df["file_id"] == file_id]
    if hit.empty:
        raise ValueError(f"file_id {file_id!r} is not in the {project} manifest ({len(df)} files)")
    return hit.iloc[0]


def _stream_s3(bucket: str, key: str) -> Iterable[bytes]:
    import boto3

    body = boto3.client("s3").get_object(Bucket=bucket, Key=key)["Body"]
    for chunk in body.iter_chunks(1 << 20):
        yield bytes(chunk)


def _stream_https(url: str) -> Iterable[bytes]:
    from urllib.request import urlopen

    with urlopen(url) as resp:  # noqa: S310 (fixed https host from the manifest)
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                return
            yield chunk


def _s3_unreadable(exc: Exception) -> bool:
    """True for the errors that mean 'fall back to HTTPS' rather than 'a bug'."""
    from botocore.exceptions import BotoCoreError, ClientError, NoCredentialsError

    if isinstance(exc, NoCredentialsError):
        return True
    if isinstance(exc, ClientError):
        err = exc.response.get("Error", {})
        code = str(err.get("Code", ""))
        status = exc.response.get("ResponseMetadata", {}).get("HTTPStatusCode")
        return code in {"AccessDenied", "NoSuchKey", "NoSuchBucket", "403", "404", "InvalidAccessKeyId", "ExpiredToken"} or status in (
            403,
            404,
        )
    return isinstance(exc, BotoCoreError)


def fetch_raw(
    project_or_manifest_row: str | pd.Series,
    file_id: str | None = None,
    dest: str | None = None,
    root: str | None = None,
) -> bytes | str:
    """Download one raw GDC file named by a manifest row (or by ``project`` + ``file_id``).

    Reads ``s3://<s3_bucket>/<s3_key>`` with a signed default boto3 client (the newer
    GDC objects are SSE-KMS and refuse anonymous reads; any AWS principal works), and
    falls back to ``gdc_download_url`` over HTTPS when S3 answers AccessDenied, 403 or
    404 or no credentials are configured. The bytes are md5-checked against the
    manifest's ``md5sum`` (NciError on mismatch). Returns the bytes, or the path
    written when ``dest`` is given (a directory gets ``file_name`` appended).
    """
    if isinstance(project_or_manifest_row, pd.Series):
        row = project_or_manifest_row
    else:
        if file_id is None:
            raise ValueError("pass a manifest row, or a project id and a file_id")
        row = _manifest_row(project_or_manifest_row, file_id, root)
    for col in ("s3_bucket", "s3_key", "md5sum", "file_name"):
        if col not in row.index:
            raise ValueError(f"manifest row lacks {col!r}; pass a row of manifest()")

    def _pull(chunks: Iterable[bytes]) -> tuple[bytes | None, str, str | None]:
        h = hashlib.md5()
        path: str | None = None
        if dest is not None:
            path = os.path.join(dest, str(row["file_name"])) if os.path.isdir(dest) else dest
            with open(path, "wb") as fh:
                for c in chunks:
                    h.update(c)
                    fh.write(c)
            return None, h.hexdigest(), path
        buf = bytearray()
        for c in chunks:
            h.update(c)
            buf.extend(c)
        return bytes(buf), h.hexdigest(), None

    source = "s3"
    try:
        data, digest, path = _pull(_stream_s3(str(row["s3_bucket"]), str(row["s3_key"])))
    except Exception as exc:
        url = row.get("gdc_download_url") if hasattr(row, "get") else None
        if not _s3_unreadable(exc) or not url:
            raise
        source = "https"
        data, digest, path = _pull(_stream_https(str(url)))
    expected = str(row["md5sum"]).lower()
    if expected and digest != expected:
        if path is not None and os.path.exists(path):
            os.remove(path)
        raise NciError(f"md5 mismatch for {row['file_name']} via {source}: got {digest}, manifest says {expected}")
    if path is not None:
        return path
    assert data is not None
    return data
