"""scigantic-nci: analysis-ready NCI cancer data (GDC projects, IDC collections) in one import.

    import scigantic_nci as nci

    nci.gdc.projects()                                   # what is mirrored
    clin = nci.gdc.clinical("TCGA-BRCA")                 # one row per case
    tpm = nci.gdc.expression("TCGA-BRCA", genes=["TP53", "ESR1"], sample_type="Primary Tumor")
    maf = nci.gdc.mutations("TCGA-BRCA", genes=["TP53"])

    nci.idc.collections()
    s = nci.idc.series("tcga_lihc", modality="CT")
    nci.idc.pull_series(s.iloc[0], "/tmp/ct")             # anonymous S3, no account

Reads a mounted Scigantic archive when one is present, otherwise the public
mirror buckets anonymously. See README for the tables and their provenance.
"""
from __future__ import annotations

from . import gdc, idc
from ._store import NciError, NotMirroredError
from ._version import __version__

__all__ = ["gdc", "idc", "NciError", "NotMirroredError", "__version__"]
