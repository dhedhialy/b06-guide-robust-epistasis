"""Burgold pilot-library parsing (Supplementary Data 1, MOESM3) and the
COLO1 100K screen (DualGuide_COLO1, Zenodo 10.5281/zenodo.17191951)."""

from __future__ import annotations

import numpy as np
import pandas as pd

PHENOTYPE = "FC_500x"
CONTROL_CLASSES = {"non-targeting", "intergenic"}

SHEET = "Library"
COLUMNS = {
    "ID": "guide_pair_id",
    "Notes": "notes",
    "Scaffold": "scaffold",
    "sgRNA1_WGE_ID": "guide_a",
    "sgRNA1_Approved_Symbol": "gene_a",
    "sgRNA1_class": "class_a",
    "sgRNA2_WGE_ID": "guide_b",
    "sgRNA2_Approved_Symbol": "gene_b",
    "sgRNA2_class": "class_b",
    "vector_class": "vector_class",
    "FC_500x": PHENOTYPE,
}

# COLO1 count columns: 3 plasmid replicates (reference arm) + 3 HT29 endpoint
# replicates (pairwise correlation separates the two triplets cleanly; the
# PositiveControls deplete under +log2(HT29/plasmid), fixing the orientation).
COLO1_REF = ["SIDM00136_CPID2437", "SIDM00136_CPID2440", "SIDM00136_CPID2443"]
COLO1_END = ["SIDM00136_CPID1020", "SIDM00136_CPID1023", "SIDM00136_CPID1026"]
COLO1_COMBO_NOTES = ("AnchorCombinations", "LibraryCombinations", "GIControlsCombinations")
COLO1_SINGLE_NOTES = ("AnchorSingletons", "LibrarySingletons", "GIControlsSingletons")


def load_design(path: str) -> pd.DataFrame:
    """Load the pilot design table; each row is a guide-pair construct."""
    df = pd.read_excel(path, sheet_name=SHEET)
    df = df[list(COLUMNS)].rename(columns=COLUMNS).copy()
    df["gene_a"] = df["gene_a"].fillna("").astype(str).str.strip()
    df["gene_b"] = df["gene_b"].fillna("").astype(str).str.strip()
    df["guide_a"] = df["guide_a"].astype(str).str.strip()
    df["guide_b"] = df["guide_b"].astype(str).str.strip()
    return df.dropna(subset=[PHENOTYPE]).reset_index(drop=True)


def load_colo1(path: str, lfc_plm: str | None = None) -> pd.DataFrame:
    """Load the COLO1 screen as the same per-construct schema as the pilot.

    Each row is one guide-pair construct. Phenotype is FC = log2(HT29/plasmid)
    derived from the two replicate triplets, or, when ``lfc_plm`` points at
    the authors' ``colo1_plm.rds`` (normalized endpoint LFC model), the
    authors' own normalized LFC (mean of the 3 endpoint replicates), matched
    by construct ID. Inert (control) guides are flagged like the pilot's
    non-targeting class so the rest of the pipeline is unchanged.
    """
    df = pd.read_csv(path, sep="\t", dtype={"sgRNA1_ID": str, "sgRNA2_ID": str, "Gene1": str, "Gene2": str})
    for col in ("sgRNA1_ID", "sgRNA2_ID", "Gene1", "Gene2"):
        df[col] = df[col].astype(str).str.strip()
    df["_id"] = df["ID"].astype(str).str.strip()

    if lfc_plm is not None:
        plm = _read_colo1_plm(lfc_plm)
        df = df.merge(plm, left_on="_id", right_index=True, how="left")
        df["_col1_fc"] = df[["__plm_0", "__plm_1", "__plm_2"]].mean(axis=1)
        ok = df["_col1_fc"].notna()
    else:
        ref = df[COLO1_REF].to_numpy(float)
        end = df[COLO1_END].to_numpy(float)
        ok = (ref.sum(axis=1) > 0) & (end.sum(axis=1) > 0)
        rmean = np.exp(np.mean(np.log(ref.clip(min=1)), axis=1))
        emean = np.exp(np.mean(np.log(end.clip(min=1)), axis=1))
        df["_col1_fc"] = np.log2(emean / rmean)
        ok = ok & np.isfinite(df["_col1_fc"])
    df = df[ok].copy()

    combo = set(
        df.loc[df["Note"].isin(COLO1_COMBO_NOTES), ["sgRNA1_ID", "sgRNA2_ID"]]
        .stack()
        .drop_duplicates()
    )
    neutral = set(
        df.loc[df["Note"] == "NegativeControls", ["sgRNA1_ID", "sgRNA2_ID"]]
        .stack()
        .drop_duplicates()
    )
    inert = set(neutral)
    # A single-note row pairs a functional guide with an inert partner: the
    # partner (the guide NOT in the combo set) becomes the matched-single
    # control for that functional guide, so it must be marked non-targeting.
    for _, r in df[df["Note"].isin(COLO1_SINGLE_NOTES)].iterrows():
        functional = r["sgRNA1_ID"] if r["sgRNA1_ID"] in combo else r["sgRNA2_ID"]
        if functional in combo:
            inert.add(r["sgRNA1_ID"] if r["sgRNA2_ID"] == functional else r["sgRNA2_ID"])

    out = pd.DataFrame(
        {
            "guide_pair_id": df["ID"].astype(str).str.strip(),
            "notes": df["Note"],
            "scaffold": "COLO1",
            "vector_class": df["Note"],
            "guide_a": df["sgRNA1_ID"],
            "gene_a": df["Gene1"],
            "class_a": np.where(df["sgRNA1_ID"].isin(inert), "non-targeting", "gene"),
            "guide_b": df["sgRNA2_ID"],
            "gene_b": df["Gene2"],
            "class_b": np.where(df["sgRNA2_ID"].isin(inert), "non-targeting", "gene"),
            PHENOTYPE: df["_col1_fc"],
        }
    )
    return out.reset_index(drop=True)


def _read_colo1_plm(path: str) -> pd.DataFrame:
    """Read the authors' normalized-endpoint LFC model (colo1_plm.rds).

    Rows are construct IDs (GI000000001 ...) with the 3 endpoint-replicate
    normalized LFC columns. Returns a frame indexed by stripped construct ID
    with three renamed columns __plm_0/1/2 (their mean is the phenotype).
    Only importable when pyreadr is installed (an optional dependency).
    """
    import pyreadr

    plm = pyreadr.read_r(path)[None]
    plm = plm.rename(columns={c: f"__plm_{i}" for i, c in enumerate(plm.columns)})
    plm.index = plm.index.astype(str).str.strip()
    return plm


def is_control_class(cls: str) -> bool:
    return cls in CONTROL_CLASSES


def condition(row: pd.Series) -> str:
    """Row-level condition: control / single_A / single_B / double."""
    ca, cb = is_control_class(row["class_a"]), is_control_class(row["class_b"])
    if ca and cb:
        return "control"
    if ca:
        return "single"
    if cb:
        return "single"
    return "double"


def mu_control(df: pd.DataFrame) -> float:
    """Reference phenotype: median over non-targeting x non-targeting pairs."""
    nt = df[
        (df["class_a"] == "non-targeting") & (df["class_b"] == "non-targeting")
    ][PHENOTYPE]
    return float(nt.median()) if not nt.empty else 0.0