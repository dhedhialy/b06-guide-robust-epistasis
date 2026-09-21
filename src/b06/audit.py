"""Stage 0 design audit: guide multiplicity, matched controls, coverage."""

from __future__ import annotations

from collections import defaultdict

import pandas as pd

from .io import CONTROL_CLASSES, condition


def singles_for(df: pd.DataFrame, guide: str) -> pd.DataFrame:
    """Rows where `guide` is the only functional guide (inert partner)."""
    inert = lambda s: s.isin(CONTROL_CLASSES)
    return df[
        ((df["guide_a"] == guide) & inert(df["class_b"]))
        | ((df["guide_b"] == guide) & inert(df["class_a"]))
    ]


def guide_usage(df: pd.DataFrame) -> pd.DataFrame:
    """Per-guide table: target gene, rows, single/double/control instances."""
    rows = []
    for pos in ("a", "b"):
        sub = df.rename(columns={f"guide_{pos}": "guide", f"gene_{pos}": "gene"})
        for guide, g in sub.groupby("guide"):
            c = g["gene"].iloc[0]
            cond = g.apply(condition, axis=1)
            rows.append(
                {
                    "guide": guide,
                    "gene": c,
                    "rows": len(g),
                    "singles": int((cond == "single").sum()),
                    "doubles": int((cond == "double").sum()),
                    "controls": int((cond == "control").sum()),
                }
            )
    out = pd.DataFrame(rows).groupby(["guide", "gene"], as_index=False).sum()
    out["has_single_control"] = out["singles"] > 0
    out["bridges_single_double"] = (out["singles"] > 0) & (out["doubles"] > 0)
    return out.sort_values("rows", ascending=False).reset_index(drop=True)


def gene_guide_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Guides per gene (unique across both positions)."""
    rows = []
    for pos in ("a", "b"):
        rows.append(df[[f"gene_{pos}", f"guide_{pos}"]].rename(
            columns={f"gene_{pos}": "gene", f"guide_{pos}": "guide"}))
    g = pd.concat(rows).drop_duplicates()
    out = g.groupby("gene")["guide"].nunique().rename("n_guides").reset_index()
    return out.sort_values("n_guides", ascending=False).reset_index(drop=True)


def gene_pair_guides(df: pd.DataFrame) -> pd.DataFrame:
    """Per gene pair: distinct guides on each side, observed guide pairs."""
    rows = []
    for gene_a, sub_a in df.groupby("gene_a"):
        for gene_b, g in sub_a.groupby("gene_b"):
            pair = tuple(sorted([gene_a, gene_b]))
            rows.append(
                {
                    "gene_a": pair[0],
                    "gene_b": pair[1],
                    "guides_a": g["guide_a"].nunique(),
                    "guides_b": g["guide_b"].nunique(),
                    "observed_pairs": len(g),
                    "max_pairs": g["guide_a"].nunique() * g["guide_b"].nunique(),
                }
            )
    out = pd.DataFrame(rows)
    out["cross_product_frac"] = out["observed_pairs"] / out["max_pairs"].clip(lower=1)
    return out


def audit_summary(df: pd.DataFrame) -> dict:
    """Machine-readable Stage 0 record for GUIDE_MULTIPLICITY_LOCK."""
    gu = guide_usage(df)
    gpu = gene_pair_guides(df)
    genes = gene_guide_counts(df)
    nnt = df[(df["class_a"] == "non-targeting") & (df["class_b"] == "non-targeting")]
    bridges = gu[gu["bridges_single_double"]]
    return {
        "guide_pairs_total": len(df),
        "unique_guides": df["guide_a"].nunique() + df["guide_b"].nunique(),
        "genes": len(genes),
        "triaged": int(df[df.vector_class.str.contains("non-targeting")].shape[0]),
        "genes_with_2plus_guides": int((genes["n_guides"] >= 2).sum()),
        "gene_pairs_observed": len(gpu),
        "gene_pairs_full_cross_product": int((gpu["cross_product_frac"] == 1).sum()),
        "gene_pairs_with_doubles": int((gpu["observed_pairs"] > 0).sum()),
        "guides_with_matched_single": int((gu["has_single_control"]).sum()),
        "guides_bridging_single_and_double": len(bridges),
        "nontargeting_control_pairs": len(nnt),
        "control_phenotype_median": float(nnt["FC_500x"].median()),
    }


def matched_guide_pairs(df: pd.DataFrame, mu: float) -> pd.DataFrame:
    """Guide pairs (a, b) with real guide-level singles on both sides.

    A guide's 'single' context is its rows paired with an inert control guide
    (non-targeting or intergenic); its 'double' context is (a, b) together.
    Returns one row per (guide_a, guide_b) with all scalars B06 needs.
    Indexed the first pass for scale (the COLO1 screen has ~79K rows).
    """
    inert = lambda s: s in CONTROL_CLASSES

    singles: dict[str, list] = defaultdict(list)
    doubles: dict[tuple, list] = defaultdict(list)
    gene_of: dict[str, str] = {}
    for _, r in df.iterrows():
        ga, gb, ca, cb = r["guide_a"], r["guide_b"], r["class_a"], r["class_b"]
        gene_of.setdefault(ga, r["gene_a"])
        gene_of.setdefault(gb, r["gene_b"])
        if inert(ca) and inert(cb):
            continue
        if inert(ca):
            singles[gb].append(r["FC_500x"])
        elif inert(cb):
            singles[ga].append(r["FC_500x"])
        else:
            key = (ga, gb) if ga <= gb else (gb, ga)
            doubles[key].append(r["FC_500x"])

    roots = {g for g, v in singles.items() if v}
    rows = [
        {
            "gene_a": gene_of[ga],
            "gene_b": gene_of[gb],
            "guide_a": ga,
            "guide_b": gb,
            "m_ab": float(sum(v) / len(v)),
            "m_a0": float(sum(singles[ga]) / len(singles[ga])),
            "m_0b": float(sum(singles[gb]) / len(singles[gb])),
            "n_double": len(v),
            "n_single_a": len(singles[ga]),
            "n_single_b": len(singles[gb]),
        }
        for (ga, gb), v in doubles.items()
        if ga in roots and gb in roots
    ]
    out = pd.DataFrame(rows)
    if not out.empty:
        out["mu"] = mu
    return out