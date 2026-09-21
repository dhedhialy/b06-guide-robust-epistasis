"""Stage 1 scout: matched-guide contrasts, rho_star, symmetric bounds."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .bounds import symmetric_bounds

RHO_GRID = (0.05, 0.10, 0.20, 0.30)

# Guide-idiosyncrasy floor: one within-pair MAD of guide-pair deltas. The
# drift envelope rho*(|d_a|+|d_b|) covers context-rate misspecification but
# not guide-level spread; calls without the floor do not reproduce under
# guide-identity holdout.
NOISE_K = 1.0


def noise_floor(deltas) -> float:
    """Robust scale of guide-pair deltas within one gene pair (MAD, delta units)."""
    x = np.asarray(deltas, dtype=float)
    if len(x) < 2:
        return 0.0
    return NOISE_K * 1.4826 * float(np.median(np.abs(x - np.median(x))))


def guide_pair_contrasts(matched: pd.DataFrame) -> pd.DataFrame:
    """Delta, rho_star and symmetric bounds for every matched guide pair."""
    if matched.empty:
        return pd.DataFrame()
    rows = []
    for _, r in matched.iterrows():
        res = symmetric_bounds(r["m_a0"], r["m_0b"], r["m_ab"], r["mu"], rho=0.0)
        gene_a, gene_b = sorted([r["gene_a"], r["gene_b"]])
        base = {
            "gene_a": gene_a,
            "gene_b": gene_b,
            "guide_a": r["guide_a"],
            "guide_b": r["guide_b"],
            "m_a0": r["m_a0"],
            "m_0b": r["m_0b"],
            "m_ab": r["m_ab"],
            "mu": r["mu"],
            "delta": res.lower,  # rho=0 collapses to delta
            "rho_star": res.rho_star,
            "n_double": int(r["n_double"]),
            "n_single_a": int(r["n_single_a"]),
            "n_single_b": int(r["n_single_b"]),
        }
        rows.append(base)
    ct = pd.DataFrame(rows)
    ct["noise"] = ct.groupby(["gene_a", "gene_b"])["delta"].transform(noise_floor)
    bound_rows = []
    for _, r in ct.iterrows():
        base = {}
        for rho in RHO_GRID:
            res = symmetric_bounds(r["m_a0"], r["m_0b"], r["m_ab"], r["mu"], rho, noise=r["noise"])
            base[f"lower_{rho}"] = res.lower
            base[f"upper_{rho}"] = res.upper
            base[f"sign_{rho}"] = res.sign
        bound_rows.append(base)
    out = ct.join(pd.DataFrame(bound_rows))
    return out


def gene_pair_calls(contrasts: pd.DataFrame, rho: float) -> pd.DataFrame:
    """Aggregate guide-pair calls to gene level.

    A gene pair is robust-POSITIVE if every guide pair with a definitive call
    agrees and at least one is callable; UNRESOLVED otherwise.
    """
    if contrasts.empty:
        return pd.DataFrame()
    cols = ["gene_a", "gene_b", f"sign_{rho}"]
    g = contrasts[cols].groupby(["gene_a", "gene_b"])
    out = g.agg(
        n_guide_pairs=(f"sign_{rho}", "count"),
        n_positive=(f"sign_{rho}", lambda s: (s == "POSITIVE").sum()),
        n_negative=(f"sign_{rho}", lambda s: (s == "NEGATIVE").sum()),
        n_unresolved=(f"sign_{rho}", lambda s: (s == "UNRESOLVED").sum()),
    ).reset_index()
    resolved = out["n_positive"].astype(int) + out["n_negative"].astype(int)
    out["gene_sign"] = "UNRESOLVED"
    out.loc[(out["n_positive"] > 0) & (out["n_negative"] == 0) & (resolved == out["n_guide_pairs"]), "gene_sign"] = "POSITIVE"
    out.loc[(out["n_negative"] > 0) & (out["n_positive"] == 0) & (resolved == out["n_guide_pairs"]), "gene_sign"] = "NEGATIVE"
    return out


def guide_conflicts(contrasts: pd.DataFrame, rho: float) -> pd.DataFrame:
    """Gene pairs whose guides imply opposite signs (Sec 21.3)."""
    if contrasts.empty:
        return pd.DataFrame()
    g = contrasts.groupby(["gene_a", "gene_b"])[f"sign_{rho}"]
    out = g.agg([("signs", lambda s: sorted(set(s))), ("n_guide_pairs", "count")]).reset_index()
    return out[(out["signs"].apply(len) > 1)]


def gene_scores(contrasts: pd.DataFrame, rho: float) -> pd.DataFrame:
    """Continuous, rankable gene-level scores.

    ``robust_score`` is the signed, envelope-normalized mean delta
    (``mean_delta / (rho*(mean|d_a|+mean|d_b|) + noise)``): magnitude is how
    deep inside the certified region the pair sits, sign is the interaction
    direction. ``conv_score`` is the raw mean delta (standard dLFC). The
    certificate's selectivity flag is ``certified = |mean_delta| > envelope``.
    """
    if contrasts.empty:
        return pd.DataFrame()
    g = contrasts.groupby(["gene_a", "gene_b"])
    out = g.agg(
        n_guide_pairs=("delta", "count"),
        mean_delta=("delta", "mean"),
        mean_abs_single_a=("m_a0", lambda s: (s - contrasts.loc[s.index, "mu"]).abs().mean()),
        mean_abs_single_b=("m_0b", lambda s: (s - contrasts.loc[s.index, "mu"]).abs().mean()),
        noise=("noise", "median"),
    ).reset_index()
    env = rho * (out["mean_abs_single_a"] + out["mean_abs_single_b"]) + out["noise"]
    out["envelope"] = env
    out["robust_score"] = out["mean_delta"] / env
    out["conv_score"] = out["mean_delta"]
    out["certified"] = out["mean_delta"].abs() > env
    return out