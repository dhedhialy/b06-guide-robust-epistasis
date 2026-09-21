"""Gene-level identified sets (intersection) and guide-level bootstrapped
sampling intervals (Sec 6.1, Sec 8)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .bounds import classify, symmetric_bounds


def gene_interval(contrasts: pd.DataFrame, rho: float) -> pd.DataFrame:
    """Intersect guide-pair intervals per gene pair (sharp gene-level set)."""
    rows = []
    for (ga, gb), g in contrasts.groupby(["gene_a", "gene_b"]):
        lows, ups = [], []
        for _, r in g.iterrows():
            res = symmetric_bounds(r["m_a0"], r["m_0b"], r["m_ab"], r["mu"], rho, noise=r["noise"])
            lows.append(res.lower)
            ups.append(res.upper)
        lo, up = max(lows), min(ups)
        rows.append(
            {
                "gene_a": ga,
                "gene_b": gb,
                "lower": lo,
                "upper": up,
                "sign": classify(lo, up),
                "n_guide_pairs": len(g),
            }
        )
    return pd.DataFrame(rows)


def block_bootstrap_interval(contrasts: pd.DataFrame, rho: float, n_iter: int, seed: int):
    """Guide-level block bootstrap of the gene-level interval per gene pair.

    Resamples guide pairs within each gene pair. Per gene pair returns the
    mean interval endpoints over bootstrap iterations and the fraction of
    iterations whose intersection interval has a definitive sign call.
    """
    rng = np.random.default_rng(seed)
    deltas = contrasts["delta"].to_numpy()
    widths = rho * (
        np.abs(contrasts["m_a0"].to_numpy() - contrasts["mu"].to_numpy())
        + np.abs(contrasts["m_0b"].to_numpy() - contrasts["mu"].to_numpy())
    ) + contrasts["noise"].to_numpy()
    lo_r = deltas - widths
    up_r = deltas + widths

    mean_lo, mean_up, robust_frac = {}, {}, {}
    for (ga, gb), idx in contrasts.groupby(["gene_a", "gene_b"]).groups.items():
        pos = contrasts.index.get_indexer(idx)
        lo = lo_r[pos]
        up = up_r[pos]
        n = len(pos)
        lo_perm, up_perm, calls = np.empty(n_iter), np.empty(n_iter), 0
        for j in range(n_iter):
            s = rng.integers(0, n, n)
            lo_perm[j] = lo[s].max()
            up_perm[j] = up[s].min()
            calls += classify(lo_perm[j], up_perm[j]) != "UNRESOLVED"
        mean_lo[(ga, gb)] = float(lo_perm.mean())
        mean_up[(ga, gb)] = float(up_perm.mean())
        robust_frac[(ga, gb)] = calls / n_iter

    gi = gene_interval(contrasts, rho).sort_values(["gene_a", "gene_b"]).reset_index(drop=True)
    order = list(mean_lo.keys())
    return (
        [mean_lo[k] for k in order],
        [mean_up[k] for k in order],
        [robust_frac[k] for k in order],
    )