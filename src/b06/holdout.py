"""Stage 2: guide-identity holdout — do robust calls predict held-out guide
behavior better than conventional point estimates? (Sec 10, Module D)"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .bounds import classify, symmetric_bounds
from .contrasts import NOISE_K


def _noise_floor_calibration(deltas, noise_k: float = NOISE_K) -> float:
    """Split-level guide-idiosyncrasy floor computed from calibration only."""
    x = np.asarray(deltas, dtype=float)
    if len(x) < 2:
        return 0.0
    return noise_k * 1.4826 * float(np.median(np.abs(x - np.median(x))))


def _robust_call(rho_bands, rho):
    """Gene-level call from the intersection of calibration intervals."""
    lows = [b.lower for b in rho_bands]
    ups = [b.upper for b in rho_bands]
    return classify(max(lows), min(ups))


def _conventional_call(deltas):
    return classify(sum(deltas) / len(deltas), sum(deltas) / len(deltas))


def evaluate_guide_holdout(contrasts: pd.DataFrame, rho: float, min_calibration: int = 2, noise_k: float = NOISE_K) -> pd.DataFrame:
    """LOOCV by guide-pair identity per gene pair.

    For each held-out guide pair: calibrate the gene-level interval on the
    remaining pairs at drift rho (robust) vs the mean-delta call (conventional),
    then score both against the held-out pair's own delta sign.

    Returns one row per gene pair with concordance metrics at matched
    coverage (only evaluated where the robust method committed to a sign).
    """
    rows = []
    for (ga, gb), g in contrasts.groupby(["gene_a", "gene_b"]):
        k = len(g)
        if k < min_calibration + 1:
            continue
        robust_hits, robust_calls, conv_hits, conv_calls = 0, 0, 0, 0
        for idx in g.index:
            held = g.loc[idx]
            cal = g.drop(idx)
            floor = _noise_floor_calibration(cal["delta"].tolist(), noise_k=noise_k)
            # Decision contract: an effect smaller than the guide-noise floor
            # cannot be reproduced by either method, so it is not scored.
            if abs(held["delta"]) <= floor:
                continue
            bands = [
                symmetric_bounds(r["m_a0"], r["m_0b"], r["m_ab"], r["mu"], rho, noise=floor)
                for _, r in cal.iterrows()
            ]
            robust = _robust_call(bands, rho)
            conventional = _conventional_call(cal["delta"].tolist())
            obs = classify(held["delta"], held["delta"])
            if obs == "UNRESOLVED":
                continue
            if robust != "UNRESOLVED":
                robust_calls += 1
                robust_hits += robust == obs
            if conventional != "UNRESOLVED":
                conv_calls += 1
                conv_hits += conventional == obs
        rows.append(
            {
                "gene_a": ga,
                "gene_b": gb,
                "rho": rho,
                "n_guide_pairs": k,
                "robust_calls": robust_calls,
                "robust_concordance": robust_hits / robust_calls if robust_calls else None,
                "conv_calls": conv_calls,
                "conv_concordance": conv_hits / conv_calls if conv_calls else None,
            }
        )
    return pd.DataFrame(rows)


def aggregate_holdout(out: pd.DataFrame) -> dict:
    """Coverage-matched comparison: robust vs conventional concordance."""
    if out.empty:
        return {"gene_pairs": 0}
    rob = out.dropna(subset=["robust_concordance"])
    conv = out.dropna(subset=["conv_concordance"])
    return {
        "gene_pairs": len(out),
        "gene_pairs_robust_callable": len(rob),
        "robust_concordance_mean": float(rob["robust_concordance"].mean()) if len(rob) else None,
        "conv_concordance_over_robust_pairs": float(conv.loc[rob.index, "conv_concordance"].mean()) if len(rob) else None,
        "conv_concordance_overall": float(conv["conv_concordance"].mean()) if len(conv) else None,
    }