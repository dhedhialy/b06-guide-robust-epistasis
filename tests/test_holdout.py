"""Guide-identity holdout (Sec 10)."""

import random

import pandas as pd

from b06.holdout import aggregate_holdout, evaluate_guide_holdout


def _synthetic(seed=0, rho=0.1):
    """Gene pairs with a real interaction + idiosyncratic guide-level noise."""
    rng = random.Random(seed)
    rows = []
    mu = 0.0
    for i in range(12):
        gamma = rng.choice([-0.8, 0.8])  # real interaction, strong vs envelope
        da, db = 0.5, 0.6
        for k in range(5):  # guide pairs per gene pair
            m_a0 = mu + da + rng.gauss(0, 0.05)
            m_0b = mu + db + rng.gauss(0, 0.05)
            m_ab = mu + gamma + rng.gauss(0, 0.05)
            rows.append({"gene_a": f"G{i}A", "gene_b": f"G{i}B",
                         "guide_a": f"g{i}a{k}", "guide_b": f"g{i}b{k}",
                         "m_a0": m_a0, "m_0b": m_0b, "m_ab": m_ab,
                         "mu": mu, "delta": m_ab - m_a0 - m_0b + mu})
    return pd.DataFrame(rows)


def test_holdout_runs_and_coverage_contracts():
    ct = _synthetic()
    out = evaluate_guide_holdout(ct, rho=0.10)
    assert len(out) == 12
    ok05 = evaluate_guide_holdout(ct, rho=0.05)
    ok20 = evaluate_guide_holdout(ct, rho=0.20)
    assert ok05["robust_calls"].sum() >= out["robust_calls"].sum() >= ok20["robust_calls"].sum()


def test_holdout_aggregate_shape():
    ct = _synthetic()
    agg = aggregate_holdout(evaluate_guide_holdout(ct, rho=0.10))
    assert agg["gene_pairs"] == 12
    assert 0 <= agg["robust_concordance_mean"] <= 1