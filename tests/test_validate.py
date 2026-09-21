"""Synthetic checks for the validate module (AUC + LOOCV ranking)."""

import numpy as np
import pandas as pd

from b06.validate import _auc, holdout_ranking_auc


def test_auc_perfect_reverse():
    # _auc contract: HIGHER score => label 1 (negative-GI)
    scores = np.array([4.0, 3.0, 2.0, 1.0])  # label-1 items score highest
    labels = np.array([1, 1, 0, 0])
    assert _auc(scores, labels) == 1.0
    assert _auc(-scores, labels) == 0.0


def test_holdout_ranking_auc_synthetic():
    """Negative-GI gene pairs must outrank positive-GI pairs under both
    scores; calibration-only floors keep the LOOCV honest."""
    rng = np.random.default_rng(7)
    rows = []
    for k in range(120):
        kind = rng.integers(0, 3)  # 0 neutral, 1 negative-GI, 2 positive-GI
        center = {0: 0.0, 1: -0.8, 2: 0.8}[kind]
        for j in range(3):
            rows.append(
                {
                    "gene_a": f"g{k}a", "gene_b": f"g{k}b",
                    "delta": float(center + rng.normal(0.0, 0.15)),
                    "m_a0": 1.2, "m_0b": 1.1, "m_ab": 2.0, "mu": 0.5, "noise": 0.2,
                }
            )
    ct = pd.DataFrame(rows)
    res = holdout_ranking_auc(ct, rho=0.10)
    assert res["instances"] > 100
    assert res["robust_auc"] > 0.9
    assert res["conv_auc"] > 0.9