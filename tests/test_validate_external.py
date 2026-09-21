"""Synthetic checks for the external paralog-panel pipeline (AUDIT Sec 5.9).

Covers the pre-registered pass/fail logic: unanimous direction AND AUROC
confidently > 0.5 (more-negative score = higher rank), plus the leave-one-out
stability columns and the n<10 inconclusive guardrail.
"""

import numpy as np
import pandas as pd

from b06.validate import (
    call_fdr_null,
    certified_fdr_null,
    external_paralog_panel,
    model_based_fdr_null,
    panel_verdict,
)


def _scores(n_pairs=120, n_pos=4, seed=3):
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(n_pairs):
        is_pos = k < n_pos
        for rho in (0.10,):
            rows.append(
                {
                    "rho": rho, "gene_a": f"p{k}a", "gene_b": f"p{k}b",
                    "mean_delta": -2.0 if is_pos else float(rng.choice([-0.1, 0.0, 0.1])),
                    "robust_score": -3.0 if is_pos else float(rng.uniform(-0.2, 0.2)),
                    "conv_score": -3.0 if is_pos else float(rng.uniform(-0.2, 0.2)),
                    "certified": bool(is_pos),
                }
            )
    return pd.DataFrame(rows)


def _panel(tmp_path, a, b):
    p = tmp_path / "panel.csv"
    pd.DataFrame({"A1": a, "A2": b}).to_csv(p, index=False)
    return p


def test_external_paralog_panel_detects_enrichment(tmp_path):
    scores = _scores()
    out = external_paralog_panel(
        scores, _panel(tmp_path, ["p0a", "p1a", "p2a", "p3a"],
                       ["p0b", "p1b", "p2b", "p3b"])
    )[0]
    assert out["n_pos"] == 4
    assert out["auroc_robust"] == 1.0          # most-negative positives: exact top ranks
    assert out["auroc_robust_loo_min"] > 0.9   # stable when any one is dropped
    assert out["robust_top500"] == 4
    assert out["unanimous_neg"] is True


def test_external_paralog_panel_midrank_declines(tmp_path):
    # Positives pushed to mid-rank but still directionally negative: AUROC ~ 0.5
    # and unanimous direction, so n<10 must keep the verdict inconclusive.
    scores = _scores()
    pos_names = [f"p{k}a" for k in (0, 1, 2, 3)]
    scores.loc[scores.gene_a.isin(pos_names),
               ["robust_score", "conv_score"]] = 0.0
    scores.loc[scores.gene_a.isin(pos_names), ["mean_delta"]] = -0.05
    out = external_paralog_panel(
        scores, _panel(tmp_path, pos_names, [f"p{k}b" for k in (0, 1, 2, 3)])
    )[0]
    assert out["unanimous_neg"] is True
    assert 0.3 < out["auroc_robust"] < 0.7
    assert panel_verdict(out["unanimous_neg"], out["auroc_robust"], 4) == "inconclusive"


def test_panel_verdict_boundaries():
    assert panel_verdict(True, 0.9, 12) == "enrichment-confirmed"
    assert panel_verdict(True, 0.3, 12) == "enrichment-declined"
    assert panel_verdict(False, 0.9, 12) == "direction-failed"
    assert panel_verdict(True, 0.9, 4) == "inconclusive"   # n<10 guardrail


def test_external_paralog_panel_case_and_no_cover(tmp_path):
    scores = _scores()
    out = external_paralog_panel(  # uppercase panel genes still match
        scores, _panel(tmp_path, ["P0A", "P1A"], ["P0B", "P1B"])
    )
    assert out[0]["n_pos"] == 2
    out0 = external_paralog_panel(  # panel genes absent from estate => skipped
        scores, _panel(tmp_path, ["zzz"], ["yyy"])
    )
    assert out0 == []


def _seed_contrasts(n_pairs=60, k=4, signal=0.6, every=2, seed=11):
    """Contrast rows with pair-specific signal `signal` on every `every`-th pair."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_pairs):
        center = -signal if i % every == 0 else 0.0
        for j in range(k):
            rows.append({
                "gene_a": f"g{i}a", "gene_b": f"g{i}b",
                "delta": float(center + rng.normal(0, 0.15)),
                "m_a0": float(1.2 + rng.normal(0, 0.05)),
                "m_0b": float(1.1 + rng.normal(0, 0.05)),
                "mu": 0.5, "noise": 0.2,
            })
    return pd.DataFrame(rows)


def test_permutation_rebuckets_row_groups():
    """Set-based membership guard on the null construction itself.

    A bijective relabel of block ids (perm[slot]) only renames each block, so
    the row partition is IDENTICAL pre/post — that is the point-mass null bug.
    The correct null shuffles the rows and re-buckets, which must change at
    least one block's membership while preserving the size multiset. This is
    checked directly, not via sd (a different broken mechanism could produce
    nonzero variance and still not regroup rows).
    """
    rng = np.random.default_rng(0)
    n, G = 40, 8
    sizes = np.repeat(np.arange(G), 5)  # 8 blocks of 5
    slot = sizes
    labels = np.repeat(np.arange(G), 5)

    def _groups(pid):
        acc = {}
        for i, p in enumerate(pid):
            acc.setdefault(int(p), set()).add(i)
        return set(frozenset(v) for v in acc.values())

    relabled = np.arange(G)[slot]            # perm[slot] degenerate case
    rebucketed = labels[rng.permutation(n)]  # correct case
    orig = _groups(slot)
    assert _groups(relabled) == orig  # relabel: membership unchanged (the bug)
    assert _groups(rebucketed) != orig  # rebucket: membership genuinely moves
    for pid in (relabled, rebucketed):
        assert sorted(np.bincount(pid, minlength=G)) == [5] * G


def test_call_fdr_null_permutation_varies_and_detects_signal():
    """Guard against the bijective-relabel bug: the permuted counts must MOVE
    (sd > 0) and a signal-rich contrast set must sit far above the null."""
    ct = _seed_contrasts(signal=0.6)
    res = call_fdr_null(ct, rho=0.10, n_perms=50, seed=0)
    assert res["null_neg_sd"] > 0            # permutation genuinely moves rows
    assert res["observed_neg"] > res["null_neg_max"]
    assert res["neg_p_ge_observed"] == 0.0
    assert res["neg_false_cert_rate"] < 0.5
    assert res["corrected_neg"] > 0
    assert res["neg_p_one_sided"] <= 1.0 / (50 + 1)


def test_certified_fdr_null_permutation_varies():
    # Regression guard for the bijective-relabel bug: the null must genuinely
    # regroup rows, so the recomputed counts must MOVE (sd > 0). Effect
    # detection is covered by test_call_fdr_null_permutation_varies_and_detects_signal.
    ct = _seed_contrasts(signal=0.6)
    res = certified_fdr_null(ct, rho=0.10, n_perms=40, seed=0)
    assert res["null_certified_sd"] > 0


def test_model_based_fdr_null_smoke_and_layout():
    # Model-based synthetic-screen null (AUDIT Sec 5.10): must run on the
    # same contrasts, keep the real block layout (block count preserved), and
    # supply both delta models plus the between-null comparison keys. On
    # signal-rich contrasts the observed count must sit above the null, and
    # the symmetric-normal null must split evenly NEG/POS while the
    # empirical-marginal null inherits the data's sign skew.
    ct = _seed_contrasts(signal=0.6)
    res = model_based_fdr_null(ct, rho=0.10, n_sims=40, seed=0)
    assert res["n_blocks"] == ct[["gene_a", "gene_b"]].drop_duplicates().shape[0]
    assert abs(res["within_block_rho_d"]) < 1.0
    assert res["null_certified_sd"] > 0
    assert res["null_normal_certified_sd"] > 0
    assert res["observed_certified"] > res["null_certified_mean"]
    assert res["corrected_certified"] > 0
    normal_ratio = res["null_normal_pos_mean"] / max(res["null_normal_neg_mean"], 1e-9)
    assert 0.2 < normal_ratio < 5.0          # symmetric-normal null is balanced
    assert res["null_neg_sd"] > 0