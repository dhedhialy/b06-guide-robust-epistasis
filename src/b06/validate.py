"""Audit-side validation of the robust certificate as a *ranking* score.

Two dependency-free checks on the already-computed guide-pair contrasts:

1. LOOCV ranking AUC — per held-out guide pair, does the calibration-derived
   robust score order observed negative-vs-positive guide interactions as
   well as (better than?) the raw mean delta? Matches the field's AUROC
   idiom for SL detection.
2. Validated-paralog rank test — where do the paper's five validated paralog
   SL pairs rank among all 17k gene pairs under each score?

Reads `bounds/guide_pair_contrasts.csv`, writes audit tables + a JSON
summary. No external gold-standard files required.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .bounds import classify
from .contrasts import NOISE_K, RHO_GRID, gene_scores

# The paper's five experimentally validated paralog SL pairs (their Fig/Suppl
# validation list), canonicalizing gene order like the pipeline.
VALIDATED_PARALOG_SL = {
    ("ARID1A", "ARID1B"),
    ("HDAC1", "HDAC2"),
    ("MAPK1", "MAPK3"),
    ("ASF1A", "ASF1B"),
    ("CNOT7", "CNOT8"),
}


def _auc(scores, labels):
    """Mann-Whitney U AUC: higher score => label 1 (negative-GI)."""
    pos = np.where(labels == 1)[0]
    neg = np.where(labels == 0)[0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    ordered = np.argsort(scores)
    ranks = np.empty(len(scores), dtype=float)
    ranks[ordered] = np.arange(len(scores)) + 1.0
    r_pos = ranks[pos].sum()
    n_pos, n_neg = len(pos), len(neg)
    return float((r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def holdout_ranking_auc(contrasts: pd.DataFrame, rho: float, min_calibration: int = 2):
    """Per-gene-pair LOOCV: predict each held-out guide pair's interaction
    sign from the calibration-derived score; labels are above-floor effects
    only (|delta| > calibration MAD), the field's decision contract."""
    scores_sl, labels_sl = [], []
    scores_conv = [] 
    for (gene_a, gene_b), idx in contrasts.groupby(["gene_a", "gene_b"]).groups.items():
        pos = contrasts.index.get_indexer(idx)
        if len(pos) < min_calibration + 1:
            continue
        for jj, j in enumerate(pos):
            mask = np.ones(len(pos), dtype=bool)
            mask[jj] = False
            cal = pos[mask]
            cal_d = contrasts["delta"].to_numpy()[cal]
            held_d = contrasts["delta"].to_numpy()[j]
            floor = NOISE_K * 1.4826 * float(np.median(np.abs(cal_d - np.median(cal_d))))
            if held_d > floor:
                label = 0  # positive interaction (not SL)
            elif held_d < -floor:
                label = 1  # negative interaction / SL-like
            else:
                continue  # sub-floor held-out effect: not scored
            env = rho * (
                np.abs(contrasts["m_a0"].to_numpy()[cal] - contrasts["mu"].to_numpy()[cal]).mean()
                + np.abs(contrasts["m_0b"].to_numpy()[cal] - contrasts["mu"].to_numpy()[cal]).mean()
            ) + floor
            mean_cal = cal_d.mean()
            # signed: more negative = stronger SL evidence. _auc expects
            # higher=SL, so store -score.
            scores_sl.append(-(mean_cal / env))
            scores_conv.append(-mean_cal)
            labels_sl.append(label)
    labels = np.array(labels_sl, dtype=int)
    return {
        "rho": rho,
        "instances": len(labels),
        "n_negative_gi": int((labels == 1).sum()),
        "n_positive_gi": int((labels == 0).sum()),
        "robust_auc": _auc(np.array(scores_sl), labels),
        "conv_auc": _auc(np.array(scores_conv), labels),
    }


def paralog_sl_ranks(contrasts: pd.DataFrame, rho: float) -> list[dict]:
    """Rank position of each validated SL pair under robust vs conv score."""
    s = gene_scores(contrasts, rho).copy()
    s["rank_robust"] = s["robust_score"].rank(ascending=True)
    s["rank_conv"] = s["conv_score"].rank(ascending=True)
    rows = []
    for (ga, gb) in VALIDATED_PARALOG_SL:
        hit = s[(s.gene_a == ga) & (s.gene_b == gb)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        rows.append(
            {
                "gene_a": ga,
                "gene_b": gb,
                "rho": rho,
                "n_guide_pairs": int(r.n_guide_pairs),
                "mean_delta": round(float(r.mean_delta), 3),
                "robust_score": round(float(r.robust_score), 3),
                "robust_certified": bool(r.certified),
                "rank_robust": int(r.rank_robust),
                "rank_conv": int(r.rank_conv),
                "n_gene_pairs": len(s),
            }
        )
    return rows


def _gene_pair_conv_precision(contrasts: pd.DataFrame) -> pd.DataFrame:
    """Per-gene-pair LOOCV mean-sign precision (conventional calibration call
    vs held-out guide-pair sign, above-floor contract). This is the metric
    the selectivity ablation compares across matched-size subsets."""
    rows = []
    d = contrasts["delta"].to_numpy()
    for (ga, gb), idx in contrasts.groupby(["gene_a", "gene_b"]).groups.items():
        pos = contrasts.index.get_indexer(idx)
        if len(pos) < 3:
            continue
        hits = total = 0
        for jj, j in enumerate(pos):
            mask = np.ones(len(pos), dtype=bool)
            mask[jj] = False
            cal = d[pos[mask]]
            held = d[j]
            floor = 1.4826 * float(np.median(np.abs(cal - np.median(cal))))
            if abs(held) <= floor:
                continue
            cal_sign = classify(cal.mean(), cal.mean())
            if cal_sign == "UNRESOLVED":
                continue
            total += 1
            hits += (cal_sign == classify(held, held))
        rows.append(
            {
                "gene_a": ga,
                "gene_b": gb,
                "precision": hits / total if total else None,
                "n_instances": total,
            }
        )
    return pd.DataFrame(rows)


def selectivity_ablation(contrasts: pd.DataFrame, rho: float) -> dict:
    """Does the `certified` set buy reproducibility beyond effect magnitude?

    Compares mean-sign precision (conventional calibration call, LOOCV) on
    the certified set vs the same-size top-|mean_delta| set. If certified
    does not beat magnitude-matched selection, the certificate's selectivity
    claim reduces to "big effects reproduce better" and must shrink.
    """
    s = gene_scores(contrasts, rho).sort_values("mean_delta", key=lambda x: x.abs(), ascending=False)
    cert = set(zip(s.loc[s.certified, "gene_a"], s.loc[s.certified, "gene_b"]))
    prec = _gene_pair_conv_precision(contrasts).set_index(["gene_a", "gene_b"])
    k = len(cert)
    mag = set(zip(s.head(k)["gene_a"], s.head(k)["gene_b"]))
    cert_rows = [(ga, gb) for ga, gb in cert if (ga, gb) in prec.index]
    mag_rows = [(ga, gb) for ga, gb in mag if (ga, gb) in prec.index]
    p_c = prec.loc[cert_rows, "precision"].dropna()
    p_m = prec.loc[mag_rows, "precision"].dropna()
    return {
        "rho": rho,
        "n_certified": len(cert),
        "n_certified_scored": int(len(p_c)),
        "n_magnitude_matched": len(mag),
        "n_magnitude_scored": int(len(p_m)),
        "certified_set_precision": float(p_c.mean()) if len(p_c) else None,
        "magnitude_matched_precision": float(p_m.mean()) if len(p_m) else None,
        "overlap": len(set(cert_rows) & set(mag_rows)),
    }


def noise_k_sweep(contrasts: pd.DataFrame, rho: float, noise_ks=(0.5, 1.0, 1.5, 2.0)) -> list[dict]:
    """Holdout concordance across the noise-floor multiplier (sensitivity)."""
    from .holdout import evaluate_guide_holdout

    rows = []
    for k in noise_ks:
        h = evaluate_guide_holdout(contrasts, rho, noise_k=k)
        rob = h.dropna(subset=["robust_concordance"])
        cvg = h.dropna(subset=["conv_concordance"])
        rows.append(
            {
                "rho": rho,
                "noise_k": k,
                "robust_callable": int(len(rob)),
                "robust_concordance": float(rob["robust_concordance"].mean()) if len(rob) else None,
                "conv_concordance": float(cvg.loc[rob.index, "conv_concordance"].mean()) if len(rob) else None,
            }
        )
    return rows


def external_paralog_panel(scores_by_rho: pd.DataFrame, panel_path: Path) -> list[dict]:
    """Pre-registered external reference test (AUDIT Sec 5.9).

    De Kegel et al. 2021 (Cell Systems, DOI 10.1016/j.cels.2021.08.006)
    validated_SLs.txt is an independent computationally-derived paralog-SL
    panel with no shared lineage to this screen. Only pairs whose genes are
    both in the matched-guide estate are evaluable (power caveat stated in
    AUDIT). Metrics fixed before reading outcomes: rank AUROC vs the rest of
    the universe, top-K recall (K=50/100/500), unanimous sign of mean_delta,
    certified fraction of the positives.
    """
    panel = pd.read_csv(panel_path)
    pos = {frozenset((r.A1.upper(), r.A2.upper())) for _, r in panel.iterrows()}
    rows = []
    for rho, d in scores_by_rho.groupby("rho"):
        d = d.copy()
        d["pair"] = [frozenset((a.upper(), b.upper())) for a, b in zip(d.gene_a, d.gene_b)]
        m = d.pair.isin(pos)
        n_pos, n_neg = int(m.sum()), int((~m).sum())
        if n_pos == 0:
            continue
        rec = {"rho": float(rho), "n_pos": n_pos, "n_neg": n_neg,
               "n_panel_covered": n_pos}
        for tag in ("robust", "conv"):
            col = "robust_score" if tag == "robust" else "conv_score"
            ordered = np.argsort(-d[col].to_numpy())  # descending: more-negative = higher = SL-like
            ranks = np.empty(len(d), dtype=float)
            ranks[ordered] = np.arange(len(d)) + 1.0
            r_pos = ranks[np.where(m.to_numpy())[0]].sum()
            rec[f"auroc_{tag}"] = round((r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg), 4)
            pos_idx = np.where(m.to_numpy())[0]
            if n_pos > 1:  # leave-one-out stability: range of AUROC dropping each positive
                loo = []
                for j in pos_idx:
                    keep = np.ones(len(d), dtype=bool)
                    keep[j] = False
                    pj = ranks[keep][m[keep].to_numpy()].sum()
                    loo.append((pj - (n_pos - 1) * n_pos / 2) / ((n_pos - 1) * n_neg))
                rec[f"auroc_{tag}_loo_min"] = round(float(min(loo)), 4)
                rec[f"auroc_{tag}_loo_max"] = round(float(max(loo)), 4)
            top = d.sort_values(col, ascending=False).head(max(50, 500))  # single sort, reuse
            for K in (50, 100, 500):
                rec[f"{tag}_top{K}"] = int(top.head(K)["pair"].isin(pos).sum())
        rec["neg_dir"] = int((d.loc[m, "mean_delta"] < 0).sum())
        rec["unanimous_neg"] = bool((d.loc[m, "mean_delta"] < 0).all())
        rec["cert_frac"] = round(float(d.loc[m, "certified"].mean()), 3)
        rows.append(rec)
    return rows


def _permuted_slot_labels(slot):
    """Correct guide-pair-label permutation null.

    A bijective relabel of the per-row factor codes only renames the gene-pair
    blocks --- every guide pair stays with its original partners, so *nothing*
    is permuted and the recomputed counts are constant (sd=0). The null must
    reshuffle the rows and re-bucket them into blocks of the same size
    multiset, so the guide-pair contrasts genuinely change group membership.
    """
    G = int(slot.max()) + 1
    sizes = np.bincount(slot, minlength=G)
    labels = np.repeat(np.arange(G), sizes)
    rng = np.random.default_rng()
    return labels[rng.permutation(len(slot))]


def _certified_counts(delta, abs_a, abs_b, noise, pid, rho: float):
    """Gene-level certified counts from per-row arrays + block assignment.

    `certified = |mean_delta| > rho*(mean|d_a| + mean|d_b|) + median(noise)`
    per block, plus the signed split. Shared by the permutation null and the
    model-based synthetic-screen null (gap 3b).
    """
    G = int(pid.max()) + 1
    cnt = np.bincount(pid, minlength=G)
    mean_d = np.bincount(pid, weights=delta, minlength=G) / cnt
    mean_a = np.bincount(pid, weights=abs_a, minlength=G) / cnt
    mean_b = np.bincount(pid, weights=abs_b, minlength=G) / cnt
    order = np.argsort(noise, kind="stable")
    sids = pid[order]
    sv = noise[order]
    starts = np.concatenate([[0], np.flatnonzero(sids[1:] != sids[:-1]) + 1])
    ends = np.concatenate([starts[1:], [len(sids)]])
    sizes_g = ends - starts
    med = (sv[starts + (sizes_g - 1) // 2] + sv[starts + sizes_g // 2]) / 2
    med_g = np.empty(G)
    med_g[sids[starts]] = med
    env = rho * (mean_a + mean_b) + med_g
    cert = np.abs(mean_d) > env
    return (int(cert.sum()), int(((mean_d < -env) & cert).sum()),
            int(((mean_d > env) & cert).sum()))


def model_based_fdr_null(contrasts: pd.DataFrame, rho: float = 0.10,
                         n_sims: int = 100, seed: int = 0) -> dict:
    """Model-based synthetic-screen null (AUDIT gap 3b).

    The permutation null (certified_fdr_null) re-buckets real rows across
    unrelated genes, so a synthetic block is a mix of main effects from
    genes that never pair in reality: the within-block `d_a`/`d_b`
    correlation of *real* gene pairs is pooled away.

    This H0 keeps the real block layout — every block's guide-pair
    multiplicity AND its real per-row covariates (|d_a|, |d_b|, noise), so
    the observed within-block main-effect correlation is preserved by
    construction — and replaces only the signal: per block, the guide-pair
    deltas are redrawn against a fitted null. Two delta models are fit:

      * "empirical": deltas sampled i.i.d. from the observed marginal
        (keeps the screen's negative-skewed delta distribution);
      * "normal": deltas drawn i.i.d. from N(0, var(delta_block)) per block
        (symmetric; shows how much the perm-null center depends on the
        delta skew inherited from the real deltas it shuffles).

    Disentangling: the permutation null mixes rows across genes AND reuses
    real (skewed) deltas; the model null keeps the correlated covariate
    blocks fixed. Agreement between the perm null and the "empirical" model
    null means the pooled-away within-block structure does not move the FCR
    reading; divergence between "empirical" and "normal" isolates the skew
    assumption.
    """
    mu = contrasts["mu"].to_numpy(float)
    abs_a = (contrasts["m_a0"] - mu).abs().to_numpy(float)
    abs_b = (contrasts["m_0b"] - mu).abs().to_numpy(float)
    noise = contrasts["noise"].to_numpy(float)
    real_delta = contrasts["delta"].to_numpy(float)

    cf = contrasts[["gene_a", "gene_b"]].apply(tuple, axis=1)
    slot, _ = pd.factorize(cf)
    gb = contrasts.groupby(["gene_a", "gene_b"])

    # observed within-block main-effect correlation (the pooled-away structure)
    da = contrasts["m_a0"].to_numpy(float) - mu
    db = contrasts["m_0b"].to_numpy(float) - mu
    rho_d = float(np.corrcoef(da, db)[0, 1])

    sizes_arr = gb["delta"].size().to_numpy()
    sigma2 = gb["delta"].var(ddof=1).fillna(0.0).to_numpy().clip(min=1e-9)
    sig = np.sqrt(np.repeat(sigma2, sizes_arr))
    rng = np.random.default_rng(seed)

    def _sims(gen):
        counts, neg, pos = [], [], []
        for _ in range(n_sims):
            delta_s = gen()
            c, ng, pg = _certified_counts(delta_s, abs_a, abs_b, noise,
                                          slot, rho)
            counts.append(c); neg.append(ng); pos.append(pg)
        return np.array(counts), np.array(neg), np.array(pos)

    c_f, neg_f, pos_f = _sims(lambda: rng.choice(real_delta, size=len(real_delta),
                                                 replace=True))
    c_n, neg_n, pos_n = _sims(lambda: rng.normal(scale=sig))

    obs_cert, obs_neg, obs_pos = _certified_counts(real_delta, abs_a, abs_b,
                                                   noise, slot, rho)
    return {
        "rho": rho, "n_sims": n_sims, "seed": seed,
        "observed_certified": obs_cert, "observed_neg": obs_neg,
        "observed_pos": obs_pos,
        "within_block_rho_d": rho_d,
        "n_blocks": int(len(sizes_arr)),
        # empirical-marginal delta model (real covariate blocks preserved)
        "null_certified_mean": float(c_f.mean()),
        "null_certified_sd": float(c_f.std()),
        "null_neg_mean": float(neg_f.mean()), "null_neg_sd": float(neg_f.std()),
        "null_pos_mean": float(pos_f.mean()), "null_pos_sd": float(pos_f.std()),
        # symmetric-normal delta model (skew-removed, same covariate blocks)
        "null_normal_certified_mean": float(c_n.mean()),
        "null_normal_certified_sd": float(c_n.std()),
        "null_normal_neg_mean": float(neg_n.mean()),
        "null_normal_pos_mean": float(pos_n.mean()),
        # FCR-scaled call counts (empirical-marginal null)
        "corrected_certified": float(obs_cert - c_f.mean()),
        "corrected_neg": float(obs_neg - neg_f.mean()),
        "corrected_pos": float(obs_pos - pos_f.mean()) if obs_pos else None,
    }


def certified_fdr_null(contrasts: pd.DataFrame, rho: float = 0.10,
                       n_perms: int = 300, seed: int = 0) -> dict:
    """Permutation null for the certified-pair count (AUDIT Sec 5.1/FDR).

    Null: reassign each guide-pair contrast to a gene-pair slot, preserving
    every slot's guide-pair multiplicity and each row's covariates (delta,
    |d_a|, |d_b|, noise floor). Recomputes only the gene-level aggregates
    behind `certified` (grouped means + grouped median of the precomputed
    noise column) plus the signed split.
    """
    delta = contrasts["delta"].to_numpy(float)
    abs_a = (contrasts["m_a0"] - contrasts["mu"]).abs().to_numpy(float)
    abs_b = (contrasts["m_0b"] - contrasts["mu"]).abs().to_numpy(float)
    noise = contrasts["noise"].to_numpy(float)
    cf = contrasts[["gene_a", "gene_b"]].apply(tuple, axis=1)
    slot, _ = pd.factorize(cf)
    G = int(slot.max()) + 1
    sizes = np.bincount(slot, minlength=G)
    labels = np.repeat(np.arange(G), sizes)
    rng = np.random.default_rng(seed)

    def _perm_certified(pid):
        return _certified_counts(delta, abs_a, abs_b, noise, pid, rho)

    obs = _perm_certified(slot)
    counts, neg, pos = [], [], []
    for _ in range(n_perms):
        c, ng, pg = _perm_certified(labels[rng.permutation(len(slot))])
        counts.append(c); neg.append(ng); pos.append(pg)
    counts = np.array(counts); neg = np.array(neg); pos = np.array(pos)
    return {
        "rho": rho, "n_perms": n_perms, "seed": seed,
        "observed_certified": obs[0], "observed_neg": obs[1], "observed_pos": obs[2],
        "null_certified_mean": float(counts.mean()),
        "null_certified_sd": float(counts.std()),
        "null_neg_mean": float(neg.mean()), "null_neg_sd": float(neg.std()),
        "null_pos_mean": float(pos.mean()), "null_pos_sd": float(pos.std()),
        "neg_false_cert_rate": float(neg.mean() / obs[1]),
        "pos_false_cert_rate": float(pos.mean() / obs[2]) if obs[2] else None,
        "neg_p_ge_observed": float((neg >= obs[1]).mean()),
        "pos_p_ge_observed": float((pos >= obs[2]).mean()),
        # (k+1)/(n+1) one-sided Monte-Carlo p (Phipson-Smyth), resolution 1/(n+1)
        "neg_p_one_sided": float(((neg >= obs[1]).sum() + 1) / (n_perms + 1)),
        "pos_p_one_sided": float(((pos >= obs[2]).sum() + 1) / (n_perms + 1)),
        "certified_perms_ge_observed": int((counts >= obs[0]).sum()),
        "neg_perms_ge_observed": int((neg >= obs[1]).sum()),
        "pos_perms_ge_observed": int((pos >= obs[2]).sum()),
        # FCR-scaled call count: expected TRUE calls among the observed set
        "corrected_neg": float(obs[1] - neg.mean()),
        "corrected_pos": float(obs[2] - pos.mean()),
        "corrected_certified": float(obs[0] - counts.mean()),
    }


def call_fdr_null(contrasts: pd.DataFrame, rho: float = 0.10,
                  n_perms: int = 1000, seed: int = 0) -> dict:
    """Permutation null for the strict Sec 5.1 POSITIVE/NEGATIVE call counts.

    The call rule is: NEGATIVE iff every guide-pair interval resolves negative
    (max upper < 0); POSITIVE iff every guide-pair interval resolves positive
    (min lower > 0). Max/min are not additive, so the recomputed counts move
    under the row-shuffle null (no bijective-relabel degeneracy).

    Modeling choice (named, not implicit): the null preserves each block's
    size multiset and each contrast row's covariates (delta, |d_a|, |d_b|,
    noise), but shuffling re-buckets rows across unrelated genes, so a
    synthetic block pools main effects from different genes. Real gene pairs
    have within-block `d_a`/`d_b` that are correlated (same gene, shared guide
    context); permuted blocks do not. The null therefore tests "counts as if
    rows carried no true pair correspondence" — a valid no-association null —
    while the main-effect correlation structure of *real* pairs is only
    approximated. Direction of the resulting FCR bias is unknown; the check
    is an order-of-magnitude one (obs 4-180x null center) that a more
    faithful model-based null (Sec 7 gap) may sharpen.
    """
    mu = contrasts["mu"].to_numpy(float)
    delta = contrasts["delta"].to_numpy(float)
    noise = contrasts["noise"].to_numpy(float)
    abs_a = (contrasts["m_a0"] - mu).abs().to_numpy(float)
    abs_b = (contrasts["m_0b"] - mu).abs().to_numpy(float)
    half = rho * (abs_a + abs_b) + noise
    lo, hi = delta - half, delta + half
    cf = contrasts[["gene_a", "gene_b"]].apply(tuple, axis=1)
    slot, _ = pd.factorize(cf)
    G = int(slot.max()) + 1
    sizes = np.bincount(slot, minlength=G)
    labels = np.repeat(np.arange(G), sizes)
    rng = np.random.default_rng(seed)

    def _calls(pid):
        hi_max = np.full(G, -np.inf); np.maximum.at(hi_max, pid, hi)
        lo_min = np.full(G, np.inf); np.minimum.at(lo_min, pid, lo)
        return int((hi_max < 0).sum()), int((lo_min > 0).sum())

    obs_neg, obs_pos = _calls(slot)
    nn, np2 = [], []
    for _ in range(n_perms):
        a, b = _calls(labels[rng.permutation(len(slot))])
        nn.append(a); np2.append(b)
    nn = np.array(nn); np2 = np.array(np2)
    return {
        "rho": rho, "n_perms": n_perms, "seed": seed,
        "observed_neg": obs_neg, "observed_pos": obs_pos,
        "null_neg_mean": float(nn.mean()), "null_neg_sd": float(nn.std()),
        "null_neg_max": int(nn.max()), "null_pos_mean": float(np2.mean()),
        "null_pos_sd": float(np2.std()), "null_pos_max": int(np2.max()),
        "neg_false_cert_rate": float(nn.mean() / obs_neg),
        "pos_false_cert_rate": float(np2.mean() / obs_pos) if obs_pos else None,
        "neg_p_ge_observed": float((nn >= obs_neg).mean()),
        "pos_p_ge_observed": float((np2 >= obs_pos).mean()),
        # (k+1)/(n+1) one-sided Monte-Carlo p (Phipson-Smyth), resolution 1/(n+1)
        "neg_p_one_sided": float(((nn >= obs_neg).sum() + 1) / (n_perms + 1)),
        "pos_p_one_sided": float(((np2 >= obs_pos).sum() + 1) / (n_perms + 1)),
        "neg_perms_ge_observed": int((nn >= obs_neg).sum()),
        "pos_perms_ge_observed": int((np2 >= obs_pos).sum()),
        # FCR-scaled call count: expected TRUE calls among the observed set
        "corrected_neg": float(obs_neg - nn.mean()),
        "corrected_pos": float(obs_pos - np2.mean()),
    }


def panel_verdict(unanimous_neg: bool, auroc_robust: float, n_pos: int) -> str:
    """Pre-registered external-panel decision rule (AUDIT Sec 5.9).

    Enrichment is only claimed when direction is unanimous AND the AUROC is
    confidently > 0.5 (AUROC as computed by external_paralog_panel: positive
    class = validated-SL pair, more-negative score = higher rank). At n<10 a
    one-pair swing is indistinguishable from a systematic failure, so the
    verdict defaults to inconclusive.
    """
    if not unanimous_neg:
        return "direction-failed"
    if n_pos < 10:
        return "inconclusive"
    return "enrichment-confirmed" if auroc_robust > 0.5 else "enrichment-declined"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contrasts", type=Path, required=True,
                    help="guide_pair_contrasts.csv from a scout run")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--external-panel", type=Path,
                    default=Path("raw/de_kegel_2021/validated_SLs.txt"),
                    help="De Kegel 2021 validated_SLs.txt (A1,A2 symbols)")
    ap.add_argument("--n-perms", type=int, default=1000,
                    help="permutations for the Sec-5.1 and certified-count nulls")
    ap.add_argument("--model-null-sims", type=int, default=100,
                    help="synthetic screens for the model-based null (gap 3b)")
    args = ap.parse_args()

    ct = pd.read_csv(args.contrasts)
    (args.out / "audit").mkdir(parents=True, exist_ok=True)

    auc_rows, sl_rows = [], []
    score_frames = []
    for rho in RHO_GRID:
        auc_rows.append(holdout_ranking_auc(ct, rho))
        sl_rows += paralog_sl_ranks(ct, rho)
        s = gene_scores(ct, rho)
        s["rho"] = rho
        score_frames.append(s)

    auc = pd.DataFrame(auc_rows)
    sl = pd.DataFrame(sl_rows)
    scores = pd.concat(score_frames, ignore_index=True)

    ablation = pd.DataFrame([selectivity_ablation(ct, rho) for rho in RHO_GRID])
    sweep = pd.DataFrame(
        [r for rho in (RHO_GRID[1],) for r in noise_k_sweep(ct, rho)]  # rho=0.10
    )
    ext = pd.DataFrame(external_paralog_panel(scores, args.external_panel))
    fdr_cert = certified_fdr_null(ct, rho=RHO_GRID[1], n_perms=args.n_perms, seed=0)
    fdr_model = model_based_fdr_null(ct, rho=RHO_GRID[1],
                                     n_sims=args.model_null_sims, seed=0)
    fdr_model["perm_comparison"] = {
        "perm_null_certified_mean": fdr_cert["null_certified_mean"],
        "perm_null_certified_sd": fdr_cert["null_certified_sd"],
        "model_null_minus_perm": round(
            fdr_model["null_certified_mean"] - fdr_cert["null_certified_mean"], 2),
        "model_over_perm_ratio": round(
            fdr_model["null_certified_mean"] / fdr_cert["null_certified_mean"], 3),
    }
    fdr_calls = pd.DataFrame(
        [call_fdr_null(ct, rho=rho, n_perms=args.n_perms, seed=0) for rho in RHO_GRID]
    )

    auc.to_csv(args.out / "audit" / "rank_auc.csv", index=False)
    if not sl.empty:
        sl.to_csv(args.out / "audit" / "validated_paralog_ranks.csv", index=False)
    scores.to_csv(args.out / "audit" / "gene_scores.csv", index=False)
    ablation.to_csv(args.out / "audit" / "selectivity_ablation.csv", index=False)
    sweep.to_csv(args.out / "audit" / "noise_k_sweep.csv", index=False)
    if not ext.empty:
        ext.to_csv(args.out / "audit" / "dekegel_paralog_sl.csv", index=False)
    pd.DataFrame([fdr_cert]).to_json(args.out / "audit" / "certified_fdr_null.json", orient="records")
    fdr_calls.to_csv(args.out / "audit" / "call_fdr_null.csv", index=False)
    pd.Series(fdr_model).to_json(args.out / "audit" / "model_based_fdr_null.json")
    fdr_model_flat = {"model_based_fdr_null_" + k: (v if not isinstance(v, dict) else json.dumps(v))
                      for k, v in fdr_model.items()}

    summary = {
        "rank_auc_by_rho": auc_rows,
        "validated_paralog_sl": {
            "set": sorted(["|".join(p) for p in VALIDATED_PARALOG_SL]),
            "found": len(sl),
            "ranks_robust": [int(r["rank_robust"]) for r in sl_rows],
            "ranks_conv": [int(r["rank_conv"]) for r in sl_rows],
            "n_gene_pairs": len(s) // len(RHO_GRID),
        },
        "selectivity_ablation": ablation.to_dict("records"),
        "noise_k_sweep_rho0_10": sweep.to_dict("records"),
        "de_kegel_external_panel": ext.to_dict("records"),
        "certified_fdr_null_rho0_10": fdr_cert,
        "model_based_fdr_null_rho0_10": {
            k: v for k, v in fdr_model.items()
        },
        "call_fdr_null_by_rho": fdr_calls.to_dict("records"),
    }
    (args.out / "audit" / "validation_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(auc.to_string(index=False))
    print()
    print("--- selectivity ablation (certified vs magnitude-matched) ---")
    print(ablation.to_string(index=False))
    print()
    print("--- certified-count permutation null @ rho=0.10 ---")
    print(fdr_cert)
    print()
    print("--- model-based synthetic-screen null @ rho=0.10 (gap 3b) ---")
    print(fdr_model)
    print()
    print("--- strict Sec 5.1 call-count permutation null (by rho) ---")
    show = fdr_calls[["rho", "n_perms", "observed_neg", "corrected_neg",
                      "null_neg_mean", "null_neg_sd", "neg_perms_ge_observed",
                      "neg_p_one_sided", "neg_false_cert_rate",
                      "observed_pos", "corrected_pos", "null_pos_mean",
                      "null_pos_sd", "pos_perms_ge_observed",
                      "pos_p_one_sided", "pos_false_cert_rate"]]
    print(show.to_string(index=False))
    print()
    print("--- noise_k sweep @ rho=0.10 ---")
    print(sweep.to_string(index=False))
    print()
    if not sl.empty:
        print(sl[["gene_a", "gene_b", "robust_score", "robust_certified", "rank_robust", "rank_conv"]].to_string(index=False))


if __name__ == "__main__":
    main()