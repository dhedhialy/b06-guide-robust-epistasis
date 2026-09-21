"""SOTA benchmark (AUDIT Sec 6): GEMINI on this exact screen vs robust/conv.

Convention (matches external_paralog_panel in b06/validate.py): rank AUROC
under "higher score = more SL-like" (Mann-Whitney U: larger rank for the
labels). Score axes:
  * robust_score, conv_score  -> negative = SL, so the SL-strength axis is
    -score (raw score passes through argsort(-col) as in validate.py)
  * gemini strong             -> positive = SL natively (CNOT7-CNOT8 strongest
    validated paralog has strong = +1.19, rank #1 of 18,920; the inverted
    reading puts every validated paralog at the bottom and is rejected)

Self-check: auroc_robust / auroc_conv must equal the stored validation values
(~0.658 / ~0.775 at rho 0.10).

Outputs:
  * results_colo1/audit/gemini_benchmark.csv   per-rho benchmark table
  * results_colo1/audit/gemini_benchmark.json  machine-readable copy
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from b06.validate import VALIDATED_PARALOG_SL


def _rank_auc(score: np.ndarray, pos: np.ndarray) -> float | None:
    """Standard rank AUC: higher score -> label 1 (as in validate._auc)."""
    pos = np.array(pos, dtype=bool)
    idx = np.where(pos)[0]
    n_pos, n_neg = int(pos.sum()), int((~pos).sum())
    if n_pos == 0 or n_neg == 0:
        return None
    ordered = np.argsort(score)
    ranks = np.empty(len(score))
    ranks[ordered] = np.arange(len(score)) + 1.0
    r_pos = ranks[idx].sum()
    return float((r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gemini", type=Path, default="benchmarks/gemini/gemini_scores.csv")
    ap.add_argument("--scores", type=Path, default="results_colo1/audit/gene_scores.csv")
    ap.add_argument("--panel", type=Path, default="raw/de_kegel_2021/validated_SLs.txt")
    ap.add_argument("--out", type=Path, default="results_colo1/audit")
    args = ap.parse_args()

    g = pd.read_csv(args.gemini)
    g["pair"] = [frozenset(p.split(";")) for p in g["gene_pair"]]

    s = pd.read_csv(args.scores)
    s["pair"] = [frozenset((a.upper(), b.upper())) for a, b in zip(s.gene_a, s.gene_b)]

    panel = pd.read_csv(args.panel)
    pos = {frozenset((r.A1.upper(), r.A2.upper())) for _, r in panel.iterrows()}

    rows = []
    for rho in sorted(s.rho.unique()):
        d = s[s.rho == rho].drop_duplicates("pair")
        joined = d.merge(g, on="pair", how="inner")
        # gemini strong can be NaN for pairs with no passing guide pairs; keep
        # only rows with a gemini score (common universe metric)
        joined = joined[joined["gemini_strong_median"].notna()].copy()
        m = joined["pair"].isin(pos).to_numpy()
        row = {"rho": float(rho), "n_shared": len(joined),
               "n_panel_covered": int(m.sum())}

        # SL-strength axes: robust/conv negative=SL -> pass -score; gemini
        # positive=SL -> pass score as-is. Higher = more SL-like either way.
        axes = {"robust": -joined["robust_score"].to_numpy(float),
                "conv": -joined["conv_score"].to_numpy(float),
                "gemini": joined["gemini_strong_median"].to_numpy(float)}
        for tag, sl in axes.items():
            row[f"auroc_{tag}"] = round(x, 4) if (x := _rank_auc(sl, m)) else None
            if tag == "gemini" and m.sum() > 1:
                loo = []
                for k in np.where(m)[0]:
                    keep = np.ones(len(m), dtype=bool)
                    keep[k] = False
                    if (x := _rank_auc(sl[keep], m[keep])) is not None:
                        loo.append(x)
                row["auroc_gemini_loo_min"] = round(float(min(loo)), 4)
                row["auroc_gemini_loo_max"] = round(float(max(loo)), 4)
            top = joined.iloc[np.argsort(sl)[::-1][:max(50, 100, 500)]]
            for K in (50, 100, 500):
                row[f"{tag}_top{K}"] = int(top.head(K)["pair"].isin(pos).sum())
            # rank position sorted by strength descending (1 = most SL-like)
            joined[f"_order_{tag}"] = np.arange(len(joined)) + 1
            joined[f"rank_{tag}"] = pd.Series(sl).rank(
                ascending=False, method="average").to_numpy()

        # Spearman between score axes on the shared universe
        sf = pd.DataFrame(
            {"robust": -joined["robust_score"], "conv": -joined["conv_score"],
             "gemini": joined["gemini_strong_median"]})
        row["corr_robust_gemini"] = round(float(sf["robust"].corr(
            sf["gemini"], method="spearman")), 3)
        row["corr_conv_gemini"] = round(float(sf["conv"].corr(
            sf["gemini"], method="spearman")), 3)

        # audit layer: does the certificate decline GEMINI's top hits?
        gem_top = joined.nlargest(500, "gemini_strong_median")
        row["gemini_top500_certified"] = int(gem_top["certified"].sum())
        row["gemini_top500_declined"] = int((~gem_top["certified"]).sum())
        cert = joined[joined["certified"]]
        row["certified_frac_in_gemini_top500"] = round(
            float(cert["pair"].isin(set(gem_top["pair"])).mean()), 3) if len(cert) else None
        fdr_hits = joined[joined["gemini_fdr"] < 0.05]
        row["gemini_fdr_hits"] = int(len(fdr_hits))
        row["gemini_fdr_hits_certified"] = int(fdr_hits["certified"].sum())
        row["gemini_fdr_hits_declined"] = int((~fdr_hits["certified"]).sum())

        for (ga, gb) in sorted(VALIDATED_PARALOG_SL):
            hit = joined[(joined.gene_a == ga) & (joined.gene_b == gb)]
            if hit.empty:
                hit = joined[(joined.gene_b == ga) & (joined.gene_a == gb)]
            if not hit.empty:
                row[f"rank_gemini_{ga}-{gb}"] = int(hit.iloc[0]["rank_gemini"])
                row[f"gemini_strong_{ga}-{gb}"] = round(
                    float(hit.iloc[0]["gemini_strong_median"]), 3)
        covered = sorted(("-".join(sorted(t)) for t in pos & set(joined["pair"])))
        row["panel_covered_pairs"] = ";".join(covered)
        rows.append(row)

    bench = pd.DataFrame(rows)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    bench.to_csv(outdir / "gemini_benchmark.csv", index=False)
    (outdir / "gemini_benchmark.json").write_text(
        json.dumps({"gemini_benchmark_by_rho": bench.to_dict("records")}, indent=2)
    )
    print(bench.to_string(index=False))

    # self-check against stored external-panel values (rho 0.10)
    stored = [0.6575, 0.7752]
    got = [round(float(x), 4) for x in
           (bench.loc[bench.rho == 0.10, ["auroc_robust", "auroc_conv"]].iloc[0])]
    assert abs(got[0] - stored[0]) < 0.01 and abs(got[1] - stored[1]) < 0.01, (got, stored)
    print(f"self-check ok: auroc_robust={got[0]} auroc_conv={got[1]} match stored")


if __name__ == "__main__":
    main()