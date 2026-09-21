"""Writes scout artifacts: CSVs, audit summary, and a Mantis-import space."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def write_artifacts(out: Path, summary, df, guides, matched, contrasts, calls, conflicts, intervals=None, boot=None, holdout=None, source="Burgold pilot library") -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "audit").mkdir(parents=True, exist_ok=True)
    (out / "bounds").mkdir(parents=True, exist_ok=True)
    (out / "mantis").mkdir(parents=True, exist_ok=True)

    with (out / "audit" / "audit_summary.json").open("w") as f:
        json.dump(summary, f, indent=2)

    guides.to_csv(out / "audit" / "guide_usage.csv", index=False)
    df[["guide_pair_id", "guide_a", "gene_a", "class_a", "guide_b", "gene_b",
        "class_b", "scaffold", "notes", "vector_class", "FC_500x"]].to_csv(
        out / "audit" / "guide_pairs_design.csv", index=False)

    if not matched.empty:
        matched.to_csv(out / "audit" / "matched_guide_pairs.csv", index=False)
    if not contrasts.empty:
        contrasts.to_csv(out / "bounds" / "guide_pair_contrasts.csv", index=False)

    if intervals:
        iv = pd.concat(
            [iv.assign(rho=rho) for rho, iv in intervals.items() if not iv.empty],
            ignore_index=True,
        )
        iv.to_csv(out / "bounds" / "gene_intervals.csv", index=False)

    if holdout:
        ho = pd.concat(
            [h.assign(rho=rho) for rho, h in holdout.items() if not h.empty],
            ignore_index=True,
        )
        ho.to_csv(out / "bounds" / "guide_holdout.csv", index=False)
    if boot is not None:
        lo, hi, frac = boot
        giv = intervals[0.10].sort_values(["gene_a", "gene_b"]).reset_index(drop=True)
        pd.DataFrame(
            {
                "gene_a": giv["gene_a"],
                "gene_b": giv["gene_b"],
                "boot_ci_lower": lo,
                "boot_ci_upper": hi,
                "boot_robust_frac": frac,
            }
        ).to_csv(out / "bounds" / "guide_bootstrap.csv", index=False)

    gene_wide = []
    for rho, c in calls.items():
        if c.empty:
            continue
        c2 = c.copy()
        c2["rho"] = rho
        gene_wide.append(c2)
    if gene_wide:
        all_calls = pd.concat(gene_wide)
        all_calls.to_csv(out / "bounds" / "gene_pair_calls.csv", index=False)

        # Mantis space: one point per gene pair, rho=0.10 layer
        base = calls[0.10].copy()
        star = (
            contrasts.groupby(["gene_a", "gene_b"])["rho_star"]
            .agg(["mean", "min", "max"])
            .reset_index()
        )
        geom = base.merge(star, on=["gene_a", "gene_b"], how="outer")
        geom["rho_star_mean"] = geom["mean"]
        geom["rho_star_min"] = geom["min"]
        geom["rho_star_max"] = geom["max"]
        geom["sign_10"] = geom["gene_sign"].fillna("UNRESOLVED")
        geom["resolved"] = geom["sign_10"] != "UNRESOLVED"
        geom = geom.sort_values(["resolved", "rho_star_mean"], ascending=[False, False])
        geom = geom[
            ["gene_a", "gene_b", "rho_star_mean", "rho_star_min", "rho_star_max",
             "sign_10", "n_guide_pairs"]
        ]
        geom.to_csv(out / "mantis" / "gene_pairs.csv", index=False)

    rows = []
    for rho, c in conflicts.items():
        if not c.empty:
            rows.append(c.assign(rho=rho))
    if rows:
        pd.concat(rows).to_csv(out / "bounds" / "guide_conflicts.csv", index=False)

    write_summary(out, summary, matched, contrasts, calls, holdout, source)


def write_summary(out: Path, summary, matched, contrasts, calls, holdout=None, source="Burgold pilot library") -> None:
    txt = [
        f"# Stage 0/1 scout summary ({source})\n",
        f"- guide pairs: {summary['guide_pairs_total']}",
        f"- unique guides: {summary['unique_guides']}",
        f"- genes: {summary['genes']} (2+ guides: {summary['genes_with_2plus_guides']})",
        f"- gene pairs observed: {summary['gene_pairs_observed']} "
        f"(full cross-product: {summary['gene_pairs_full_cross_product']})",
        f"- NTC control pairs: {summary['nontargeting_control_pairs']}, "
        f"control phenotype median: {summary['control_phenotype_median']:.3f}",
        f"- guides with a matched single control: {summary['guides_with_matched_single']}",
        f"- guides bridging single and double contexts: {summary['guides_bridging_single_and_double']}",
        "",
        f"- matched guide pairs computable: {len(matched)}",
        f"- guide-pair contrasts computable: {len(contrasts)}",
        "",
        "### Sign calls vs rho (gene-level, resolvable pairs only)\n",
        "| rho | gene pairs with calls | positive | negative | unresolved | guide conflicts |",
        "|-----|----------------------|----------|----------|------------|-----------------|",
    ]
    for rho, c in calls.items():
        n = len(c)
        pos = int(c.gene_sign.eq("POSITIVE").sum()) if not c.empty else 0
        neg = int(c.gene_sign.eq("NEGATIVE").sum()) if not c.empty else 0
        unres = int(c.gene_sign.eq("UNRESOLVED").sum()) if not c.empty else 0
        txt.append(f"| {rho:.2f} | {n} | {pos} | {neg} | {unres} | |")
    if holdout:
        txt.append("\n### Guide-identity holdout (LOOCV by guide pair)\n")
        txt.append("| rho | gene pairs | robust callable | robust concordance | conventional (same pairs) |")
        txt.append("|-----|------------|-----------------|--------------------|---------------------------|")
        for rho in sorted(holdout):
            h = holdout[rho]
            rob = h.dropna(subset=["robust_concordance"])
            cvg = h.dropna(subset=["conv_concordance"])
            if rob.empty:
                txt.append(f"| {rho:.2f} | {len(h)} | 0 | - | - |")
                continue
            rc = rob["robust_concordance"].mean()
            cc = cvg.loc[rob.index, "conv_concordance"].mean()
            txt.append(f"| {rho:.2f} | {len(h)} | {len(rob)} | {rc:.3f} | {cc:.3f} |")
    (out / "scout_summary.md").write_text("\n".join(txt) + "\n")