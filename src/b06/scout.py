"""Stage 0 + Stage 1 scout on the Burgold pilot library.

Usage: python3 -m b06.scout --design raw/MOESM3_design.xlsx
Writes audit + rho* artifacts under results/."""

from __future__ import annotations

import argparse
from pathlib import Path

from .audit import (
    audit_summary,
    gene_guide_counts,
    matched_guide_pairs,
)
from .confidences import block_bootstrap_interval, gene_interval
from .contrasts import RHO_GRID, gene_pair_calls, guide_conflicts, guide_pair_contrasts
from .holdout import evaluate_guide_holdout
from .io import load_colo1, load_design, mu_control
from .reporting import write_artifacts


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="raw/MOESM3_design.xlsx")
    ap.add_argument("--colo1", help="run on the COLO1 100K screen counts file")
    ap.add_argument("--lfc-rds", help="authors' colo1_plm.rds normalized LFC as phenotype (needs pyreadr)")
    ap.add_argument("--out", default="results")
    ap.add_argument("--n-boot", type=int, default=200)
    args = ap.parse_args(argv)

    if args.colo1:
        df = load_colo1(args.colo1, lfc_plm=args.lfc_rds)
        source = "COLO1 100K screen (HT29 vs plasmid, authors' nLFC)" if args.lfc_rds else "COLO1 100K screen (HT29 vs plasmid)"
    else:
        df = load_design(args.design)
        source = "Burgold pilot library"
    mu = mu_control(df)
    summary = audit_summary(df)
    matched = matched_guide_pairs(df, mu)
    contrasts = guide_pair_contrasts(matched)
    calls = {rho: gene_pair_calls(contrasts, rho) for rho in RHO_GRID}
    conflicts = {rho: guide_conflicts(contrasts, rho) for rho in RHO_GRID}
    intervals = {rho: gene_interval(contrasts, rho) for rho in RHO_GRID}
    if len(contrasts) >= 5:
        boot = block_bootstrap_interval(contrasts, rho=0.10, n_iter=args.n_boot, seed=0)
    else:
        boot = None
    holdout = {rho: evaluate_guide_holdout(contrasts, rho) for rho in RHO_GRID}

    out = Path(args.out)
    write_artifacts(out, summary, df, gene_guide_counts(df), matched, contrasts,
                    calls, conflicts, intervals, boot, holdout, source=source)


if __name__ == "__main__":
    main()