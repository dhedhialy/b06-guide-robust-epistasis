"""Render the pilot result figures from results/ CSVs (reproducible: nothing
derived here, only plotted)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def rho_scout(calls: pd.DataFrame, out: Path) -> None:
    resolved = calls[calls.gene_sign != "UNRESOLVED"].groupby("rho").size()
    unresolved = calls[calls.gene_sign == "UNRESOLVED"].groupby("rho").size()
    rhos = sorted(calls.rho.unique())
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    ax.bar(rhos, [resolved.get(r, 0) for r in rhos], width=0.05, color="#2e7d32", label="resolved")
    ax.bar(rhos, [unresolved.get(r, 0) for r in rhos], bottom=[resolved.get(r, 0) for r in rhos],
           width=0.05, color="#bdbdbd", label="unresolved")
    ax.set_xlabel("context drift rho"); ax.set_ylabel("gene pairs (44 total)")
    ax.set_title("rho-scout: resolvable gene pairs shrink as drift widens")
    ax.legend(); plt.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def holdout_curve(holdout: pd.DataFrame, out: Path) -> None:
    ho = holdout.sort_values("rho")
    robust = ho.groupby("rho")["robust_concordance"].mean()
    conventional = ho.groupby("rho")["conv_concordance"].mean()
    callable = ho.groupby("rho")["robust_concordance"].apply(lambda s: s.notna().sum())
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(callable.values, robust.values, "o-", color="#2e7d32", label="robust (lower/upper)")
    ax.plot(callable.values, conventional.values, "s--", color="#757575", label="conventional (point est.)")
    for r, x, y in zip(robust.index, callable.values, robust.values):
        ax.annotate(f"rho={r:.2f}", (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_xlabel("callable gene pairs (coverage)"); ax.set_ylabel("sign concordance")
    ax.set_title("guide-identity holdout: robust beats point estimates\nat matched coverage")
    ax.legend(); ax.set_ylim(0.5, 1.0); plt.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def sign_map(calls: pd.DataFrame, out: Path, max_pairs: int = 250) -> None:
    c = calls[calls.rho == 0.10].sort_values(
        ["gene_sign", "n_guide_pairs", "gene_a", "gene_b"],
        ascending=[True, False, True, True],
    )
    # 17K gene pairs can't be drawn as a strip; show the resolved calls only.
    resolved = c[c.gene_sign != "UNRESOLVED"].head(max_pairs)
    if len(resolved) < len(c) * 0.999:
        resolved = pd.concat(
            [resolved, c[c.gene_sign == "UNRESOLVED"].head(max(0, max_pairs - len(resolved)))]
        )
    c = resolved
    order = {"POSITIVE": 0, "NEGATIVE": 1, "UNRESOLVED": 2}
    cmap = {"POSITIVE": "#2e7d32", "NEGATIVE": "#c62828", "UNRESOLVED": "#cfd8dc"}
    fig, ax = plt.subplots(figsize=(8, 0.32 * len(c) + 1.0))
    for i, (_, row) in enumerate(c.iterrows()):
        ax.add_patch(plt.Rectangle((0, i), 1, 1, color=cmap[row.gene_sign], ec="none"))
        ax.text(1.02, i + 0.5, f"{row.gene_a}-{row.gene_b}", va="center", fontsize=8)
    ax.set_xlim(0, 1); ax.set_ylim(0, len(c)); ax.set_yticks([]); ax.set_xticks([])
    ax.set_title("gene-pair sign map at rho=0.10 (green=positive, red=negative, grey=unresolved)")
    ax.invert_yaxis(); plt.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def rho_star_dist(mantis: pd.DataFrame, out: Path) -> None:
    ax = mantis["rho_star_min"].dropna()
    fig, h = plt.subplots(figsize=(5.5, 3.5))
    h.hist(ax, bins=40, color="#1565c0")
    h.axvline(0.05, color="#ef6c00", ls="--", label="rho=0.05")
    h.axvline(0.10, color="#c62828", ls="--", label="rho=0.10")
    h.axvline(0.20, color="#4a148c", ls="--", label="rho=0.20")
    h.set_xlabel("rho* (min over guide pairs; larger = more drift-tolerant)")
    h.set_ylabel("gene pairs")
    h.set_title("perturbation-strength frontier: rho* distribution")
    h.legend(); plt.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results", help="results/ dir")
    args = ap.parse_args(argv)
    r = Path(args.results)
    fig_dir = r / "figures"
    fig_dir.mkdir(exist_ok=True)
    calls = pd.read_csv(r / "bounds" / "gene_pair_calls.csv")
    rho_scout(calls[calls.rho.isin([0.05, 0.10, 0.20])], fig_dir / "rho_scout.png")
    holdout_curve(pd.read_csv(r / "bounds" / "guide_holdout.csv"), fig_dir / "holdout_coverage_precision.png")
    sign_map(calls, fig_dir / "sign_map_rho10.png")
    rho_star_dist(pd.read_csv(r / "mantis" / "gene_pairs.csv"), fig_dir / "rho_star_distribution.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())