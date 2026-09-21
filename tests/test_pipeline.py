"""End-to-end pipeline on the real pilot library (regression)."""

from pathlib import Path

from b06.audit import audit_summary, matched_guide_pairs
from b06.confidences import block_bootstrap_interval, gene_interval
from b06.contrasts import gene_pair_calls, guide_pair_contrasts
from b06.io import load_design, mu_control

DESIGN = Path(__file__).resolve().parent.parent / "raw" / "MOESM3_design.xlsx"


def get_pipeline():
    df = load_design(str(DESIGN))
    mu = mu_control(df)
    return df, mu, matched_guide_pairs(df, mu)


def test_audit_plausible():
    df, _, _ = get_pipeline()
    s = audit_summary(df)
    assert s["guide_pairs_total"] > 8000
    assert s["guides_with_matched_single"] > 1000
    assert s["genes_with_2plus_guides"] >= 500


def test_matched_guide_pairs_produce_contrasts():
    df, mu, matched = get_pipeline()
    assert len(matched) >= 100
    ct = guide_pair_contrasts(matched)
    assert len(ct) >= 100
    assert (ct["delta"].notna()).all()
    assert (ct["rho_star"] >= 0).all()


def test_rho_monotonicity():
    """Larger drift envelope only ever moves calls to UNRESOLVED."""
    df, mu, matched = get_pipeline()
    ct = guide_pair_contrasts(matched)
    c05, c10, c20 = (gene_pair_calls(ct, r) for r in (0.05, 0.10, 0.20))
    n_resolved = lambda c: int((c.gene_sign != "UNRESOLVED").sum()) if not c.empty else 0
    assert n_resolved(c05) >= n_resolved(c10) >= n_resolved(c20)


def test_known_paralog_hit():
    """A genuinely reproducible interaction must still be called once the
    guide-idiosyncrasy floor is in place (old PRMT1-PRMT5 call was an artifact
    of dispersion the floor now correctly rejects)."""
    df, mu, matched = get_pipeline()
    ct = guide_pair_contrasts(matched)
    c = gene_pair_calls(ct, 0.10)
    row = c[(c.gene_a == "MTBP") & (c.gene_b == "MYC")].iloc[0]
    assert row.gene_sign == "POSITIVE"


def test_confidence_covered():
    df, mu, matched = get_pipeline()
    ct = guide_pair_contrasts(matched.iloc[:20])
    n_genes = len(gene_interval(ct, 0.10))
    lo, hi, frac = block_bootstrap_interval(ct, rho=0.10, n_iter=50, seed=0)
    assert len(lo) == n_genes
    assert all(0.0 <= f <= 1.0 for f in frac)
    # intersection may be empty (lo > hi) when guide pairs conflict -> not identified