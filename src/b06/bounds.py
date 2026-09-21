"""Closed-form partial-identification bounds for combinatorial CRISPR screens.

Implements the math-only core of the B06 plan. No data needed; the
experimental data only feeds the same scalar arguments.

* Theorem 4.1 (stable-efficacy guardrail): under the matched-guide bilinear
  model with fixed positive efficacy, the factorial contrast is exactly
  ``delta = s_a*s_b*gamma`` with ``s_a, s_b > 0``, so the interaction sign is
  unchanged.
* Proposition 5.1 (sharp guide-pair bounds): with context ratios confined to a
  rectangle, ``N(r_a, r_b) = m_ab - mu - r_a*d_a - r_b*d_b`` is affine, so its
  extrema sit at the rectangle corners and the interval is sharp.
* Theorem 5.2 (symmetric context drift): for ``r_a, r_b in [1-rho, 1+rho]`` the
  sign is identified iff ``|delta| > rho*(|d_a| + |d_b|)``.

All reported intervals are on the ``q_ab * gamma`` scale (``q_ab > 0``
product of double-context strengths); only the *sign* of ``gamma`` is
identified, which is what the intervals certify.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

POSITIVE = "POSITIVE"
NEGATIVE = "NEGATIVE"
UNRESOLVED = "UNRESOLVED"

_EPS = 1e-12


@dataclass(frozen=True)
class BoundResult:
    lower: float
    upper: float
    sign: str
    rho_star: float | None = None


def classify(lower: float, upper: float) -> str:
    """Sign call from an identified interval (boundary zeros are unresolved)."""
    if lower > 0:
        return POSITIVE
    if upper < 0:
        return NEGATIVE
    return UNRESOLVED


def delta_ab(m_a0: float, m_0b: float, m_ab: float, mu: float) -> float:
    """Conventional matched-guide factorial contrast: ``m_ab - m_a0 - m_0b + mu``."""
    return m_ab - m_a0 - m_0b + mu


def main_effects(m_a0: float, m_0b: float, mu: float) -> tuple[float, float]:
    """Single-guide effects relative to control: ``(d_a, d_b)``."""
    return m_a0 - mu, m_0b - mu


def context_corrected_numerator(
    m_ab: float, mu: float, d_a: float, d_b: float, r_a: float, r_b: float
) -> float:
    """``N_ab = m_ab - mu - r_a*d_a - r_b*d_b``; sign matches ``gamma``."""
    return m_ab - mu - r_a * d_a - r_b * d_b


def rho_star(m_a0: float, m_0b: float, m_ab: float, mu: float) -> float:
    """Robustness radius ``|delta| / (|d_a| + |d_b|)`` (Sec 5.4)."""
    d_a, d_b = main_effects(m_a0, m_0b, mu)
    denom = abs(d_a) + abs(d_b)
    if denom <= _EPS:
        return float("inf")
    return abs(delta_ab(m_a0, m_0b, m_ab, mu)) / denom


def symmetric_bounds(
    m_a0: float, m_0b: float, m_ab: float, mu: float, rho: float, noise: float = 0.0
) -> BoundResult:
    """Theorem 5.2: context ratios in ``[1-rho, 1+rho]`` giving the closed-form
    interval ``[delta - rho*(|d_a|+|d_b|), delta + rho*(|d_a|+|d_b|)]``.

    ``noise`` widens the interval by a guide-level idiosyncrasy floor
    (per-gene-pair MAD of guide-pair deltas) that the drift envelope does not
    capture; without it the bounds certify signs that guide noise alone can
    flip and the calls do not reproduce under guide identity holdout.
    """
    center = delta_ab(m_a0, m_0b, m_ab, mu)
    half = rho * sum(abs(d) for d in main_effects(m_a0, m_0b, mu)) + noise
    lower, upper = center - half, center + half
    return BoundResult(lower, upper, classify(lower, upper), rho_star(m_a0, m_0b, m_ab, mu))


def guide_pair_bounds(
    m_a0: float,
    m_0b: float,
    m_ab: float,
    mu: float,
    r_a_lo: float,
    r_a_hi: float,
    r_b_lo: float,
    r_b_hi: float,
) -> BoundResult:
    """Proposition 5.1: sharp interval over arbitrary rectangle of context
    ratios. Affine in ``(r_a, r_b)``, so extrema occur at rectangle corners."""
    d_a, d_b = main_effects(m_a0, m_0b, mu)
    corners = product((r_a_lo, r_a_hi), (r_b_lo, r_b_hi))
    vals = [context_corrected_numerator(m_ab, mu, d_a, d_b, ra, rb) for ra, rb in corners]
    return BoundResult(min(vals), max(vals), classify(min(vals), max(vals)))


def bound_interaction(
    singles: tuple[float, float],
    double: float,
    mu: float,
    rho: float,
) -> BoundResult:
    """Plan API entry point (Sec 16.1) for the symmetric context-drift model.

    ``singles = (m_a0, m_0b)``, ``double = m_ab``, ``mu`` = control mean.
    Not yet wired: efficacy_bounds and off_target_bounds.
    """
    return symmetric_bounds(singles[0], singles[1], double, mu, rho)