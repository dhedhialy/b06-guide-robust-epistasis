"""Theorem 5.2 / Proposition 5.1: analytic bounds vs brute-force enumeration."""

import math
import random

from b06.bounds import (
    NEGATIVE,
    POSITIVE,
    UNRESOLVED,
    classify,
    context_corrected_numerator,
    guide_pair_bounds,
    main_effects,
    symmetric_bounds,
)


def brute_force(m_a0, m_0b, m_ab, mu, rho, steps=200):
    d_a, d_b = main_effects(m_a0, m_0b, mu)
    lo, hi = float("inf"), float("-inf")
    for i in range(steps + 1):
        r_a = 1 - rho + 2 * rho * i / steps
        for j in range(steps + 1):
            r_b = 1 - rho + 2 * rho * j / steps
            v = context_corrected_numerator(m_ab, mu, d_a, d_b, r_a, r_b)
            lo = min(lo, v)
            hi = max(hi, v)
    return lo, hi


def test_symmetric_bound_matches_brute_force():
    rng = random.Random(1)
    for _ in range(15):
        m_a0, m_0b, m_ab, mu = (rng.uniform(-3, 3) for _ in range(4))
        rho = rng.uniform(0.0, 0.5)
        res = symmetric_bounds(m_a0, m_0b, m_ab, mu, rho)
        lo, hi = brute_force(m_a0, m_0b, m_ab, mu, rho)
        assert abs(res.lower - lo) < 1e-9, (res.lower, lo)
        assert abs(res.upper - hi) < 1e-9, (res.upper, hi)


def test_symmetric_matches_rectangle():
    """Theorem 5.2 is the rectangle of Prop 5.1 with r in [1-rho, 1+rho]."""
    rng = random.Random(2)
    for _ in range(15):
        m_a0, m_0b, m_ab, mu = (rng.uniform(-3, 3) for _ in range(4))
        rho = 0.3
        sym = symmetric_bounds(m_a0, m_0b, m_ab, mu, rho)
        rect = guide_pair_bounds(m_a0, m_0b, m_ab, mu, 1 - rho, 1 + rho, 1 - rho, 1 + rho)
        assert math.isclose(sym.lower, rect.lower, abs_tol=1e-12)
        assert math.isclose(sym.upper, rect.upper, abs_tol=1e-12)


def test_robustness_criterion():
    rng = random.Random(3)
    for _ in range(200):
        m_a0, m_0b, m_ab, mu = (rng.uniform(-2, 2) for _ in range(4))
        rho = 0.2
        res = symmetric_bounds(m_a0, m_0b, m_ab, mu, rho)
        if res.sign != UNRESOLVED:
            assert res.lower > 0 or res.upper < 0


def test_boundary_case_is_unresolved():
    """Effect-size stress test (Sec 21.4): |delta| == rho*(|d_a|+|d_b|)."""
    mu, rho = 0.0, 0.5
    d_a, d_b, gamma = 2.0, 1.0, 3.0  # |delta| = 5, rho*(|d_a|+|d_b|) = 1.5 != 5
    # construct delta exactly on the boundary
    half = rho * (abs(d_a) + abs(d_b))
    center = half  # |delta| == half -> interval [0, 2*half]
    m_a0, m_0b = mu + d_a, mu + d_b
    m_ab = center + m_a0 + m_0b - mu
    res = symmetric_bounds(m_a0, m_0b, m_ab, mu, rho)
    assert res.lower == 0.0
    assert res.sign == UNRESOLVED
    assert abs(res.rho_star - rho) < 1e-12


def test_sign_calls_match_interval():
    assert classify(1.0, 2.0) == POSITIVE
    assert classify(-2.0, -1.0) == NEGATIVE
    assert classify(-1.0, 1.0) == UNRESOLVED
    assert classify(0.0, 2.0) == UNRESOLVED
    assert classify(-2.0, 0.0) == UNRESOLVED