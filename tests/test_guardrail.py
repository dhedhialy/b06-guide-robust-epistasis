"""Theorem 4.1: fixed positive matched-guide efficacy cannot flip the sign."""

import math
import random

from b06.bounds import delta_ab


def test_stable_efficacy_cancellation():
    rng = random.Random(0)
    mu = 1.0
    for _ in range(500):
        sa = rng.uniform(0.01, 5.0)
        sb = rng.uniform(0.01, 5.0)
        alpha_a = rng.uniform(-2.0, 2.0)
        alpha_b = rng.uniform(-2.0, 2.0)
        gamma = rng.uniform(-3.0, 3.0)

        m_a0 = mu + sa * alpha_a
        m_0b = mu + sb * alpha_b
        m_ab = mu + sa * alpha_a + sb * alpha_b + sa * sb * gamma

        d = delta_ab(m_a0, m_0b, m_ab, mu)
        assert math.isclose(d, sa * sb * gamma, abs_tol=1e-9)
        if abs(gamma) > 1e-9:
            assert (d > 0) == (gamma > 0), "sign flipped by stable efficacy"


def test_zero_interaction_stays_zero():
    mu = 2.0
    d = delta_ab(mu + 1, mu - 1, mu, mu)
    assert d == 0.0