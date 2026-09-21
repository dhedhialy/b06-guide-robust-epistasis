# ASSUMPTIONS

Every model used here is only claimable under these. Extending them = new record.

1. Matched-guide bilinear model (Sec 4.1): `m_ab = mu + s_a*alpha + s_b*beta + s_a*s_b*gamma`, `s_a, s_b > 0`.
2. Linear-in-effective-dose main effects (Sec 5.1): `m_ab = mu + r_a*d_a + r_b*d_b + q_ab*gamma`, `q_ab > 0`.
3. Context ratios bounded: `r_a, r_b in [1-rho, 1+rho]` (symmetric) or a stated rectangle (general).
4. No resolution of guide-specific direct/off-target effects (future layer, Sec 12).
5. No nonlinear dose response (future layer, Sec 13).
6. Population means only; sampling uncertainty folded in separately (Sec 8).