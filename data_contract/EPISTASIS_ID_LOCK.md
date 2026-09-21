# EPISTASIS_ID_LOCK

Fixed at first scout (pilot phase); change = new record.
Emitted 2026-09-04.

| # | Field | Value |
|---|-------|-------|
| 1 | Interaction estimand | Gene-level genetic interaction on the matched-dual-guide factorial contrast `delta = m_ab - m_a0 - m_0b + mu`, per gene `log2` fold-change phenotype |
| 2 | Phenotype | `FC_500x` (pooled viability, log2 fold change D14 vs D3), pilot screen, HT-29 |
| 3 | Control definition | non-targeting x non-targeting guide pairs; `mu` = median FC_500x |
| 4 | Sign convention | POSITIVE = synthetic sick/lethal; NEGATIVE = recovery/suppression; UNRESOLVED else |
| 5 | Guides matched singles<->doubles | Same guide identity with inert partner = single context; guide pairs where either side lacks a matched single are excluded |
| 6 | Allowable efficacy uncertainty | Retained as symmetric context-drift envelope rho in {0.05, 0.10, 0.20} |
| 7 | Allowable context drift | r_a, r_b in [1-rho, 1+rho] (Theorem 5.2) |
| 8 | Off-target assumptions | Not yet modeled; per-guide `sgRNA_Off_Target` metadata available for the sensitivity layer |
| 9 | Nonlinear-dose assumptions | Linear-in-effective-dose model; nonlinear sensitivity deferred (Sec 13) |