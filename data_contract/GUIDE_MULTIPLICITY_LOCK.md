# GUIDE_MULTIPLICITY_LOCK — Burgold pilot library

Dataset: Burgold et al., Nat Commun 17:561 (2026), Supplementary Data 1
(MOESM3). Pilot design + counts, HT-29 cells, dual-guide tRNA system.
Emitted 2026-09-04 from `b06.scout`.

| # | Field | Value |
|---|-------|-------|
| 1 | Dataset / screen | Burgold pilot (GI_PILOT_LIB2), HT-29 |
| 2 | Guide pairs (rows with FC) | 8,869 |
| 3 | Unique guide identities | 3,146 |
| 4 | Genes targeted | 538 (all with >= 2 guides) |
| 5 | Gene pairs observed | 1,249 |
| 6 | Gene pairs, complete guide cross-product | 521 |
| 7 | Guides with a matched single control | 1,071 |
| 8 | Guides bridging single and double contexts | 214 |
| 9 | Matched guide pairs usable for B06 | 238 |
| 10 | Gene pairs entering the rho* scout | 81 |
| 11 | Non-targeting control pairs | 595 (median FC_500x = 0.949) |
| 12 | Replicates | 500x/100x/PCR500x, D3 and D14, 3 reps |

## Audit verdict

Pilot has enough structure to *demonstrate* the method but **not** enough
guide multiplicity / matched singles to make strong gene-level claims
(only 214/3,146 guides bridge single and double) — exactly the
GUIDE_MULTIPLICITY_LOCK failure the plan (Sec 20.2) expects. The large-scale
100,136-pair library (EGA EGAD00001015754, access-controlled) is the primary
dataset for real claims; its design scan remains to be reconstructed.