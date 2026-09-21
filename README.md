# B06 — Guide-Efficacy-Robust Epistasis

What combinatorial CRISPR screens can identify about genetic interaction when
perturbation strength is uncertain. Closed-form partial identification: the
output is an *identified set* for each interaction sign, not a point estimate.

## Status

Stage 0 + Stage 1 **run on real data** (Burgold et al. pilot library, HT-29):

| rho | gene pairs called | POSITIVE | NEGATIVE | UNRESOLVED |
|-----|-------------------|----------|----------|------------|
| 0.05 | 44 | 21 | 3 | 20 |
| 0.10 | 44 | 17 | 2 | 25 |
| 0.20 | 44 | 8 | 2 | 34 |

Scope: 81 gene pairs have a matched single; only 44 have >=2 matched guide
pairs (the multiplicity needed to certify robustness) — the single-guide-pair
designs are excluded from gene-level calls.

Sanity (corrected noise-floor model): only genuinely reproducible signals
survive the guide-idiosyncrasy floor — MTBP-MYC stays POSITIVE, PRMT1-PRMT5
(tested earlier as POSITIVE) is honestly UNRESOLVED because its 16 guide
deltas span both signs (floor 0.80 > effect). Audit verdict: the pilot lacks
the guide multiplicity needed for strong gene-level claims (214/3146 guides
bridge single<->double); the 100,136-pair Zenodo library is the primary
dataset for real claims.

**Guide-identity holdout (Stage 2, LOOCV by guide pair)** — with the
above-floor scoring contract the pilot retains only 6 gene pairs; it is
underpowered and its earlier "robust wins" spread was a small-sample
artifact:

| rho | callable pairs | robust concordance | conventional |
|-----|----------------|--------------------|--------------|
| 0.05 | 6 | 0.692 | 0.787 |
| 0.10 | 6 | 0.692 | 0.787 |
| 0.20 | 6 | 0.769 | 0.787 |
| 0.30 | 6 | 0.614 | 0.787 |

The 100K screen is the decisive test.

## COLO1 (100,136-pair screen) — run on real data

Source: `COLO1_RUNMERGED_EXACT_ANNOTATED.txt`, Zenodo 10.5281/zenodo.17191951
(open). 78,896 constructs; FC = log2(HT29/plasmid) from 3+3 replicate
triplets (triplet structure confirmed by pairwise correlation; orientation
by control depletion). 67,949 matched guide pairs reach the contrast stage —
the multiplicity the pilot lacked.

Intervals carry a data-driven guide-idiosyncrasy noise floor (one within-pair
MAD of guide-pair deltas); without it the bounds certify signs that guide
noise alone can flip and no longer reproduce (this was the root-cause of the
first failed holdout — PRMT1-PRMT5 is the canonical example).

| rho | gene pairs called | POSITIVE | NEGATIVE | UNRESOLVED |
|-----|-------------------|----------|----------|------------|
| 0.05 | 16999 | 127 | 1114 | 15758 |
| 0.10 | 16999 | 61 | 897 | 16041 |
| 0.20 | 16999 | 21 | 597 | 16381 |
| 0.30 | 16999 | 10 | 398 | 16591 |

Biology validates: all five paralog pairs the paper lists as validated GIs
(ARID1A-ARID1B, HDAC1-HDAC2, MAPK1-MAPK3, ASF1A-ASF1B, CNOT7-CNOT8) are called
NEGATIVE at rho=0.05 and 0.10.

**Honest headline — the certificate's value is selectivity, not higher
transfer.** Guide-identity holdout at the above-floor scoring contract (both
methods scored only where |delta| exceeds the guide-noise floor):

| rho | robust callable | robust | conventional (same pairs) |
|-----|-----------------|--------|---------------------------|
| 0.05 | 15779 | 0.566 | 0.620 |
| 0.10 | 15140 | 0.575 | 0.627 |
| 0.20 | 13422 | 0.593 | 0.643 |
| 0.30 | 11528 | 0.606 | 0.654 |

The naive point estimate transfers as well or slightly better on every
matched subset, so the Stage-2 over-claim (robust *predicts* held-out guides
better) is rejected by its own decisive test (KILL test #1). What the
certificate measurably buys: it refuses to certify ~5,500 of 17,000 gene
pairs (guide disagreement = false-positive control), and the certified set
alone lifts even the naive predictor from 0.620 (blanket) to 0.654
(certified) — reproducibility concentrates inside the robustness set. Robust
transfer (0.606) is comparable to conventional inside that set once
operating at rho>=0.20-0.30 rather than the over-lax grid bottom. The next
honest step, per the plan's Sec 14-21, is magnitude-adaptive rho and a
decision-rule that emits robust's certified-mean call — not more multiplicity.

**Position vs current SOTA (AUDIT Sec 6).** GEMINI (sellerslab/gemini, Bioc
1.24.0) run on this exact screen ranks the 4 covered De Kegel paralogs
better than B06 (panel AUROC **0.881**, LOO 0.844-0.968, vs robust 0.657 /
conv 0.775) and puts the validated paralogs at the top (CNOT7-CNOT8 #1,
ASF1A-ASF1B #2, MAPK1-MAPK3 #3, HDAC1-HDAC2 #8 of 18,920) — **B06 is not
SOTA on ranking.** But GEMINI emits **zero FDR-significant hits on this
screen** (min-FDR floors at 0.83 across the library) and rank-agrees with
B06 only weakly (Spearman 0.15-0.22); the certificate's audit-layer value
stands, and per-pair GEMINI vs certificate agreement is the next check.
`benchmarks/gemini/run_gemini.R` + `bench_gemini.py` are pinned.

Re-run: `PYTHONPATH=src python3 -m b06.scout --colo1
raw/colo1_release/.../COLO1_RUNMERGED_EXACT_ANNOTATED.txt --out results_colo1
--n-boot 200`; figures in `results_colo1/figures/`.

## Results (figures)

`results/figures/` (regenerate: `PYTHONPATH=src python3 -m b06.figures`):

| figure | shows |
|--------|-------|
| `rho_scout.png` | resolved vs unresolved gene pairs as rho widens |
| `holdout_coverage_precision.png` | robust vs conventional transfer at matched coverage (Fig 5 style) |
| `sign_map_rho10.png` | per-gene-pair sign at rho=0.10 |
| `rho_star_distribution.png` | how many pairs sit on the perturbation-strength frontier |

## Modules

```
data_contract/   EPISTASIS_ID_LOCK, GUIDE_MULTIPLICITY_LOCK (emitted), ASSUMPTIONS, KILL_TESTS
raw/             Burgold Supplementary Data 1 (MOESM3, click-through license)
src/b06/
  bounds.py      closed-form guardrail + sharp bounds + rho* (no deps)
  io.py          pilot-library parsing
  audit.py       guide multiplicity, matched singles, coverage
  contrasts.py   matched-guide deltas + gene-level calls
  confidences.py gene-level intersection + guide-level bootstrap
  holdout.py     guide-identity LOOCV evaluator (robust vs conventional)
  validate.py    audit: LOOCV ranking-AUC, validated-paralog rank test, selectivity ablation, NOISE_K sweep, De Kegel external-panel test (LOO-stable AUROC), certified-count + Sec-5.1 call-count permutation-FDR null, model-based synthetic-screen null
  figures.py     reproducible result figures from results/ CSVs
  scout.py       pipeline entry point (--lfc-rds = authors' colo1_plm.rds phenotype)
  provenance.py  sha256 + versions + seed audit manifest
tests/           guardrail, symmetric bound, sharpness, end-to-end, monotonicity, holdout, synthetic io regression, external-panel verdict + LOO, permutation-FDR nulls
notebooks/       01_rho_scout.ipynb, 02_guide_holdout.ipynb (Mantis Coding / Jupyter)
`results/`         audit + bounds CSVs + figures + Mantis-import space (results/mantis/gene_pairs.csv)
`results_colo1/`   same pipeline outputs on the COLO1 100K screen (HT29, Zenodo)
`results_colo1_plm/` phenotype cross-check run (authors' nLFC, Sec 5.6 of AUDIT.md)
`benchmarks/gemini/`  GEMINI driver + benchmark join (Sec 6; needs framework R 4.5.3 arm64 + Bioc gemini 1.24.0)
`AUDIT.md`         master audit: full methodology + all results + known gaps
```

## Use

```bash
PYTHONPATH=src python3 -m b06.scout --design raw/MOESM3_design.xlsx --n-boot 200
PYTHONPATH=src python3 -m pytest   # 25 tests
PYTHONPATH=src python3 benchmarks/gemini/bench_gemini.py   # Sec 6 SOTA benchmark join
```

```python
from b06.bounds import bound_interaction
bound_interaction(singles=(m_a0, m_0b), double=m_ab, mu=mu, rho=0.1)
# BoundResult(lower, upper, sign, rho_star); sign in POSITIVE|NEGATIVE|UNRESOLVED
```

## Mantis merge path

`results/mantis/gene_pairs.csv` (one row per gene pair: rho*, sign at
rho=0.10, n guide pairs, confidence from `guide_bootstrap.csv`) imports as a
Mantis space via the CSV importer; the notebook runs in Mantis Coding;
locks + summaries live as workspace files. MantisAPI (this repo's backend)
already exposes space creation/ingestion endpoints for a full programmatic
push once the 100K library is ingested.

## Next

1. Second independent caller (GRAPE/Orthrus) + per-pair joint-agreement
   between the certificate and GEMINI hits (GEMINI itself: done, AUDIT Sec 6).
2. External-panel upgrade: the full ~126-pair De Kegel training set (gated,
   needs a DepMap 20Q2 pull + four combinatorial-screen tables) with a
   paralog-matched negative background to power the ranking test past n=4.
3. ~~Model-based null (GEMINI-style synthetic screens) for the certified-count
   claim, beyond the permutation null (Sec 5.1/7).~~ **Done** — Sec 5.10:
   within-block `d_a`/`d_b` correlation measured −0.042; the FCR reading is
   delta-model-dependent (permutation 779 vs empirical-marginal 4,208 vs
   symmetric-normal 2,744 of 4,545 observed) and is reported as a bracket
   without restating the headline; the certificate's interaction content is
   corroborated independently by the GEMINI cross-check.
4. ~~Seed-encoded off-target sensitivity using the sgRNA_Off_Target
   columns.~~ **Closed as a data dead-end** — no `sgRNA_Off_Target`/guide
   metadata ships in this release (verified on the uploaded raw file). Would
   need re-acquiring annotations from the depositor. (FC orientation vs
   `colo1_plm.rds` was confirmed in Sec 5.6; counts-derived FC remains
   primary.)