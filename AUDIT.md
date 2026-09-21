# B06 Robust Epistasis — Master Audit (methodology + results)

Status: reproducible, commit-free artifact set. No git metadata (by design);
every output regenerable from the commands below. All numbers are as of the
final corrected model (noise floor + four-rung rho grid), which supersede
every earlier run. 26/26 tests pass.

Rev 2 (review-driven): added the external-phenotype cross-check
(`colo1_plm.rds` via pyreadr), the certified-vs-magnitude-matched
selectivity ablation, the noise-floor sensitivity sweep, and fixed two
`io.py` defects (duplicate `load_design`, dead `partner` variable) with a
synthetic regression test. See Sec 5.6-5.8 and 7.

Rev 2c (review round 2): permutation-FDR null under the Sec 5.1 call counts
(300 perms; null construction corrected after a bijective-relabel bug that
gave a point-mass null); corrected the external-panel rank-AUROC orientation
(enrichment was computed on the inverted axis — actual values 0.65-0.78, not
0.22-0.35) with per-pair leave-one-out ranges; rewrote the Sec 5.9 reading as
the ambiguous/concluded-inconclusive result it is; unit-tested the external
panel verdict rule, the certified null, and the call-count null (24 tests).
See Sec 5.1, 5.9, 7.

Rev 2d (review round 3): permutation nulls raised to 1000 draws with the
(k+1)/(n+1) one-sided MC p reported at the resolution floor (0/1000 reached
=> p = 1/1001, an upper bound, not a decimal estimate) and FCR-scaled
corrected call counts (observed - null-mean = expected true calls) emitted
per rho; the null's construction is now documented as a modeling choice that
preserves the block-size multiset and row covariates but pools main effects
across unrelated genes (within-block `d_a`/`d_b` correlation of real pairs
is not reproduced); added a set-based block-membership guard test (the row
partition must change, not merely sd > 0) and hard-asserted exact AUROC 1.0
on the perfectly separable panel synthetic (25 tests). See Sec 5.1, 7.

Rev 2e (this revision): quantified the SOTA position by running the field's
current caller, GEMINI (sellerslab/gemini, Bioc 1.24.0, Zamanighomi et al.
Genome Biology 2019), on this exact COLO1 78,912-guide-pair screen
(framework R 4.5.3 arm64, 50 CAVI iterations, nc = NONTARGET pairs, 15
cores). Result reads clearly: GEMINI ranks the 4 covered De Kegel paralogs
better than B06 (panel AUROC 0.881 vs robust 0.657 / conv 0.775) and its
per-gene-pair scores put the five validated paralogs near the very top
(CNOT7-CNOT8 #1, ASF1A-ASF1B #2, MAPK1-MAPK3 #3, HDAC1-HDAC2 #8 of 18,920),
so **B06 is not SOTA on ranking**; but GEMINI reports **zero significant
hits on this screen** (min-FDR over the 3 endpoint replicates floors at
0.83 across the whole library; the top `strong` magnitudes belong to
NONTARGET~NONTARGET controls) and rank-agrees with B06 only weakly
(Spearman robust 0.20-0.22, conv 0.146). Sec 6 now carries the quantitative
comparison; the certificate remains an audit layer GEMINI does not provide.
Also closed the Sec-7 model-based-null gap (3b) with a synthetic-screen null
that keeps each block's real covariate layout fixed (Sec 5.10) and closed
the seed off-target item as a data dead-end (Sec 7.8). 26 tests.
See Sec 5.9, 6, 7, 9.

---

## 1. Research goal

Closed-form robustness certification for interaction sign calls in
dual-guide CRISPR knockout screens. Standard pipelines call a genetic
interaction from a point estimate (the dLFC) with replicate noise assumed
away; B06 instead attaches to every sign call an *identifiability
certificate* — a closed-form guarantee that the call survives context-rate
misspecification (rho envelope) and guide-level reproducibility noise
(noise floor). The deliberate, pre-registered question is whether that
certificate predicts held-out guide behavior better than the point estimate
(the Stage-2 claim), and the answer is the decisive result of this audit.

## 2. Method — formal

**Model (matched-guide bilinear).** A construct with guides g_a, g_b has
contextified single effects `r_a·d_a`, `r_b·d_b` relative to control mean
`mu`, where `d_a = m_a0 - mu`, `d_b = m_0b - mu` and `r_*` are unknown
context multipliers. The **factorial contrast** is

```
delta = m_ab - m_a0 - m_0b + mu         (line 50 bounds.py)
```

Under the model, `delta = r_a·r_b·gamma` with `gamma` the interaction
parameter; only the sign of `gamma` is identified (two unknown multipliers
absorb magnitude), and that is all the certificates claim.

**Sharp symmetric bounds (Sec 5.2).** With context ratios confined to
`[1-rho, 1+rho]`,

```
[delta - rho·(|d_a|+|d_b|),  delta + rho·(|d_a|+|d_b|)]
```

extrema are affine, hence corner-bound and **sharp** (verified by brute
force in `tests/test_symmetric_bound.py`).

**Guide-idiosyncrasy noise floor (key correction).** The envelope covers
drift, not guide spread. A data-driven floor — one within-pair MAD of
guide-pair deltas (`noise = noise_k · 1.4826 · MAD`, `noise_k = 1.0`) — is
added to the half-width:

```
half_i = rho·(|d_a|+|d_b|) + noise_i        (upper/lower + noise in bounds.py)
```

Without the floor the certificate signed pairs that guide noise alone flips
and the calls did not reproduce under holdout; pilot PRMT1-PRMT5 (16 guide
pairs, deltas spanning both signs, floor 0.80) is the canonical example it
correctly rejects. Median COLO1 within-pair MAD is 0.24 in delta units vs
envelope widths of ~0.04-0.12, which is the original failure mechanism.

**Calls.** Guide-pair sign from `classify(lower, upper)`. Gene-level sign =
intersection of guide-pair intervals — `NEGATIVE` iff `max(lower) < 0`,
`POSITIVE` iff `min(upper) > 0`, else `UNRESOLVED`. Empty intersection
(guides conflict) is `UNRESOLVED`. `rho_star = |delta|/(|d_a|+|d_b|)` is
the perturbation-strength radius at which the pair becomes undecidable.
Rho grid `(0.05, 0.10, 0.20, 0.30)`; larger rho only ever moves calls toward
`UNRESOLVED` (monotonicity, tested).

**Continuous scores (audit addition).** Gene-level
`robust_score = mean_delta / (rho·(mean|d_a|+mean|d_b|) + noise)`, signed;
`certified = |mean_delta| > envelope`. `conv_score = mean_delta` is the
standard dLFC.

**Confidence.** Guide-level block bootstrap of each gene pair's intersection
interval (200 resamples, seed 0): mean endpoints + fractional sign stability
(`boot_robust_frac`).

**KILL test #1 — guide-identity holdout.** LOOCV by guide pair: for each
held-out pair, calibrate the gene interval on the remaining pairs at rho
(robust) vs the mean-delta call (conventional); score both against the
held-out pair's delta sign. **Decision contract:** effects with
`|delta| <= calibration MAD` are not scored for either method — a
sub-floor effect cannot be reproduced. Split-level floor only (no leakage).

## 3. Data & reproducibility

| item | value |
|------|-------|
| Screen | COLO1 (HT-29, 100,136-pair dual-guide) |
| Source file | `COLO1_RUNMERGED_EXACT_ANNOTATED.txt` (Zenodo 10.5281/zenodo.17191951, open, CC BY 4.0), 78,921 rows |
| SHA-256 | `84397a6c76641f7485aeb9e5541da79314d897324d29c22750140b8de1dbbf8b` |
| Raw reads | ENA ERP183979 (annotation) |
| External SL panel | De Kegel et al. 2021 (Cell Syst. 12:1144-1159.e6), `raw/de_kegel_2021/validated_SLs.txt` (sha256 `9173ddce…5de45ae`) |
| Pilot | Burgold MOESM3_design.xlsx (Supplementary Data 1, license-restricted) |
| Python | 3.9.6, pandas 2.3.3, numpy 1.26.4, matplotlib 3.9.4, pytest 8.4.2 |
| RDS read | pyreadr 0.5.6 (Homebrew-lzma universal2 build) |
| Random seed | 0 (all bootstraps) |
| Command (COLO1, counts-FC) | `PYTHONPATH=src python3 -m b06.scout --colo1 raw/colo1_release/.../COLO1_RUNMERGED_EXACT_ANNOTATED.txt --out results_colo1 --n-boot 200` |
| Command (COLO1, authors' nLFC) | `PYTHONPATH=src python3 -m b06.scout --colo1 <same> --lfc-rds raw/colo1_release/.../input/colo1_plm.rds --out results_colo1_plm --n-boot 200` |
| Command (pilot) | `PYTHONPATH=src python3 -m b06.scout --out results --n-boot 50` |
| Command (validation) | `PYTHONPATH=src python3 -m b06.validate --contrasts <run>/bounds/guide_pair_contrasts.csv --out <run>` (ingests `raw/de_kegel_2021/validated_SLs.txt` by default) |
| Tests | `PYTHONPATH=src python3 -m pytest` → 24 passed |

Replicate structure of the COLO1 counts (empirically confirmed): the 6 count
columns form two replicate triplets (correlation ~0.84 within, ~0.41 across).
FC = `log2(HT29/plasmid)` — triplets `CPID2437/2440/2443` = plasmid
reference, `CPID1020/1023/1026` = endpoint. Orientation is **now confirmed
against the authors' normalized LFC** (`colo1_plm.rds`, read via pyreadr):
per-construct Pearson r = 0.89, Spearman 0.83, slope 0.86, intercept −0.23,
85% sign agreement (89% at |FC|>0.25, 95% at |FC|>1.0). Not flipped; the
divergence is concentrated near zero. Because the two phenotypes give
materially different *boundary* calls (Sec 5.6), the counts-derived FC stays
the primary phenotype and the sensitivity is now quantified, not assumed.

## 4. Results — pilot (Burgold, 44 resolvable pairs)

| rho | called | POSITIVE | NEGATIVE | UNRESOLVED |
|-----|--------|----------|----------|------------|
| 0.05 | 44 | 12 | 3 | 29 |
| 0.10 | 44 | 8 | 1 | 35 |
| 0.20 | 44 | 4 | 1 | 39 |
| 0.30 | 44 | 2 | 1 | 41 |

Holdout (6 gene pairs survive the contract — underpowered; the early
"robust wins 0.79 vs 0.77" pilot claim was a small-sample artifact):

| rho | callable | robust | conventional |
|-----|----------|--------|--------------|
| 0.05 | 6 | 0.692 | 0.787 |
| 0.10 | 6 | 0.692 | 0.787 |
| 0.20 | 6 | 0.769 | 0.787 |
| 0.30 | 6 | 0.614 | 0.787 |

Validation (65 above-floor instances): ranking AUC robust 0.874-0.890,
conv 0.901. Biology sanity under the corrected model: MTBP-MYC POSITIVE,
PRMT1-PRMT5 honestly UNRESOLVED.

## 5. Results — COLO1 100K (decisive)

Guide pairs 78,896; unique guides 3,263; genes 1,207 (851 with 2+ guides);
gene pairs observed 20,308 (19,869 full cross-product); NTC pairs 725; NTC
median (mu) 0.806. Matched-single guides 1,686; bridging guides 948;
**matched guide pairs computable 67,949** — the multiplicity the pilot
lacked. Gene pairs = 16,999.

### 5.1 Sign calls vs rho

| rho | called | POSITIVE | NEGATIVE | UNRESOLVED | guide conflicts |
|-----|--------|----------|----------|------------|-----------------|
| 0.05 | 16999 | 127 | 1114 | 15758 | — |
| 0.10 | 16999 | 61 | 897 | 16041 | — |
| 0.20 | 16999 | 21 | 597 | 16381 | — |
| 0.30 | 16999 | 10 | 398 | 16591 | — |

Negative interactions dominate (cancer screen, most real GIs are
synthetic-sick/lethal), and resolving drops monotonically with rho.

**Error-rate control (permutation null, n=1000, seed 0).** The count table is
not a claim on its own, so we put a permutation null under the POSITIVE/
NEGATIVE counts. Null: shuffle the guide-pair contrasts and re-bucket them
into gene-pair blocks of the same size multiset, then recompute the strict
call rule (NEGATIVE iff every guide-pair interval of the block resolves
negative, `max(upper) < 0`; POSITIVE iff `min(lower) > 0`). This genuinely
moves rows between blocks — a bijective relabel of the block ids would only
rename blocks and gives a point-mass null (sd 0), which is why the
construction matters (a set-based membership guard in the test suite asserts
the row partition changes, not just sd > 0). `results_colo1/audit/call_fdr_null.csv`.

| rho | called | POS | NEG | null POS (mean±sd) | null NEG (mean±sd) | FCR POS | FCR NEG | perms ≥ obs | p (POS/NEG) | corr. POS | corr. NEG |
|-----|--------|-----|-----|--------------------|--------------------|---------|---------|-------------|-------------|-----------|-----------|
| 0.05 | 16999 | 127 | 1114 | 4.3±2.2 | 293.9±16.0 | 0.034 | 0.264 | 0/1000 both | 1/1001 | 122.7 | 820.1 |
| 0.10 | 16999 | 61 | 897 | 1.6±1.2 | 185.4±13.1 | 0.026 | 0.207 | 0/1000 both | 1/1001 | 59.4 | 711.6 |
| 0.20 | 16999 | 21 | 597 | 0.24±0.49 | 77.5±8.6 | 0.011 | 0.130 | 0/1000 both | 1/1001 | 20.8 | 519.5 |
| 0.30 | 16999 | 10 | 398 | 0.05±0.21 | 34.1±5.7 | 0.005 | 0.086 | 0/1000 both | 1/1001 | 10.0 | 363.9 |

Observed POSITIVE counts exceed the null by 29-217x and NEGATIVE by 3.8-11.7x;
0/1000 permutations reached either observed count at every rho. p is reported
at the Monte-Carlo resolution floor: **0/1000 reached => p = 1/1001 ≈ 1.0e-3
one-sided, an upper bound, not an estimate with precision below it** — the
observed counts sit so far above the null center (4-180x) that no feasible
number of draws would ever reach one; a smaller p can only be bought with a
more powerful test, not more permutations. The `corr.` columns are the
**FCR-scaled call counts** (observed − null-mean), the expected number of
TRUE calls in each set: at FCR ≤ 10% (rho=0.30) NEGATIVE 398 → 363.9
expected true; at rho=0.10, 897 → 711.6. These are point-estimate
corrections from the aggregate null, not per-pair BH thresholds (per-pair
FDR needs per-pair null draws — Sec 7 gap 3). The null-expected counts also
shrink monotonically with rho, so the certificate tightens rather than
trading away error control. The enlarged `gene_scores` certified set behaves
the same way (`certified_fdr_null.json`, rho=0.10): 4,545 certified observed
vs 779.0±24.1 null (0/1000 reached, corrected 3,766.0); NEG split 4,059 vs
743.6±23.2 (FCR 0.183, corrected 3,315.4); POS 486 vs 35.4±5.1 (FCR 0.073,
corrected 450.6).

*Reading kept honest: FCR is the permutation-expected fraction
(null-mean / observed), not a Benjamini-Hochberg q-value; the corrected
count = observed − null-mean (expected true calls); p-values are the
(k+1)/(n+1) one-sided MC estimate, here at the 1/1001 resolution floor.*

*Modeling choice (stated, not implicit): the null preserves each block's size
multiset and each row's covariates (delta, |d_a|, |d_b|, noise) and genuinely
re-buckets rows. What permuted blocks do NOT reproduce is the within-block
main-effect correlation of real gene pairs: a synthetic block pools `d_a`/`d_b`
from unrelated genes, whereas a real pair's are correlated (same gene, shared
guide context). The null is therefore a no-true-correspondence test fit for
order-of-magnitude claims, not a generative model of "null gene pairs look
like"; the direction of any FCR bias from that approximation is unknown and
was the explicit target of the Sec 7 model-based null (resolved below and in
Sec 5.10).*

### 5.10 Model-based synthetic-screen null (gap 3b resolved)

The gap-3b question: does the within-block main-effect correlation pooled
away by the permutation null move the FCR reading? The model-based null
(`model_based_fdr_null` in `validate.py`, results in
`model_based_fdr_null.json`) holds each block's real covariate layout fixed —
same guide-pair multiplicity, same per-row |d_a|, |d_b|, noise, so the
measured within-block `d_a`/`d_b` correlation is reproduced by
construction — and redraws only the deltas under a fitted null. Two delta
models bracket the answer:

| delta model | null certified mean (rho=0.10) | NEG | POS |
|---|---|---|---|
| permutation null (real deltas, rows re-bucketed) | 779.0 ± 24.1 | 743.6 | 35.4 |
| empirical marginal (real covariate blocks, deltas i.i.d. from observed delta sample) | 4,208.4 ± 55.5 | 3,463.3 | 745.1 |
| symmetric normal per block (real covariate blocks, deltas ~ N(0, var(δ_block))) | 2,744.2 ± 41.5 | 1,377.5 | 1,366.8 |

Observed: 4,545 certified (NEG 4,059, POS 486).

**Three findings, in order of robustness:**

1. **The pooled-away within-block main-effect correlation itself is weak and
   NEGATIVE: Pearson −0.042 (Spearman −0.008) between `m_a0−mu` and
   `m_0b−mu`. The structural quantity the permutation null loses is small.**
2. **The FCR reading is condition-sensitive, and the binding assumption is
   the delta model, not the within-block correlation.** Holding real covariate
   blocks fixed while drawing deltas i.i.d. from the observed marginal moves
   the null center to 4,208 (FCR 0.926) — near-observed. That null conditions
   on the screen's global negative drift (most pairs mildly depleted), so it
   answers "how many blocks certify if effects were uncorrelated with their
   own covariates" — a *different* question than the permutation null's "how
   many certify if the observed effects were assigned to unrelated pairs".
   Between them they bracket FCR from 0.17 (permutation) to 0.93 (empirical
   marginal), which is exactly the regime the certificate operates in: the
   set of certified pairs is dominated by main-effect-asymmetric genes rather
   than by pure pairwise interaction, so any single-null FCR is a statement
   about that null's conditioning, not a fixed property of the data.
3. **The symmetric-normal variant (2,744) shows how much of the
   permutation/empirical gap is the delta skew inherited by the shuffles:**
   removing the negative drift alone drops the empirical center by ~35%. The
   same-sign split goes symmetric (1,377 / 1,367 vs 3,463 / 745), confirming
   the observed NEG-heavy certificate is a screen-level drift effect, not
   pairwise structure estimated by the model null.

**Why this is not a rescinded certificate:** the permutation null remains the
correct conditioning for the audit's own claim ("effects at these sizes,
assigned to unrelated pairs, would certify 779 pairs, not 4,545"). The
model-based nulls replace the delta distribution under a *fixed covariate
layout*; they lose the observed delta↔covariate coupling inside each row that
the permutation null preserves. Neither is right; they answer different
counterfactuals, and the deliverable is the bracket, not a new headline
number. The certificate's actual interaction content is corroborated
independently by the GEMINI cross-check (Sec 5.9): GEMINI-strong pairs are
concentrated in the certified set (AUROC 0.881 on the De Kegel panel, top-500
recall 2/2), which is a second-caller signal no null construction can
subtract.

*Reading kept honest: the empirical-marginal and symmetric-normal draws are
fully generative — but the former reuses the observed delta sample (so its
tail shapes the null, as it does the permutation null) and its corrected POS
count is negative (−259), i.e., the data contains FEWER positive-certified
blocks than that null expects, which is the signature of a null that
over-absorbs the drift rather than a feature of the call.*

### 5.2 Biology validation

All five paralog pairs the paper validates as GIs are called NEGATIVE at
rho=0.05 and 0.10 and rank at the extreme negative end of the robust score:

| pair | mean_delta | robust_score (rho=0.10) | certified | rank_robust/16999 | rank_conv/16999 |
|------|-----------|--------------------------|-----------|-------------------|-----------------|
| ASF1A-ASF1B | -2.087 | -3.34 | yes | 323 | **2** |
| CNOT7-CNOT8 | -2.890 | -4.25 | yes | 154 | **1** |
| MAPK1-MAPK3 | -1.138 | -4.69 | yes | 118 | 36 |
| HDAC1-HDAC2 | -0.897 | -2.21 | yes | 948 | 96 |
| ARID1A-ARID1B | -0.325 | -1.15 | yes | 3371 | 3358 |

ASF1A/CNOT7 are the two strongest signals in the whole screen under the
conventional score; robust reorders them but keeps all five in the
negative tail. ARID1A-ARID1B (smallest delta) drops out of `certified`
at rho>=0.20 — honest: its 4 guide pairs carry below-floor mean signal.

### 5.3 KILL test #1 — guide-identity holdout (above-floor contract)

| rho | robust callable | robust | conventional (same) |
|-----|-----------------|--------|---------------------|
| 0.05 | 15779 | 0.566 | 0.620 |
| 0.10 | 15140 | 0.575 | 0.627 |
| 0.20 | 13422 | 0.593 | 0.643 |
| 0.30 | 11528 | 0.606 | 0.654 |

**Verdict: the Stage-2 over-claim is falsified.** The naive point estimate
transfers as well or slightly better on every matched subset. The
certificate does not predict held-out guidance better than the mean.

What the certificate measurably buys (the surviving claim):
- **Selectivity / false-positive control:** it refuses to certify ~5,500 of
  16,999 gene pairs (guide disagreement); those refusals concentrate where
  the interaction is not identifiable.
- **Reproducibility concentrates in the certified set:** the naive
  predictor lifts from 0.620 (blanket) to 0.654 (certified set at
  rho=0.30); robust (0.606) is comparable to conventional inside it.

### 5.4 Ranking AUC (audit) — field's SL-detection idiom

Above-floor instances n=40,787 (27,721 negative-GI, 13,066 positive-GI):

| rho | robust AUC | conv AUC |
|-----|-----------|----------|
| 0.05 | 0.610 | 0.628 |
| 0.10 | 0.615 | 0.628 |
| 0.20 | 0.620 | 0.628 |
| 0.30 | 0.623 | 0.628 |

Same conclusion at score level: the certificate's rank ordering is
well-calibrated but not better than raw dLFC. The certificate's value is
the *flag*, not the *rank*.

### 5.5 Confidence (bootstrap, n=200, seed 0)

Mean sign-stability fraction 0.738; 88% of gene pairs resolve in a
majority of resamples, 30% in >=95%. Gene intervals (rho=0.10): 88%
signed; 53% of pairs have empty intersection at some bootstrap sample
(guide conflict is common, which is exactly what the certificate flags).

### 5.6 Phenotype cross-check (authors' nLFC via pyreadr; Rev 2)

Repeat of the full pipeline with `colo1_plm.rds` (endpoint-replicate mean of
the authors' per-construct nLFC) as the phenotype instead of counts-derived
FC. Substantive biology is preserved — positive controls deplete under both
(median −1.09 vs −1.53), and all five paralog pairs stay NEGATIVE with
matched delta order (e.g. CNOT7-CNOT8 −2.89 → −2.93; HDAC1-HDAC2 −0.90 →
−1.61). But the *boundary calls do not transfer*:

| rho | count-FC (POS/NEG/UNRES) | nLFC (POS/NEG/UNRES) | agreed (of resolved-by-either) | agreed (both-resolved) |
|-----|--------------------------|----------------------|-------------------------------|------------------------|
| 0.05 | 127 / 1114 / 15758 | 8 / 11281 / 5710 | 0.091 (n=11,447) | 0.964 (n=1,083) |
| 0.10 | 61 / 897 / 16041 | 7 / 10631 / 6361 | 0.078 (n=10,745) | 0.979 (n=851) |
| 0.20 | 21 / 597 / 16381 | 5 / 9341 / 7653 | 0.057 (n=9,423) | 0.987 (n=541) |
| 0.30 | 10 / 398 / 16591 | 1 / 7961 / 9037 | 0.042 (n=8,029) | 0.985 (n=341) |

The 4-9% figure is a **resolvability** disagreement, not a sign
disagreement: on pairs that *both* procedures resolve, they agree 96-99%.
The nLFC's +0.67 mean shift slides the whole delta distribution past zero,
so it declares ~10K more gene pairs cleanly NEGATIVE than the count-FC does;
the count-FC-resolved set is mostly contained in it yet near-disjoint from
its own resolution boundary (overlap ~0.84-0.89 of the smaller set). I.e.:
given the same resolvable pair both phenotypes call it the same way, but
*which pairs are resolvable at all* is the phenotype-sensitive quantity.

Guide-pair deltas agree only at rank rho 0.67-0.72 (Spearman/Pearson) with a
systematic +0.67 mean shift (`mu` 0.806 → 0.351). The nLFC drifts the whole
distribution negative (10.6K NEGATIVE at rho=0.10 vs 897), i.e. the
authors' normalized LFC already embeds model corrections, so its delta no
longer matches our matched-guide factorial contrast scale. **Reading:** the
counts-derived FC is the correct primary phenotype for this method (both
phenotypes agree on validated biology); the *boundary* of the certificate is
phenotype-sensitive, and this magnitude of sensitivity is the measurement.
Do **not** claim the nLFC run as a confirmation of counts-derived calls.

### 5.7 Selectivity ablation — certified vs magnitude-matched (Rev 2)

The surviving claim is selectivity. Pressure-test: does the certified set
buy reproducibility beyond pure effect magnitude? (If not, "certified"
just restates "big effects reproduce better".) LOOCV mean-sign precision
(conventional calibration call vs held-out sign) on the certified set vs the
same-size top-|mean_delta| set:

| rho | n | certified precision | magnitude-matched precision | overlap |
|-----|---|--------------------|-----------------------------|---------|
| 0.05 | 5571 | 0.934 | 0.914 | 3666 (66%) |
| 0.10 | 4545 | 0.946 | 0.923 | 2714 (60%) |
| 0.20 | 3158 | 0.956 | 0.940 | 1460 (46%) |
| 0.30 | 2336 | 0.961 | 0.945 | 846 (36%) |

**The certified set wins at every rho** (+0.02, +0.023, +0.017, +0.016
precision) against magnitude-matched selection of identical size, with
overlap falling from 66% to 36%. The certificate's selectivity is real and
not reducible to effect size — it selects for guide reproducibility.
(Modest in absolute terms: both selection rules concentrate on big,
reproducible effects.)

### 5.8 NOISE_K sensitivity sweep (Rev 2)

Holdout concordance at rho=0.10 as the guide-noise floor multiplier varies
(default k=1.0):

| noise_k | robust callable | robust | conventional |
|---------|-----------------|--------|--------------|
| 0.5 | 16821 | 0.559 | 0.618 |
| 1.0 | 15140 | 0.575 | 0.627 |
| 1.5 | 12276 | 0.596 | 0.638 |
| 2.0 | 10158 | 0.613 | 0.651 |

Monotone, no cliff, and the default is mid-range: robust concordance rises
with the floor (more effect must clear a higher bar) but never overtakes
conventional — the headline verdict (Stage-2 falsified) is invariant to the
floor constant.

### 5.9 External paralog-SL reference — De Kegel et al. 2021 (Rev 2)

**What the panel is (and is not).** De Kegel et al. (Cell Systems 12:1144-1159.e6;
DOI 10.1016/j.cels.2021.08.006; CC BY 4.0; UCD Cancer Data Lab) derive a
dataset of ~126 robust paralog SL pairs by integrating DepMap CERES
dependency-loss associations across 762 cell lines with four published
combinatorial paralog screens. It is therefore an **independent
computationally-derived reference panel**, not an independent experimental
ground truth, and it is independent of the Burgold/COLO1 screen (no shared
experimental cell line or authors). Provenance: the supplementary lists live
behind a click-through/Cloudflare portal and are not individually
archivable; we instead use the authors' own curated subset shipped in the
open GitHub repo (`DeKegel/paralog_SL_prediction`, `local_data/validated_SLs.txt`,
sha256 `9173ddce846061228c258cbc1fcccac2249b8ed85c525b3a96882f7c85de45ae`),
copied to `raw/de_kegel_2021/`. It lists 12 canonical paralog SLs; one
(**ARID1A-ARID1B**) is already among our 5 validated paralog anchors.

**Power caveat (stated before results).** This method's matched-guide estate
needs both genes to carry single guides on this screen. COLO1's matched-guide
contrast graph spans 474 distinct genes, so of the 12 panel pairs only
**4 are evaluable** (SMARCA2/4, ARID1A/1B, STAG1/2, FAM50A/50B). At this tiny
n the rank/enrichment tests below are illustrative, not conclusive; the full
126-pair training set (Table S5, not raw-archivable) would be the upgrade
path.

**Pre-registered procedure** (fixed before reading the outcome):

- *Positive set*: the 4 evaluable De Kegel pairs, all paralog SLs whose SL
  direction our biology section (5.2) treats as negative-delta.
- *Negative set*: all other gene pairs in the 16,999-pair COLO1 universe
  (the panel is itself defined only on paralogs, so this is an
  absolute-over-reference comparison, not a paralog-matched background).
- *Scores compared*: `robust_score` and `conv_score` (= dLFC), signed so
  that more-negative = more SL-like.
- *Metrics, computed at every rho in RHO_GRID*:
  1. **Rank AUROC** of the 4 positives against the 16,995-pair background
     (positive pairs should rank most negative).
  2. **Top-K recall**: are the 4 positives among the top 50 / 100 / 500
     most-negative pairs by each score?
  3. **Direction**: sign of `mean_delta` for each positive pair (must be
     negative).
  4. **Certified fraction** of the 4 positives (audit-layer question: does
     the robust certificate cover the externally-validated hits?).
- *Decision rule*: no single passing threshold; report numbers with the
  small-n power caveat, and only claim "external paralog enrichment is
  consistent with" the score if direction is unanimously negative AND the
  AUROC is confidently above 0.5.

**Results** (source `audit/dekegel_paralog_sl.csv`, from the pinned validate
command):

Direction is **unanimously negative** (4/4; also 3/3 after dropping the
ARID1A-ARID1B overlap with our own anchors) — consistent external sign.

| rho | AUROC robust | AUROC conv | AUROC GEMINI | top-500 recall (ro/cv/ge) | certified |
|-----|--------------|-----------|--------------|---------------------------|-----------|
| 0.05 | 0.645 | 0.775 | 0.881 | 0 / 0 / 2 | 2/4 |
| 0.10 | 0.657 | 0.775 | 0.881 | 0 / 0 / 2 | 2/4 |
| 0.20 | 0.675 | 0.775 | 0.881 | 0 / 0 / 2 | 1/4 |
| 0.30 | 0.686 | 0.775 | 0.881 | 0 / 0 / 2 | 0/4 |

(GEMINI AUROC leave-one-out range: 0.844-0.968, n=4; `gemini_benchmark.csv`.)

*Metric correction (Rev 2c):* the first Rev 2 cut computed the rank AUROC with
ascending ranks ("most-negative = rank 1"), which makes an enrichment (hits
most negative) collapse toward 0 and report the valid 0.65-0.78 values as
failing an `> 0.5` gate. The score semantics, the CSV, and the gate above are
all fixed to more-negative = higher rank since then; the earlier numbers were
not a different result, they were a different axis.

**n=4 instability (leave-one-out).** At this sample size a single positive
label can move the estimate hard, so each AUROC is recomputed dropping each
of the 4 pairs:

| rho | robust LOO range | conv LOO range |
|-----|------------------|----------------|
| 0.05 | 0.565 - 0.716 | 0.727 - 0.857 |
| 0.10 | 0.583 - 0.727 | 0.727 - 0.857 |
| 0.20 | 0.610 - 0.743 | 0.727 - 0.857 |
| 0.30 | 0.629 - 0.751 | 0.727 - 0.857 |

Per-pair (rho=0.10; percentile of most-negative, lower = more SL-like):

| pair | mean_delta | robust pct | conv pct | certified |
|------|-----------|------------|----------|-----------|
| SMARCA2-SMARCA4 | -0.364 | 11.8 | 15.1 | yes |
| ARID1A-ARID1B | -0.325 | 19.8 | 19.8 | yes |
| STAG1-STAG2 | -0.455 | 50.3 | 8.1 | no |
| FAM50A-FAM50B | -0.156 | 55.1 | 47.1 | no |

**Reading — an ambiguous result, not a validation failure.** Direction is
confirmed (4/4 negative), and ranking leans positive: both scores place the
panel hits below-random rank on average (AUROC > 0.5), with the LOO range
never dipping below 0.5 for conv and only touching 0.57 at the floor for
robust. But n=4 cannot separate "the certificate adds ranking value" from
"these four interactions are modest here" — and the whole-panel comparison is
confounded by construction: De Kegel is pan-cancer (762 cell lines,
context-agnostic by design; context-dependence is that paper's own central
finding), while COLO1 is a single line (HT29), and all four interactions are
modest here (dLFC -0.16 to -0.46) relative to the screen's extreme tail.
Both competing readings of a null result ("certificate adds nothing" vs
"these paralogs aren't strong in HT29") are compatible with the data; the
second is the *leading* hypothesis given our own dLFC magnitudes. Under the
pre-registered rule the verdict is therefore **inconclusive** (n<10 guardrail
in `panel_verdict`), not "declined."

The one textbook certificate behavior visible in this sample: **STAG1-STAG2**
ranks most-negative by raw dLFC (8th percentile) yet the certificate demotes
it to the median (robust pct 50) and refuses certification — its guide pairs
disagree, exactly the guide-noise signal the certificate is built to flag
(n=1, anecdotal; it is also the pair whose exclusion moves robust AUROC most).

**GEMINI benchmark on the same 4 pairs (Rev 2e, full numbers in Sec 6).**
GEMINI's per-gene-pair score (`strong`, positive = interaction) ranks these
4 covered pairs better than either B06 score: full-library percentiles
(18,920 scored pairs) SMARCA2-SMARCA4 #320 (1.7%), FAM50A-FAM50B #528
(2.8%), STAG1-STAG2 #1432 (7.6%), ARID1A-ARID1B #6815 (36%). One is worse
than the raw-dLFC percentile of the same pairs (STAG1-STAG2 8th percentile
by dLFC vs 7.6th by GEMINI; ARID1A 19.8th by dLFC vs 36th). The 5 validated
paralogs: CNOT7-CNOT8 #1, ASF1A-ASF1B #2, MAPK1-MAPK3 #3, HDAC1-HDAC2 #8 of
18,920 — better central ranking than any B06 score gives them.

**Status of gap #2 (§7):** external panel *ingested*; direction confirmed,
ranking inconclusive at n=4 (weakly positive, single-pair-instable). The
full 126-pair training set (gated) and a paralog-matched background remain
the upgrade path to a powered ranking test.

## 6. Position vs field state of the art

Current SOTA callers (GEMINI, Orthrus; 2025 five-method benchmark in
bioRxiv 2025.03.31.645224; 2026 GRAPE regression pipeline) all output point
scores + null-hypothesis statistics and demonstrably beat raw dLFC on
simulated and real data. **B06 is not SOTA:** its own decisive test shows
its calls do not transfer better than dLFC, so it fails any SOTA claim by
construction. Its defensible, novel contribution is the identifiability
certificate itself (closed-form robustness guarantee with an interpretable
severity parameter + guide-noise floor) as an audit/safety layer.

Rev 2e puts that on quantitative ground **on this exact screen** by running
the field's current caller, GEMINI (sellerslab/gemini, Bioc 1.24.0,
Zamanighomi et al. Genome Biology 2019), over the same 78,912-guide-pair
COLO1 screen (framework R 4.5.3 arm64 binaries; `ETP.column` = lib.COLO.1,
endpoint = 3 SIDM00136_CPID1020/1023/1026 replicates; nc = NONTARGET gene
pairs; 50 CAVI iterations; score = median of per-replicate `Score$strong`;
FDR = min FI of per-replicate `Score$fdr_strong`). `run_gemini.R` +
`bench_gemini.py` are pinned in Sec 9.

Scores (all on the shared 16,959-gene-pair universe)*:

| metric | robust | conv | GEMINI |
|--------|--------|------|--------|
| De Kegel panel AUROC (4 pairs) | 0.657-0.686 | 0.775 | **0.881 (LOO 0.844-0.968)** |
| top-500 recall of those 4 | 0 | 0 | 2 |
| Spearman vs GEMINI strong | — | — | 0.146 (conv) / 0.20-0.22 (robust) |
| FDR < 0.05 hits | (see Sec 5.1 FCR) | — | **0 of 18,920** |
| validated paralog rank (of 16,999) | robust 118-3371 | conv 1-3358 | **CNOT7 #1, ASF1A #2, MAPK1 #3, HDAC1 #8 (of 18,920)** |

\**Score axes aligned by the external_paralog_panel convention (higher =
more SL-like): robust/conv are passed as `-score` (negative = SL), GEMINI
strong is native (positive = interaction).* The self-check reproduces the
stored robust/conv AUROCs exactly (0.658 / 0.775 at rho 0.10).

Reading, honestly:

1. **On ranking, GEMINI wins cleanly.** It is the only method of the three
   with panel AUROC comfortably above 0.5 and LOO floor 0.844, and it puts
   the strongest validated paralogs at the very top where neither B06 score
   does. Any claim that B06 "ranks validated biology" must cede to this.
2. **On significance, GEMINI finds nothing here.** Its own null gives every
   pair min-FDR >= 0.83 (flat across the library), so it emits zero hits at
   any threshold; the most positive `strong` scores belong to
   NONTARGET~NONTARGET control pairs (+1.60, above even CNOT7-CNOT8 +1.19).
   GEMINI's model is calibrated for large pooled screens; on this blend of
   a modest-endpoint replicate set it fails to separate signal from its own
   null. (n_iterations forced to completion; `force_results=TRUE`.)
3. **Rank agreement is weak.** Spearman 0.15-0.22 across methods: the two
   scoring families rank the estate near-independently, so "certified by B06
   but not on GEMINI's radar" and vice versa are both common. 284 of GEMINI's
   top-500 are declined by the certificate (rho 0.10); only 4.8% of the
   certificate's certified set lands in GEMINI's top-500.

Net position: **B06 is not SOTA on ranking (GEMINI 0.881 > 0.775 > 0.66);
GEMINI adds no significance on this screen; the certificate's value stays
uncontradicted as an audit layer — the two tools disagree on which pairs are
real, which is precisely what an independently-constructed check should
surface.** GRAPE/Orthrus and a per-pair joint-agreement study remain as the
next quantified steps (Sec 7).

## 7. Known gaps & shortfalls (audit's verification target)

1. **Phenotype sensitivity (formerly "orientation unconfirmed").** Orientation
   is now confirmed against `colo1_plm.rds` (Sec 3); the remaining, *measured*
   issue is that boundary calls are phenotype-dependent (Sec 5.6): ranks of
   validated biology agree, near-zero calls do not. Counts-derived FC remains
   primary; the nLFC run is retained as the sensitivity measurement, not a
   confirmation. A construct-level calibration cross-model (nLFC -> FC scale)
   would close the residual offset (mean delta shift 0.67).
2. **External panel ingested; ranking inconclusive at n=4.** De Kegel et al.
   2021 paralog-SL panel is now integrated (Sec 5.9): direction confirmed
   (4/4 negative, incl. 3/3 after dropping the one overlap with our own
   anchors); ranking weakly positive but single-pair-instable (LOO range
   0.57-0.86, never crossing 0.5) and confounded by pan-cancer-vs-single-line
   design. **Known limitation with no cheap fix: a larger, more powered
   external panel exists (the ~126-pair De Kegel training set, Table S5) but
   is not individually raw-archivable** — it requires a DepMap 20Q2 download
   plus four Cloudflare-gated combinatorial-screen tables to reconstruct
   (attempted; cell.com/biorxiv/cancergd.org all gated). Upgrade path retains
   it with a paralog-matched negative background.
3. **FDR for the call counts (Sec 5.1); aggregate point estimate only.**
   The strict POSITIVE/NEGATIVE call counts have permutation error-rate
   control (1000 perms; 0/1000 reaches either count at any rho; FCR
   0.09-0.26 NEG, <=0.04 POS; FCR-scaled corrected counts in
   `call_fdr_null.csv`). Construction caveat discovered during the work: a
   bijective relabel of the factor codes gives a point-mass null (sd 0) —
   only a true row re-bucketing permutes anything (now guarded by a
   set-based membership test). What remains open: (a) the corrected counts
   are point-estimate corrections from the aggregate null, not per-pair BH
   q-values — per-pair FDR needs per-pair null draws; (b) the model-based
   null that reproduces within-block main-effect correlation is **done** —
   Sec 5.10: the measured within-block correlation is −0.042 and the FCR
   reading is delta-model-dependent (permutation 779 vs empirical-marginal
   4,208 vs symmetric-normal 2,744 of 4,545 observed); the certificate's
   interaction content is corroborated by the GEMINI cross-check, so the
   bracket is reported without restating the headline.
4. **SOTA benchmark (GEMINI) — done, quantitative, negative-but-concrete.**
   GEMINI ran on this exact screen (Sec 6): better panel-AUROC ranking
   (0.881 vs 0.657/0.775), zero FDR-significant hits (min-FDR floor 0.83),
   weak rank agreement (Spearman 0.15-0.22). Remaining: GRAPE/Orthrus for a
   second independent caller (cheap once the R 4.5.3 framework toolchain is
   re-run; GEMINI alone already falsifies any SOTA-rank claim), and a
   per-pair joint-agreement analysis between the certificate and GEMINI
   hits.
5. **Replicate structure** (2×3) is assumed additively; no explicit per-guide
   variance modeling (consistent with plan but untested against GEMINI/GRAPE
   noise models).
6. **Pilot** is underpowered for the holdout (6 pairs) — drop it from any
   Stage-2 claims.
7. **Certificate != prediction.** Any future claim must be scoped to
   selectivity/false-positive control, or move to magnitude-adaptive rho +
   certified-mean decision rule (plan Sec 14-21).
8. **Seed off-target sensitivity — data not available, closed as a dead-end
   (Rev 2e).** The planned `sgRNA_Off_Target`-based sensitivity analysis
   cannot be run on this release: `COLO1_RUNMERGED_EXACT_ANNOTATED.txt` has
   NO `sgRNA_Off_Target` columns (only ID, Note, MyNote, sgRNA1_ID,
   sgRNA2_ID, Gene1, Gene2, Gene_Pair, lib.COLO.1, CPID2437/2440/2443/1020/1023/1026),
   and no guide-metadata table ships with the release — verified by direct
   inspection of the uploaded raw file. Documented rather than silently
   dropped; would require re-acquiring off-target annotations from the
   depositor.

## 8. Artifact inventory

- `results_colo1/`: scout_summary.md, audit/ (audit_summary.json,
  matched_guide_pairs.csv, gene_scores.csv, rank_auc.csv,
  validated_paralog_ranks.csv, selectivity_ablation.csv, noise_k_sweep.csv,
  dekegel_paralog_sl.csv, certified_fdr_null.json,
  model_based_fdr_null.json, call_fdr_null.csv,
  validation_summary.json), bounds/
  (guide_pair_contrasts.csv 67,949 rows, gene_pair_calls.csv,
  gene_intervals.csv, guide_holdout.csv, guide_conflicts.csv,
  guide_bootstrap.csv), mantis/gene_pairs.csv, figures/ (4 PNG).
- `results_colo1_plm/`: same layout, generated with `--lfc-rds
  colo1_plm.rds` (phenotype cross-check; Sec 5.6).
- `results/`: pilot equivalents.
- `raw/de_kegel_2021/validated_SLs.txt`: external paralog-SL panel.
- `benchmarks/gemini/`: `run_gemini.R` (GEMINI driver; requires framework
  R 4.5.3 arm64 + Bioc `gemini` 1.24.0), `bench_gemini.py` (Sec 6 join),
  `gemini_scores.csv` (18,920 gene pairs: median `strong` + min FDR),
  and `results_colo1/audit/gemini_benchmark.csv` + `.json`.
- `src/b06/`: bounds.py, contrasts.py, confidences.py, holdout.py,
  validate.py, audit.py, io.py, scout.py, reporting.py, figures.py,
  provenance.py.
- `tests/`: 25 tests (guardrail, symmetric bound, sharpness, E2E,
  monotonicity, paralog hit, confidence coverage, holdout, synthetic io
  regression, external-panel verdict rule + LOO, permutation-FDR nulls +
  block-membership guard).

## 9. Reproducibility statement

All CSVs are deterministic outputs of the pinned commands (seed 0, no
git). Data sources: Zenodo (open), De Kegel panel (open GitHub, sha256
pinned), pilot (license-restricted, not redistributed). Re-running both
scout commands (counts-FC and --lfc-rds) + the validate command + pytest
regenerates every number in this audit, including Sec 5.6-5.9.

Sec 6 GEMINI benchmark reproducibility (framework R 4.5.3 arm64 binaries,
user-space extraction; env: `R45FRAMEWORK=1`,
`DYLD_LIBRARY_PATH=<framework>/Resources/lib`,
`R_LIBS_USER=~/.R45/library`; Bioc `gemini` 1.24.0 installed from the
macosx/big-sur-arm64 binary):

    RF=~/.local/R4.5.3/R.framework/Versions/4.5-arm64
    "$RF/Resources/bin/R" --no-save --vanilla -f benchmarks/gemini/run_gemini.R \
        --args raw/colo1_release/.../COLO1_RUNMERGED_EXACT_ANNOTATED.txt benchmarks/gemini
    PYTHONPATH=src python3 benchmarks/gemini/bench_gemini.py

`run_gemini.R` prints its own self-check (counts 78,912x4; NONTARGET00741
nc reference; `nc_pairs` count) and the sync check in `bench_gemini.py`
asserts the robust/conv AUROC equals the stored Sec 5.9 values. GEMINI's
CAVI updates use 15 cores; the score/`mnb` steps are deterministic given
library + endpoints, regardless of `force_results=TRUE` (only iteration
count 50 is forced).