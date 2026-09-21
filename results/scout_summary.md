# Stage 0/1 scout summary (Burgold pilot library)

- guide pairs: 8869
- unique guides: 3146
- genes: 538 (2+ guides: 538)
- gene pairs observed: 1249 (full cross-product: 521)
- NTC control pairs: 595, control phenotype median: 0.949
- guides with a matched single control: 1071
- guides bridging single and double contexts: 214

- matched guide pairs computable: 163
- guide-pair contrasts computable: 163

### Sign calls vs rho (gene-level, resolvable pairs only)

| rho | gene pairs with calls | positive | negative | unresolved | guide conflicts |
|-----|----------------------|----------|----------|------------|-----------------|
| 0.05 | 44 | 12 | 3 | 29 | |
| 0.10 | 44 | 8 | 1 | 35 | |
| 0.20 | 44 | 4 | 1 | 39 | |
| 0.30 | 44 | 2 | 1 | 41 | |

### Guide-identity holdout (LOOCV by guide pair)

| rho | gene pairs | robust callable | robust concordance | conventional (same pairs) |
|-----|------------|-----------------|--------------------|---------------------------|
| 0.05 | 6 | 6 | 0.692 | 0.787 |
| 0.10 | 6 | 6 | 0.692 | 0.787 |
| 0.20 | 6 | 6 | 0.769 | 0.787 |
| 0.30 | 6 | 6 | 0.614 | 0.787 |
