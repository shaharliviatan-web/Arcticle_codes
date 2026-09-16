# v3 vs v2 — hits and top signal, all traits x all configs

Generated 2026-08-16 from `v2_vs_v3/hit_counts_compare.tsv`.

**How to read this.** Each version is first scored at its own Bonferroni line
(v3 = 6.0454 on 111,017 pruned SNPs; v2 = 6.7712 on 590,462). The last column
rescores the **v2** run at the **v3** threshold — compare it against `v3 hits`
to see what the new correction did, with the threshold change held constant.

## Betaglucan

| config | v3 top -log10p | v3 hits | v2 top -log10p | v2 hits | v2 hits @6.045 |
|---|---|---|---|---|---|
| BLUP x 3 PC (headline) | 7.602 | 20 | 8.906 | 12 | 51 |
| BLUP x 5 PC | 7.515 | 22 | 8.735 | 11 | 44 |
| BLUP x 10 PC | 7.435 | 15 | 8.340 | 14 | 49 |
| BLUE x 3 PC | 7.494 | 20 | 8.784 | 13 | 54 |
| BLUE x 5 PC | 7.549 | 22 | 8.618 | 12 | 47 |
| BLUE x 10 PC | 7.506 | 15 | 8.226 | 15 | 51 |

Strongest v3 signal for betaglucan: **BLUP x 3 PC** at -log10p = **7.602** (20 hits).

## Fiber

| config | v3 top -log10p | v3 hits | v2 top -log10p | v2 hits | v2 hits @6.045 |
|---|---|---|---|---|---|
| BLUP x 3 PC (headline) | 7.245 | 24 | 7.146 | 2 | 5 |
| BLUP x 5 PC | 8.323 | 34 | 7.519 | 3 | 5 |
| BLUP x 10 PC | 8.021 | 8 | 7.601 | 2 | 5 |
| BLUE x 3 PC | 7.294 | 33 | 7.236 | 3 | 6 |
| BLUE x 5 PC | 8.656 | 40 | 7.711 | 3 | 6 |
| BLUE x 10 PC | 8.340 | 12 | 7.896 | 3 | 5 |

Strongest v3 signal for fiber: **BLUE x 5 PC** at -log10p = **8.656** (40 hits).

## Protein

| config | v3 top -log10p | v3 hits | v2 top -log10p | v2 hits | v2 hits @6.045 |
|---|---|---|---|---|---|
| BLUP x 3 PC (headline) | 6.433 | 2 | 6.840 | 1 | 2 |
| BLUP x 5 PC | 6.626 | 2 | 6.859 | 1 | 2 |
| BLUP x 10 PC | 6.565 | 1 | 6.971 | 1 | 2 |
| BLUE x 3 PC | 6.367 | 1 | 6.778 | 1 | 2 |
| BLUE x 5 PC | 6.558 | 2 | 6.796 | 1 | 2 |
| BLUE x 10 PC | 6.489 | 1 | 6.899 | 1 | 2 |

Strongest v3 signal for protein: **BLUP x 5 PC** at -log10p = **6.626** (2 hits).

## Starch

| config | v3 top -log10p | v3 hits | v2 top -log10p | v2 hits | v2 hits @6.045 |
|---|---|---|---|---|---|
| BLUP x 3 PC (headline) | 6.437 | 6 | 6.444 | 0 | 4 |
| BLUP x 5 PC | 6.879 | 4 | 6.691 | 0 | 3 |
| BLUP x 10 PC | 6.972 | 5 | 7.006 | 2 | 3 |
| BLUE x 3 PC | 6.441 | 6 | 6.466 | 0 | 4 |
| BLUE x 5 PC | 6.859 | 4 | 6.682 | 0 | 3 |
| BLUE x 10 PC | 6.955 | 4 | 6.999 | 2 | 3 |

Strongest v3 signal for starch: **BLUP x 10 PC** at -log10p = **6.972** (5 hits).

## Totals across all 24 cells

| scoring | hits |
|---|---|
| v3 at its own threshold (6.0454) | 303 |
| v2 at its own threshold (6.7712) | 103 |
| v2 rescored at the v3 threshold (6.0454) | 360 |

Note: v2 rescored at the v3 line yields **360** hits versus v3's **303**. The headline jump from 103 to 303 is therefore driven by the easier threshold, not by the new correction. See the README for why the correction is still the better-calibrated model (lambda).

