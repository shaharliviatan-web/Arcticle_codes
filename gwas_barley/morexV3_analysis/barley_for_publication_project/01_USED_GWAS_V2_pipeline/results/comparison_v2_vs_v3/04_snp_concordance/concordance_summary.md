# SNP concordance: v2 flags vs v3

Generated 2026-08-16. Headline config: **BLUP x 3 PCs**.

## Coordinates are identical by construction

Both runs used the same 7,110,996-SNP tped; only the covariates changed. Verified against `snp_map.bim`: **coord_match TRUE for 24 / 24** flagged SNPs. No SNP moved position, and none is missing from either run.

So the question is purely whether each SNP keeps its **class**:

- **significant**: -log10p >= 6.0454 (v3), 6.7712 (v2)
- **below_threshold**: everything else

**No `marginal` class is assigned to v3.** In v2 that set was curated by hand, not by a rule, so the v3 marginals have to be chosen manually -- the raw `v3_nlp` is given here so that selection can be made.

`v3_best` = strongest v3 signal for that SNP across all 6 configs, so a SNP that fades at 3 PCs but survives elsewhere is still visible.

## Forward: every SNP v2 flagged, looked up in v3

### Betaglucan

| SNP | bp | v2 -log10p | v2 class | v3 -log10p | v3 class | delta | v3 best | best config |
|---|---|---|---|---|---|---|---|---|
| 2H:41996626 |  41996626 | 8.9060 | significant | 7.6023 | significant | -1.304 | 7.6023 | BLUP x 3 PC |
| 4H:34858354 |  34858354 | 7.9777 | significant | 7.3349 | significant | -0.643 | 7.3912 | BLUE x 3 PC |
| 3H:537247620 | 537247620 | 7.9751 | significant | 7.1876 | significant | -0.788 | 7.5490 | BLUE x 5 PC |
| 4H:34971609 |  34971609 | 7.7266 | significant | 7.0864 | significant | -0.640 | 7.1239 | BLUE x 3 PC |
| 4H:34861659 |  34861659 | 7.2971 | significant | 6.9020 | significant | -0.395 | 6.9438 | BLUE x 3 PC |
| 6H:475144663 | 475144663 | 7.1682 | significant | 6.2708 | significant | -0.897 | 6.2894 | BLUE x 3 PC |
| 4H:34825152 |  34825152 | 7.1073 | significant | 6.7759 | significant | -0.331 | 6.8122 | BLUE x 3 PC |
| 2H:602246902 | 602246902 | 7.0568 | significant | 5.8087 | below_threshold | -1.248 | 5.8630 | BLUE x 3 PC |
| 1H:479493867 | 479493867 | 6.9044 | significant | 6.2281 | significant | -0.676 | 6.2978 | BLUE x 3 PC |
| 4H:528254769 | 528254769 | 6.8965 | significant | 6.5307 | significant | -0.366 | 6.9308 | BLUE x 10 PC |
| 1H:483379118 | 483379118 | 6.8860 | significant | 5.5181 | below_threshold | -1.368 | 5.5425 | BLUE x 3 PC |
| 3H:537301368 | 537301368 | 6.7886 | significant | 6.2963 | significant | -0.492 | 6.5316 | BLUE x 5 PC |
| 7H:144534576 | 144534576 | 6.7694 | marginal | 5.7404 | below_threshold | -1.029 | 5.7566 | BLUE x 3 PC |
| 5H:224439424 | 224439424 | 6.6805 | marginal | 6.3553 | significant | -0.325 | 6.6028 | BLUE x 5 PC |
| 6H:545947347 | 545947347 | 6.6540 | marginal | 6.1146 | significant | -0.539 | 6.1146 | BLUP x 3 PC |

### Fiber

| SNP | bp | v2 -log10p | v2 class | v3 -log10p | v3 class | delta | v3 best | best config |
|---|---|---|---|---|---|---|---|---|
| 7H:573606306 | 573606306 | 7.1456 | significant | 7.2451 | significant | 0.099 | 7.6180 | BLUP x 10 PC |
| 4H:24707376 |  24707376 | 6.9727 | significant | 7.0209 | significant | 0.048 | 8.6560 | BLUE x 5 PC |
| 7H:14817657 |  14817657 | 6.7307 | marginal | 6.7984 | significant | 0.068 | 7.5240 | BLUE x 5 PC |
| 1H:344520079 | 344520079 | 6.5484 | marginal | 6.8521 | significant | 0.304 | 6.9568 | BLUE x 5 PC |

### Protein

| SNP | bp | v2 -log10p | v2 class | v3 -log10p | v3 class | delta | v3 best | best config |
|---|---|---|---|---|---|---|---|---|
| 3H:106623911 | 106623911 | 6.84 | significant | 6.433 | significant | -0.407 | 6.6259 | BLUP x 5 PC |

### Starch

| SNP | bp | v2 -log10p | v2 class | v3 -log10p | v3 class | delta | v3 best | best config |
|---|---|---|---|---|---|---|---|---|
| 6H:525776080 | 525776080 | 6.4439 | marginal | 6.1175 | significant | -0.326 | 6.1385 | BLUE x 3 PC |
| 7H:151110354 | 151110354 | 6.3624 | marginal | 5.9304 | below_threshold | -0.432 | 5.9678 | BLUE x 3 PC |
| 3H:546433616 | 546433616 | 6.2818 | marginal | 6.3441 | significant |  0.062 | 6.4405 | BLUE x 3 PC |
| 7H:573606460 | 573606460 | 6.1949 | marginal | 6.4367 | significant |  0.242 | 6.9719 | BLUP x 10 PC |

## Class transitions (forward)

| v2 class | v3 class | n SNPs |
|---|---|---|
| marginal | significant |  7 |
| marginal | below_threshold |  2 |
| significant | significant | 13 |
| significant | below_threshold |  2 |

Significant in v3 (headline config): **20 / 24**.
Significant in at least one of the 6 v3 configs: **20 / 24**.

## Reverse: every v3-significant lead SNP, looked up in v2

LD-clumped at +/- 188 kb (the window v2 used).

### Betaglucan (11 lead SNPs)

| SNP | bp | v3 -log10p | v2 -log10p | v2 class | delta | in v2 tables |
|---|---|---|---|---|---|---|
| 2H:41996626 |  41996626 | 7.6023 | 8.9060 | significant | -1.304 |  TRUE |
| 4H:34858354 |  34858354 | 7.3349 | 7.9777 | significant | -0.643 |  TRUE |
| 3H:537247620 | 537247620 | 7.1876 | 7.9751 | significant | -0.787 |  TRUE |
| 4H:528254769 | 528254769 | 6.5307 | 6.8965 | significant | -0.366 |  TRUE |
| 5H:224439424 | 224439424 | 6.3553 | 6.6805 | below_threshold | -0.325 |  TRUE |
| 6H:475144663 | 475144663 | 6.2708 | 7.1682 | significant | -0.897 |  TRUE |
| 1H:479493867 | 479493867 | 6.2281 | 6.9044 | significant | -0.676 |  TRUE |
| 6H:477651616 | 477651616 | 6.1919 | 6.2945 | below_threshold | -0.103 | FALSE |
| 4H:35062236 |  35062236 | 6.1768 | 6.5171 | below_threshold | -0.340 | FALSE |
| 6H:545947347 | 545947347 | 6.1146 | 6.6540 | below_threshold | -0.539 |  TRUE |
| 3H:475967335 | 475967335 | 6.0576 | 5.9526 | below_threshold |  0.105 | FALSE |

### Fiber (19 lead SNPs)

| SNP | bp | v3 -log10p | v2 -log10p | v2 class | delta | in v2 tables |
|---|---|---|---|---|---|---|
| 7H:573606306 | 573606306 | 7.2451 | 7.1456 | significant | 0.099 |  TRUE |
| 4H:24707376 |  24707376 | 7.0209 | 6.9727 | significant | 0.048 |  TRUE |
| 1H:344520079 | 344520079 | 6.8521 | 6.5484 | below_threshold | 0.304 |  TRUE |
| 7H:14817657 |  14817657 | 6.7984 | 6.7307 | below_threshold | 0.068 |  TRUE |
| 7H:586569512 | 586569512 | 6.6754 | 6.0722 | below_threshold | 0.603 | FALSE |
| 7H:415651802 | 415651802 | 6.4915 | 4.9711 | below_threshold | 1.520 | FALSE |
| 1H:330239605 | 330239605 | 6.4843 | 5.2193 | below_threshold | 1.265 | FALSE |
| 5H:437005013 | 437005013 | 6.4594 | 5.7084 | below_threshold | 0.751 | FALSE |
| 5H:377309369 | 377309369 | 6.3720 | 5.4120 | below_threshold | 0.960 | FALSE |
| 3H:543690146 | 543690146 | 6.3664 | 5.8271 | below_threshold | 0.539 | FALSE |
| 1H:330625182 | 330625182 | 6.3071 | 5.6511 | below_threshold | 0.656 | FALSE |
| 3H:544084309 | 544084309 | 6.2422 | 5.7766 | below_threshold | 0.466 | FALSE |
| 1H:440573355 | 440573355 | 6.2189 | 4.7991 | below_threshold | 1.420 | FALSE |
| 1H:329707484 | 329707484 | 6.1816 | 5.2125 | below_threshold | 0.969 | FALSE |
| 1H:440022799 | 440022799 | 6.1584 | 4.7281 | below_threshold | 1.430 | FALSE |
| 3H:545255371 | 545255371 | 6.1351 | 5.3760 | below_threshold | 0.759 | FALSE |
| 5H:208944021 | 208944021 | 6.0799 | 5.4179 | below_threshold | 0.662 | FALSE |
| 4H:470989624 | 470989624 | 6.0644 | 5.7484 | below_threshold | 0.316 | FALSE |
| 5H:462564771 | 462564771 | 6.0506 | 4.6864 | below_threshold | 1.364 | FALSE |

### Protein (2 lead SNPs)

| SNP | bp | v3 -log10p | v2 -log10p | v2 class | delta | in v2 tables |
|---|---|---|---|---|---|---|
| 3H:106623911 | 106623911 | 6.4330 | 6.8400 | significant | -0.407 |  TRUE |
| 3H:159799839 | 159799839 | 6.0456 | 5.6464 | below_threshold |  0.399 | FALSE |

### Starch (5 lead SNPs)

| SNP | bp | v3 -log10p | v2 -log10p | v2 class | delta | in v2 tables |
|---|---|---|---|---|---|---|
| 7H:573606460 | 573606460 | 6.4367 | 6.1949 | below_threshold |  0.242 |  TRUE |
| 3H:546433616 | 546433616 | 6.3441 | 6.2818 | below_threshold |  0.062 |  TRUE |
| 3H:530660057 | 530660057 | 6.2104 | 5.8744 | below_threshold |  0.336 | FALSE |
| 1H:330432546 | 330432546 | 6.1296 | 5.5839 | below_threshold |  0.546 | FALSE |
| 6H:525776080 | 525776080 | 6.1175 | 6.4439 | below_threshold | -0.326 |  TRUE |

## Reverse summary

| v2 class of v3-significant SNPs | n |
|---|---|
| below_threshold | 28 |
| significant |  9 |

New in v3 (not in the v2 significant or marginal tables): **21 / 37**.

