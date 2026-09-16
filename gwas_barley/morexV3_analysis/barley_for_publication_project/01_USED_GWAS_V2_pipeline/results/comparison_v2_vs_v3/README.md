# v2 vs v3 population-structure correction — comparison set

Created 2026-08-16.

## What changed

The PCA + aIBS kinship used to correct for population structure is built on an
LD-pruned marker set. Only the **pruning window** changed between the two runs.

| | v2 | v3 |
|---|---|---|
| PLINK call | `--indep-pairwise 50 5 0.2` | `--indep-pairwise 1000kb 1 0.2` |
| Window | 50 SNPs | 1,000 kb (physical) |
| Markers kept | 590,462 (8.3%) | **111,017 (1.56%)** |
| Bonferroni α=0.10 | 6.7712 | **6.0454** |

Rationale: barley SNP density here averages 1,701 SNPs/Mb (max 7,825), so a fixed
50-SNP window spans a wildly variable physical distance. It removed only very local
LD and left long-range LD blocks over-represented — inflating both the leading PCs
and the kinship off-diagonal, and making the "independent tests" denominator larger
than it should be.

Everything else is identical: same 290 samples, same full 7,110,996-SNP tped, same
24-cell grid (4 traits × {BLUP, BLUE} × {3, 5, 10} PCs), same EMMAX version.
**Any difference in these figures is attributable to the correction alone.**

## Directory layout

Each topic has a `v3_only/` view (read the new run on its own terms) and a
`v2_vs_v3/` view (direct comparison).

```
01_pc_selection/
  v3_only/    scree_plot, pca_scatter_PC1_PC2 / PC1_PC3 / PC2_PC3,
              lambda_vs_npcs, pc_variance_table.tsv, lambda_vs_npcs_{long,wide}.tsv
  v2_vs_v3/   scree_compare, scree_cumulative_compare, scree_elbow_compare,
              pca_scatter_compare, pc_score_correlation, pc_crosscorrelation,
              lambda_vs_npcs_compare, pc_variance_compare.tsv,
              pc_score_correlation.tsv, pc_best_match.tsv
02_qq_lambda/
  v3_only/    qq_grid_BLUP, qq_grid_BLUE, lambda_v3, lambda_v3.tsv
  v2_vs_v3/   qq_overlay_BLUP, qq_overlay_BLUE, lambda_shift, lambda_compare.tsv
03_manhattan/
  v2_vs_v3/   mirror_<pheno>_pc<N> (v3 up / v2 down), hit_counts_compare.tsv
```

Every figure is written as both `.pdf` (vector, TAG spec) and `.png` (300 dpi).

## Colour convention (consistent across every figure)

- **grey** = v2, 50-SNP window
- **blue** = v3, 1000 kb window
- **red** = reference line (y = x on QQ, λ = 1, Bonferroni threshold)

## Headline findings

**Genomic inflation improved.** 21 of 24 cells moved closer to λ = 1; mean |λ−1|
halved, 0.0227 → 0.0136. v2 sat uniformly *below* 1 (0.969–0.987) — systematic
over-correction from an inflated kinship matrix. v3 spans 0.974–1.013.

**Kinship is less inflated.** Mean off-diagonal 0.742 → 0.671, diagonal 0.952 → 0.939.

**The population structure itself is unchanged.** PC1–PC3 sample scores correlate
|r| = 0.99 / 0.98 / 0.98 between versions. What dropped is the *share of variance*
those PCs absorb (cumulative at PC10: 33.3% → 23.4%) — redundant LD blocks removed,
not structure lost. The total eigenvalue trace is essentially identical
(581.2 → 580.4).

**PC4+ appear uncorrelated only because they reorder.** Rank-matched correlations
look poor from PC4 on, but the cross-correlation matrix shows v3_PC4 ↔ v2_PC5
(|r| = 0.80) and v3_PC9 ↔ v2_PC10 (|r| = 0.68). PCs with near-equal eigenvalues swap
rank between marker sets; this is expected and is not instability.
See `01_pc_selection/v2_vs_v3/pc_crosscorrelation.png` and `pc_best_match.tsv`.

**Hit counts changed in both directions.** Scoring v2 at the *same* 6.0454 threshold
isolates the correction's effect from the easier threshold (BLUP × 3 PCs):

| Trait | v3 @6.045 | v2 @6.771 | v2 @6.045 |
|---|---|---|---|
| betaglucan | 20 | 12 | 51 |
| fiber | 24 | 2 | 5 |
| protein | 2 | 1 | 2 |
| starch | 6 | 0 | 4 |

Fiber gained substantially, betaglucan lost, protein/starch roughly unchanged. The
QQ overlays explain why: v2's betaglucan tail was inflated above the null, while
fiber's genuine polygenic signal was being flattened by over-correction.

Three of four lead SNPs are **identical** to v2 (betaglucan 2H:41996626,
fiber 7H:573606306, protein 3H:106623911). Starch's lead moved to 7H:573606460 —
154 bp from fiber's lead, the same locus.

## Caveat to carry into the methods

Fiber's v3 QQ curve departs from the null relatively early (expected ≈ 1.5–2).
With λ = 0.975 it is not globally inflated, so this reads as polygenic architecture
rather than residual structure — and the lift shrinks as PCs increase. Worth a
sentence in the methods if fiber becomes a headline result.

## Reproducing

| Script | Produces |
|---|---|
| `scripts/08c_pc_selection_diagnostics.R` | `01_pc_selection/v3_only/` (via `results/pc_selection/`) |
| `scripts/08d_pc_selection_compare_v2_v3.R` | `01_pc_selection/v2_vs_v3/` |
| `scripts/06c_qq_lambda_v3_grid.R` | `02_qq_lambda/v3_only/` |
| `scripts/06b_qq_compare_v2_v3.R` | `02_qq_lambda/v2_vs_v3/` |
| `scripts/07b_manhattan_compare_v2_v3.R` | `03_manhattan/v2_vs_v3/` |

The v2 reference data lives in `results/_archive_v2_win50snp_2026-08-16/` and
`intermediates/_archive_v2_win50snp_2026-08-16/`. **Both are read-only archives —
nothing in this comparison writes into them.**
