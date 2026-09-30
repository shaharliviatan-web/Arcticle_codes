# 01_block4_2019_interblock_correlations

Inter-block Pearson correlations of the four NIR grain traits (protein, starch,
β-glucan, fiber; grain content, %) in the wild-barley common garden, per season.

Run on 2026-09-29 with R 4.1.2 (2021-11-01).

**Start here → [`results/tables/results_numbers.txt`](results/tables/results_numbers.txt)**:
the key numbers in plain prose, ready to quote.

## Purpose

To document why block 4 of the 2019–20 season was left out of the phenotype
analysis (`00_THIN_Generate_Plots_For_Publication/01_phenotypic_analysis_no_GxE_v3_streamlined.R`,
section 2, "Filter 2: drop Block 4 of season 2019 (low inter-block concordance)").
These results support one sentence and possibly one Online Resource in the
Materials and methods of the TAG manuscript.

## Input (read-only)

`00_THIN_Generate_Plots_For_Publication/all years barley.csv`: 3,918 rows. Columns used:
`short_Tag`, `season` (2019 = 2019–20, 2020 = 2020–21, 2021 = 2021–22), `Block_no`,
`ProteinAsis.`, `StarchAsis.`, `BetaglucansAsis.`, `FiberAsis.`.

## Filters

These are the same as section 2 of the 00_THIN script (checked against the script):

| filter | applied here |
|---|---|
| Drop cultivated checks (`short_Tag` starting with c / C / M: Morex, Clipper) | yes |
| Drop site 04 (Hachola) accessions (`short_Tag` starting with `HS04`) | yes |
| Drop block 4 of 2019 | **no**, because block 4 is what this check tests |
| Starch + 8.85 in 2021 | **no** (see below) |
| Per-season mean centring | **no** for the raw analysis (see below) |

**Why the starch +8.85 shift and the centring are not needed.** Both add or subtract
one constant for all values within a season. A constant within a season does not
change a within-season correlation, so the raw-value correlations would be
identical with or without these steps. The ±3 SD sensitivity step centres within
season itself, which removes the 2021 starch shift, so the set of dropped
outliers would also be the same with or without it.

After filtering there are 3,770 records and 290 accessions. Each season × block
has 290 records, with one plant per accession per block (the script checks
that no accession appears twice in a block). Seasons 2019–20 and 2020–21 have
4 blocks each; 2021–22 has 5.

## Method

For each season × trait:

1. Reshape the data to an accession × block matrix.
2. For every pair of blocks, compute the Pearson correlation (`cor.test`) across
   the accessions that have values in both blocks. Report n, r, the 95% CI
   (Fisher z) and the two-sided P.
3. Run two analyses:
   - **raw** (primary): raw values.
   - **sd3** (sensitivity): the project's per-trait ±3 SD outlier removal, which
     reproduces the logic of `clean_trait_local()`. Values are centred within
     season, pooled across seasons, NAs are removed, and values outside
     mean ± 3 SD are set to NA. One difference from the 00_THIN script is that
     the season means and the pooled SD here include block 4 of 2019, because
     that block is kept in this check. Values removed per trait (of 3,453
     non-NA values each): protein 13, starch 14, β-glucan 19, fiber 26
     (from the run log).
4. Summarise per season × trait × analysis: the mean, min and max r of the pairs
   that involve block 4, compared with the pairs among blocks 1–3. In 2021–22,
   which has 5 blocks, the block-4 pairs are 1-4, 2-4, 3-4 and 4-5. Pairs 1-5, 2-5
   and 3-5 fall into neither group; they are listed in `interblock_correlations_all.tsv`.

## How to run

```bash
cd /mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/12_USED_adittional_checks/01_block4_2019_interblock_correlations
export TMPDIR=/mnt/data/shahar/.tmp
/usr/bin/Rscript scripts/01_interblock_correlations.R > logs/01_interblock_correlations.log 2>&1
```

Packages: ggplot2, patchwork, ragg (for the figure only). The script sets TMPDIR itself.

## Outputs

| file | content |
|---|---|
| `results/tables/interblock_correlations_all.tsv` | every season × trait × analysis × block pair: n, r, ci_low, ci_high, p (176 rows) |
| `results/tables/interblock_summary.tsv` | per season × trait × analysis: n_blocks; count and mean/min/max r of block-4 pairs; count and mean/min/max r of pairs among blocks 1–3; n range |
| `results/tables/Table_OnlineResource_2019_block_correlations.tsv` | 2019–20, raw values: 6 block pairs × 4 traits (r and n); draft Online Resource |
| `results/tables/results_numbers.txt` | key numbers in prose |
| `results/figures/Fig_interblock_correlation_heatmaps.png` | block × block heatmaps of r (raw), 600 dpi, 174 × 170 mm |
| `logs/01_interblock_correlations.log` | full run log |

**Figure caption (draft).** Pearson correlations between blocks for grain protein,
starch, β-glucan and fiber content (raw values) in the wild-barley common garden.
**a** 2019–20, **b** 2020–21, **c** 2021–22. Each cell gives r across the
accessions measured in both blocks. The font is Liberation Sans, which has the
same metrics as Arial (Arial is not installed on this server). Text is 8–12 pt
at final size.

## Results in brief

All numbers come from `interblock_summary.tsv` and `results_numbers.txt`.

**2019–20, raw.** This table compares r for block 4 against blocks 1–3 with r
among blocks 1–3 (range, with the mean in brackets):

| trait | block 4 vs blocks 1–3 | among blocks 1–3 |
|---|---|---|
| Protein | 0.327–0.368 (0.352) | 0.343–0.426 (0.376) |
| Starch | 0.461–0.503 (0.487) | 0.524–0.594 (0.568) |
| β-glucan | 0.268–0.344 (0.301) | 0.605–0.641 (0.625) |
| Fiber | 0.154–0.300 (0.215) | 0.215–0.312 (0.259) |

Each pair has n = 283–288 accessions. P values for the 24 correlations range
from 4.26e-34 to 0.00929.

**2019–20, sd3 (sensitivity).**

| trait | block 4 vs blocks 1–3 | among blocks 1–3 |
|---|---|---|
| Protein | 0.305–0.368 (0.345) | 0.329–0.403 (0.363) |
| Starch | 0.466–0.516 (0.497) | 0.524–0.591 (0.567) |
| β-glucan | 0.240–0.320 (0.272) | 0.605–0.641 (0.625) |
| Fiber | 0.152–0.370 (0.242) | 0.270–0.313 (0.289) |

In both analyses, the mean r of the block-4 pairs is lower than the mean r
among blocks 1–3 for all four traits. The ranges of the two groups do not
overlap for starch and β-glucan, but they do overlap for protein and fiber.

**Reference seasons, raw.** These are mean r values: block-4 pairs / pairs among blocks 1–3.

| trait | 2020–21 (4 blocks) | 2021–22 (5 blocks) |
|---|---|---|
| Protein | 0.409 / 0.398 | 0.211 / 0.277 |
| Starch | 0.630 / 0.645 | 0.501 / 0.510 |
| β-glucan | 0.517 / 0.510 | 0.189 / 0.259 |
| Fiber | 0.425 / 0.501 | 0.153 / 0.209 |

The ranges and the sd3 values are in `interblock_summary.tsv`. In 2021–22, each
pair has n = 171–223, because fewer grain samples were measured by NIR.

**Sanity check against older drafts.** Older drafts quoted, without a source,
β-glucan r = 0.328–0.391 for block 4 and 0.670–0.682 among blocks 1–3 in
2019–20. This analysis does **not** reproduce those values. It gives
0.268–0.344 and 0.605–0.641 (raw), or 0.240–0.320 and 0.605–0.641 (sd3).
Nothing was tuned to try to match them.
