# 07_fiber_starch_tradeoff_direction

Direction test for the shared **fiber / starch** locus on chromosome **7H (~573.6 Mb)**,
wild-barley grain-quality GWAS (TAG paper). 290 wild barley accessions (*Hordeum vulgare*
ssp. *spontaneum*, Southern Levant), Morex V3. This step sits **on top of the crosshap
haplotype analysis** (`../04_runs/candidate_genes_1000bp_V1`) and answers two specific
questions about the three genes that passed the statistical filter for **both** traits.
No biological interpretation here -- numbers only.

## Why this step exists

Earlier steps established that at 7H ~573.6 Mb the fiber and starch GWAS signals are the
same locus: the fiber lead SNP `7H:573606306` (**significant**, -log10p 7.146) and the
starch lead SNP `7H:573606460` (**marginal**, -log10p 6.195) are 154 bp apart. In the FDR
candidate-gene annotation master (`../../05_USED_gene_annotation_analysis/08_USED_annotation_master`)
**three genes carry both a fiber row and a starch row** -- they passed the statistical
filter for both traits:

| gene_id | annotation (final_call) | tier |
|---|---|---|
| HORVU.MOREX.r3.7HG0729020 | Probable high-affinity nitrate transporter 2.4 | HIGH |
| HORVU.MOREX.r3.7HG0729090 | Pentatricopeptide repeat-containing protein At3g61360 | MEDIUM |
| HORVU.MOREX.r3.7HG0729100 | NifU-like protein 3, chloroplastic | HIGH |

Dual association alone does not prove a tradeoff. This step asks whether the **same
haplotype** that raises fiber **lowers** starch.

## Questions tested

- **Q1 -- Identical grouping?** Are the crosshap haplotype groups the same in the fiber run
  and the starch run? crosshap clusters on **genotype only** (phenotype is never used to
  build groups), so they must be identical; this verifies it per accession.
- **Q2 -- Inverse direction?** Within each gene's shared haplotype partition, is the
  highest-mean-**fiber** group the lowest-mean-**starch** group?

## Inputs (read-only)

- `../04_runs/candidate_genes_1000bp_V1/Cache/{fiber,starch}/<gene>/MGmin_2/HapObject.rds`
  -- crosshap objects. `$HapObject[[<Haplotypes_MGmin#_E#>]]$Indfile` gives, per accession:
  `Ind`, `hap` (`0` = unassigned; `A`,`B`,... = haplotype group), `Pheno` (BLUP of that
  run's trait).
- `../04_runs/candidate_genes_1000bp_V1/Stats/per_test_stats.csv` -- source of the per-gene
  **best** (most significant) crosshap setting (`is_best_within_gene == TRUE`). For all three
  genes fiber and starch agree on the same best setting (they must -- the partition is shared).
- `/mnt/data/shahar/gwas_barley/data/inputs/{fiber,starch}_corrected_V3.pheno` -- BLUP
  phenotypes, 290 accessions (used only for the overall correlation anchor).

Per-gene haplotype partition analysed (the per-gene best, MGmin = 2):

| gene | setting | groups (sizes) |
|---|---|---|
| 7HG0729020 | eps 0.80 | 4 (17 \| 42 \| 49 \| 111) |
| 7HG0729090 | eps 0.85 | 2 (116 \| 120) |
| 7HG0729100 | eps 0.85 | 4 (10 \| 15 \| 88 \| 93) |

`hap == "0"` (unassigned individuals) is dropped from the group means; n_assigned is reported.

## Run

```bash
bash scripts/run_all.sh
```

`scripts/01_tradeoff_direction.R` (R 4.1.2, `/usr/bin/Rscript`). TMPDIR pinned to
`/mnt/data/shahar/.tmp`; absolute paths; ASCII-only output; nothing written to `/tmp` or
`$HOME`. Log: `logs/01_tradeoff_direction.log`.

## Output tables (`results/tables/`)

- **shared_genes.tsv** -- the 3 genes, annotation + tier, fiber/starch lead SNP and class,
  MGmin/epsilon used.
- **haplotype_group_means.tsv** -- one row per (gene x haplotype group): `n`, `mean_fiber`,
  `sd_fiber`, `median_fiber`, `mean_starch`, `sd_starch`, `median_starch`. Sorted by
  descending mean_fiber.
- **tradeoff_summary.tsv** -- one row per gene: `n_individuals`, `n_assigned`, `n_groups`,
  `grouping_identical`, `n_mismatch`, `kw_p_fiber`, `kw_p_starch` (Kruskal-Wallis on the
  shared partition), `top_fiber_hap`, `top_fiber_starch_rank` (1 = lowest starch),
  `direction_inverse`, `pearson_groupmeans`, `spearman_groupmeans` (corr of the per-group
  mean_fiber vs mean_starch).
- **overall_phenotypic_correlation.txt** -- fiber vs starch BLUP correlation across all 290
  accessions (anchor for the locus-level direction).

## Results (this run)

**Q1 -- grouping identical: YES.** `grouping_identical = TRUE`, `n_mismatch = 0` / 290 for
all three genes. The fiber and starch runs partition the accessions into exactly the same
haplotype groups; only the overlaid BLUP (and thus the KW p) differs. So "shared locus" is
the *same haplotype partition*, not merely overlapping windows.

**Q2 -- inverse direction: YES for all three genes.** The highest-mean-fiber haplotype is
the lowest-mean-starch haplotype (`top_fiber_starch_rank = 1`, `direction_inverse = TRUE`).

| gene | n_groups | KW p fiber | KW p starch | Pearson (group means) | Spearman |
|---|---|---|---|---|---|
| 7HG0729020 (nitrate transporter 2.4) | 4 | 1.0e-09 | 1.7e-09 | -0.976 | -0.40 |
| 7HG0729090 (PPR At3g61360)           | 2 | 1.1e-06 | 1.8e-10 | -1.000 | -1.00 |
| 7HG0729100 (NifU-like)               | 4 | 1.1e-04 | 7.0e-08 | -0.965 | -1.00 |

Per-group means (high-fiber group in **bold**; full values in `haplotype_group_means.tsv`):

- **7HG0729020:** group **C** (n=42) mean fiber **+0.304** / mean starch **-6.767**;
  B +(-0.073)/+2.059, A -0.076/+0.889, D -0.156/+1.510. The 3 low-fiber groups are not
  perfectly monotone among themselves, so Spearman is only -0.40 although Pearson is -0.976;
  the dominant high-fiber/low-starch contrast is unambiguous (starch rank 1 of 4).
- **7HG0729090:** group **A** (n=120) fiber **+0.094** / starch **-2.994**; group B (n=116)
  fiber -0.115 / starch +2.390. Clean two-group inversion (Pearson and Spearman = -1).
- **7HG0729100:** group **A** (n=93) fiber **+0.107** / starch **-3.090**; B (n=88)
  -0.112/+1.962; C (n=15) -0.123/+4.228; D (n=10) +0.011/-0.976. Monotone inverse for the
  major groups.

**Anchor:** overall phenotypic correlation fiber vs starch (BLUP, 290 accessions) =
**Pearson -0.780, Spearman -0.718**. The locus-level direction matches the genome-wide
phenotypic relationship.

## Scope / what this step does and does not show

- Shows: a single genotype-defined haplotype partition at 7H ~573.6 Mb is significantly
  associated with both fiber and starch, and the high-fiber haplotype is the low-starch
  haplotype in every one of the three shared genes.
- Does not show: which of the three genes is causal (they lie in one LD block and yield the
  same partition, so they are statistically indistinguishable here), nor whether the inverse
  direction is fully independent of the compositional coupling of the two traits. Those are
  outside this step.

## Reproducibility

R 4.1.2; `/usr/bin/Rscript`. All paths absolute. `TMPDIR=/mnt/data/shahar/.tmp`. Tables are
plain TSV; re-running `run_all.sh` regenerates them deterministically from the cached
crosshap objects.
