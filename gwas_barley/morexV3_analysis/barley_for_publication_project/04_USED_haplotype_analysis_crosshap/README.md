# 04_USED_haplotype_analysis_crosshap

crosshap haplotype analysis of the step-03 candidate genes: per trait, does haplotype
group at a gene predict the grain-quality phenotype? 290 wild barley accessions
(*Hordeum vulgare* ssp. *spontaneum*, Southern Levant), Morex V3, 4 traits.

**Rebuilt 2026-09-08; current run `loci_LDspan_eps06_V4` (2026-09-09).** Earlier runs are
kept intact under `04_runs/`; the retired v1 scripts are in
`01_scripts/_ARCHIVE_v1_epsilon_sweep_2026-09-08/`.

**Start here → [`04_runs/loci_LDspan_eps06_V4/Stats/results_chapter_numbers.txt`](04_runs/loci_LDspan_eps06_V4/Stats/results_chapter_numbers.txt)**

### Run history (each is a different step-01 locus definition, not a method change)

| run | step-01 rule | genes in | tested | significant | on disk |
|---|---|---|---|---|---|
| `candidate_genes_1000bp_V1` | 200 kb single-linkage (retired) | 108 | 76 | 48 | kept — the original submission, and `06_publication_figures/` + `07_fiber_starch_*/` still read it |
| `loci_LDspan_eps06_V2` | clump 1 Mb, gap 50 kb, +2 peaks | 64 | 30 | 21 | **deleted 2026-09-10** |
| `loci_LDspan_eps06_V3` | clump 2 Mb, gap 60 kb, +2 peaks | 90 | 48 | 37 | **deleted 2026-09-10** |
| **`loci_LDspan_eps06_V4`** | **clump 2 Mb, gap 50 kb, iterative, no peaks** | **55** | **30** | **23** | **live** |

V2 and V3 were intermediate step-01 locus definitions, superseded and depended on by
nothing. To recreate either, set the step-01 parameters in that row, re-run step 03,
then `01_scripts/00_build_gene_windows.R` → `01` → `02` → `03_run_crosshap_pipeline.R`
with a matching `run_id`.

---

## The result

| | |
|---|---|
| Candidate genes in (from step 03) | 55 |
| **Genes tested** | **30** |
| Genes not testable | 25 |
| **Significant (BH q ≤ 0.05)** | **23 genes, in 12 loci** |
| Median haplotype groups per gene | 3 |

| trait | genes in | tested | not testable | significant | loci with a hit |
|---|---|---|---|---|---|
| beta-glucan | 7 | 6 | 1 | 6 | 4 of 4 |
| fiber | 32 | 19 | 13 | 12 | 5 of 6 |
| **protein** | **6** | **0** | **6** | **0** | **0** |
| starch | 10 | 5 | 5 | 5 | 3 of 3 |

**Protein drops out entirely in this run.** Its 2 significant loci yield 6 candidate genes,
and none is testable: 3 have **zero** SNPs in window, 1 has a single SNP, and 2 fail to
cluster at ε = 0.6 (5 and 2 SNPs). The 2 sub-threshold protein peaks that previously
supplied protein's only hit were dropped from step 01 by decision on 2026-09-09.

---

## What changed from v1, and why

v1 swept **7 epsilon × 2 MGmin** per gene, collapsed duplicate results, applied a
**within-gene Holm** correction, took the smallest Holm p as the gene's result, and then
applied **both BH and Bonferroni** across genes. All of that is gone. Three reasons:

1. **Double correction.** v1 fed Holm-*adjusted* p-values into BH. BH applied to
   FWER-adjusted values does not control FDR at any interpretable rate, so v1's
   `trait_fdr_p` was not an FDR. Here BH receives the **raw** p.
2. **Arbitrary penalty.** The Holm penalty a gene paid equalled its number of unique
   epsilon results — 1 to 13, a property of the gene's SNP structure, not of the
   hypothesis. 12 genes paid no penalty at all; 8 of those were called significant.
3. **Selection on the outcome.** Choosing the epsilon with the smallest p is a forking
   path. For **3 of v1's 4 publication genes the epsilon choice changed nothing** —
   they were flat across the entire sweep. Only PHT was sensitive, and its curated
   ε = 0.85 was the value that maximised its significance, off its stable plateau.

**Now: one grouping per gene, one test per gene, one correction.**

### What was kept, deliberately

The engineering is v1's and is good: caching, event logging, failure capture, and
especially the **raw/imputed split** in `R/run_crosshap.R` — LD is computed from the
**imputed** per-gene VCF (complete r² matrix), while haplotypes are called on the
**raw** VCF (observed genotypes), the two intersected on CHROM:POS. The plotting code
is untouched; the figures are the real product of this step.

---

## The two fixed parameters

Both are declared in [`00_config/config.yaml`](00_config/config.yaml) and must not be
tuned per gene.

**MGmin = 2** (DBSCAN `minPts`). All four v1 publication genes used it.

**epsilon = 0.6** (DBSCAN `eps`). Two things to know:

> **epsilon is not an r² cutoff.** crosshap calls `dbscan::dbscan(LD, eps, minPts)`
> where `LD` is a plain PLINK `--r2 square` **matrix**, not a distance object. So each
> SNP is a point whose coordinates are its r² against every other SNP, and epsilon is a
> **Euclidean radius in that profile space**. Its effective stringency therefore depends
> on how many SNPs are in the window.
>
> **Not stated in the manuscript** (user decision, 2026-09-23). This README previously said
> "State this in the methods"; that instruction is withdrawn. The Methods report ε = 0.6 and
> MGmin = 2 as fixed parameters and the genotype-only grounds for choosing them, without the
> parameter-space explanation. The point remains recorded here because it governs how any
> future run must interpret ε.

It was chosen on **genotype-only** criteria — crosshap never uses the phenotype to build
haplotype groups, so this selection is not a forking path. Evidence in
[`04_runs/loci_LDspan_eps06_V4/Diagnostics/epsilon_choice_supplementary.tsv`](04_runs/loci_LDspan_eps06_V4/Diagnostics/epsilon_choice_supplementary.tsv),
measured on **this run's own 55 candidate genes**:

| ε | genes testable | accessions assigned | median groups |
|---|---|---|---|
| 0.2 | 23 / 55 | 70.0% | 2 |
| 0.4 | 27 / 55 | 69.3% | 3 |
| **0.6 (used)** | **30 / 55** | **65.0%** | **3** |
| 0.8 | 34 / 55 | 60.9% | 4 |
| 1.0 | 33 / 55 | 60.3% | 4 |

*(measured on this run's own 55 candidate genes — see
`04_runs/loci_LDspan_eps06_V4/Diagnostics/epsilon_choice_supplementary.tsv`, with the per-gene
data in `epsilon_choice_per_gene.tsv` (55 genes × 5 ε); both are written by
`01_scripts/diagnostics/parameter_choice_tables.R` (2026-09-30), which re-runs crosshap per gene
and ε and computes no p-value. The ε = 0.6 rows reproduce the published run exactly. This table
is an Online Resource of the manuscript (M&M, Haplotype analysis).)*

Assignment rate and gene coverage trade off monotonically against each other: a low ε
assigns more accessions but in fewer genes (and with only ~2 groups), a high ε makes more
genes testable but leaves more accessions unassigned. **0.6 is the midpoint**, not the
maximum of either — that is the claim, and it is the one to use in the methods text.

A wider MGmin × epsilon grid scan (`01_scripts/diagnostics/param_grid_scan.R`) exists to
re-examine this on genotype-only criteria. It reports coverage and assignment only, never
p-values, so it cannot be used to tune toward significance.

**To change it:** edit `epsilon_vector` in `00_config/config.yaml` and re-run
`01_scripts/03_run_crosshap_pipeline.R` (~6 min). The pipeline refuses to start if more
than one epsilon or MGmin is given — parameter exploration belongs in the archived v1
script, not here.

---

## Statistics

One **Kruskal–Wallis** test per gene: BLUP phenotype ~ haplotype group, dropping
unassigned accessions (`hap 0`) and missing phenotypes, requiring ≥ 2 groups.
**Benjamini–Hochberg** on the raw p-values, **within each trait**, α = 0.05.

**Why BH is valid even though genes within a locus are in LD:** BH controls FDR under
positive regression dependence (Benjamini & Yekutieli 2001), which LD-induced
correlation satisfies. The problem LD creates is one of **reporting**, not validity —
significant genes cluster in loci, so `locus_summary.tsv` reports loci alongside genes
and the prose summary carries an explicit counting caveat.

**Genes with no valid grouping are excluded from the BH denominator** (no test was
performed) but are listed with a reason in `genes_not_tested.tsv`, so the real search
space stays visible:

| reason | genes | fixable by a different ε? |
|---|---|---|
| no marker groups at ε = 0.6 | 14 | possibly |
| no variants in the window | 7 | no |
| fewer than MGmin variants after raw/imputed intersection | 3 | no |
| crosshap internal error | 1 | no |

**10 of the 25 fail for lack of data**, which no parameter can repair.

**Effect sizes are reported alongside every p-value.** A KW p-value at a locus that was
selected for association with this very phenotype is close to guaranteed; η² and the
top-vs-bottom group difference in phenotype SD are what a reader can actually judge.

---

## Pipeline

```bash
Rscript 01_scripts/00_build_gene_windows.R      # step-03 genes -> gene_windows.tsv (+ 290-sample keep list)
bash    01_scripts/01_make_raw_per_gene_vcfs.sh      # ~18 s
bash    01_scripts/02_make_imputed_per_gene_vcfs.sh  # ~17 s
Rscript 01_scripts/03_run_crosshap_pipeline.R   # ~6 min: crosshap + KW + BH + figures
Rscript 01_scripts/04_results_summary.R         # paper-ready prose + parameters table
Rscript 01_scripts/05_collect_significant_genes.R  # the hits, isolated and ranked
Rscript 01_scripts/06_named_gene_figures.R         # named, merged per-gene PDFs by trait
```

Scripts 00-05 are the pipeline. **`06_named_gene_figures.R` runs after step 05
annotation**, because it needs each gene's functional name; re-run it whenever the
annotations change.

## Directory map

| path | contents |
|---|---|
| `00_config/` | `config.yaml` (**all parameters**), `gene_windows.tsv`, sample keep-list |
| `01_scripts/` | the 5 pipeline scripts + `R/` helpers — see `01_scripts/README.md` |
| `01_scripts/diagnostics/` | not part of the pipeline: `parameter_choice_tables.R` (the ε-choice table of the manuscript, 2026-09-30); older scans `eps_coverage_scan.R`, `param_grid_scan.R`, `protein_epsilon_rescue_scan.R` |
| `02_imputation/` | imputation provenance (untouched) |
| `03_per_gene_vcfs/` | per-gene raw + imputed VCFs, with manifests |
| `04_runs/loci_LDspan_eps06_V4/` | **the live run** — `Stats/`, `Significant_genes/`, `CombinedPDF/`, `Heatmaps/`, `Cache/`, `Logs/` |
| `04_runs/candidate_genes_1000bp_V1/` | the v1 run, frozen (v2/v3 deleted — see run history) |
| `05_results/`, `06_publication_figures/`, `07_fiber_starch_tradeoff_direction/` | **stale — see below** |

## Output tables (`04_runs/loci_LDspan_eps06_V4/Stats/`)

| file | one row per | notes |
|---|---|---|
| **`gene_results.tsv`** | tested gene | **the main table.** Grouping, KW, BH q, η², top-vs-bottom SD, figure paths. Replaces v1's `gene_shortlist.csv` + `gene_summary.csv` + `per_test_stats.csv` |
| `haplotype_groups.tsv` | gene × haplotype group | n, mean, median, SD of the phenotype — the effect-size detail |
| `genes_not_tested.tsv` | untestable gene | with a classified reason |
| `locus_summary.tsv` | locus | genes tested / significant, best gene |
| `per_trait_summary.tsv` | trait | |
| `analysis_parameters.tsv` | parameter | every constant, for the methods section |
| `results_chapter_numbers.txt` | — | all of the above as prose |
| `failures_notes.tsv` | failure event | diagnostics |

Plus **`Significant_genes/`** — a self-contained, browsable folder holding only the genes
that passed BH: the ranked table, per-group means, and both PDFs per gene with
rank-prefixed filenames (`01` = smallest q). Rebuilt by `01_scripts/05_collect_significant_genes.R`.

---

## ⚠ Downstream directories are now stale

These were built against v1 and **will not run** until repointed. Each carries a
`BREAKING_CHANGES.md` with the specifics.

| directory | what breaks |
|---|---|
| `05_results/` | v1 tables **archived 2026-09-10** to `05_results/_ARCHIVE_v1_tables_2026-09-10/`; nothing current here |
| `06_publication_figures/` | hard-codes `candidate_genes_1000bp_V1`; reads `gene_shortlist.csv` |
| `07_fiber_starch_tradeoff_direction/` | hard-codes the v1 run and `per_test_stats.csv`; builds the HapObject key from a per-gene epsilon |
| ~~`01_scripts/05_publication_plots*`, `06*.R`, `07_build_publication_table.R`~~ | **archived 2026-09-10** to `01_scripts/_ARCHIVE_v1_epsilon_sweep_2026-09-08/` — they read tables that no longer exist; superseded by step 05 and `06_named_gene_figures.R` |

**`00_config/final_genes.tsv` is also stale.** Its four curated genes were chosen from
v1's 108-gene set; three of them (Pho, AP2/ERF, BAHD) are not in the current 64-gene
candidate set, because step 03 now uses the LD locus span with no flanking window.
Only PHT survives. This was accepted deliberately — see step 03's README.
