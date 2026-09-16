# 01_scripts/

Five pipeline scripts, run in order. All parameters live in
[`../00_config/config.yaml`](../00_config/config.yaml) and nowhere else.

| # | script | does | runtime |
|---|---|---|---|
| 0 | `00_build_gene_windows.R` | step-03 `candidate_genes.tsv` → `gene_windows.tsv` (gene ± `window_bp`), plus the 290-sample keep-list | s |
| 1 | `01_make_raw_per_gene_vcfs.sh` | `bcftools view` per gene window from the **raw** source VCF | ~18 s |
| 2 | `02_make_imputed_per_gene_vcfs.sh` | same from the **imputed** VCF (staged bgzip, one-time) | ~17 s |
| 3 | `03_run_crosshap_pipeline.R` | **the engine** — crosshap, figures, one KW test per gene, per-trait BH | ~6 min |
| 4 | `04_results_summary.R` | parameters table + paste-ready prose for the manuscript | s |
| 5 | `05_collect_significant_genes.R` | isolate the BH-significant genes into `Significant_genes/`, ranked, with both PDFs each | s |
| 6 | `06_named_gene_figures.R` | **runs after step-05 annotation** — re-render each significant gene's figures with its functional name in the header and merge violin+heatmap into one PDF under `Significant_genes/by_trait/<trait>/` | ~1 min |

## `R/` — helpers, unchanged from v1

| file | role |
|---|---|
| `run_crosshap.R` | **the proven core.** LD from the *imputed* per-gene VCF (complete r² matrix); haplotyping on the *raw* VCF (observed genotypes); variants intersected on CHROM:POS; PLINK called with `--keep-allele-order`. Do not "simplify" this split — it is deliberate. |
| `utils.R` | logging, event log, failure table, gene-window reader |
| `plot_combined_pdf.R`, `plot_heatmaps.R` | the tree+violin and heatmap renderers |

`03_run_pipeline_logged.sh` is a thin timestamped wrapper around
`03_run_crosshap_pipeline.R`, kept for long runs under `screen`.

**Archived 2026-09-10:** the v1 downstream scripts (`05_publication_plots*`, `06*`,
`07_build_publication_table.R`) moved to `_ARCHIVE_v1_epsilon_sweep_2026-09-08/`. They
read `gene_shortlist.csv` / `gene_annotation_review.csv`, which no longer exist.
Functional annotation is now the separate top-level step
`../../05_USED_gene_annotation_analysis/`.

## `03_run_crosshap_pipeline.R` — what to know

- It **refuses to run** if `config.yaml` supplies more than one `epsilon_vector` or
  `mgmin_values` entry. This pipeline is single-configuration by design; sweeping
  parameters and then picking the best p-value is the exact failure mode it was
  rewritten to remove. Use the archived v1 script for exploration.
- Results are **cached** per gene at `04_runs/<run_id>/Cache/<trait>/<gene>/MGmin_2/HapObject.rds`.
  Delete the cache or set `overwrite_cache: true` to force a re-run.
- Failures are **classified**, not lumped: `no_marker_groups_at_epsilon_*` (DBSCAN found
  no clusters — a parameter-coverage outcome) is distinguished from
  `no_variants_in_window` and `fewer_than_MGmin_variants_*`.
- `kw_test_with_effects()` is local to this script rather than in `utils.R`, so
  `utils.R` stays exactly as v1 left it.

## `diagnostics/`

`eps_coverage_scan.R` — **not part of the pipeline.** Sweeps epsilon and reports how
many genes yield a usable grouping, using genotype-only criteria. Kept as the evidence
behind the ε = 0.6 choice.

> Caveat found while writing it: passing crosshap a *vector* of epsilons makes a failure
> at any one value abort the whole gene, which under-counts coverage. Measure a fixed
> epsilon by running the real pipeline at that value, not by vector-sweeping.

## `_ARCHIVE_v1_epsilon_sweep_2026-09-08/`

The superseded v1 scripts and config. Read-only, kept reproducible.
