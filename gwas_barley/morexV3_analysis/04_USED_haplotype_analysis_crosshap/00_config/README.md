# 04_USED_haplotype_analysis_crosshap

CrossHap haplotype analysis of the step-03 candidate genes, to find -- per trait -- genes whose
haplotype groups differ significantly in phenotype. This is the statistical filter applied before
the manual biological filtering; the surviving genes get publication plots (script 05) and feed the
later wild+elite comparison.

## Context

Wild-barley grain-quality GWAS for a TAG paper. 290 wild barley accessions (Hordeum vulgare ssp.
spontaneum, Southern Levant), Morex V3 reference, 4 traits (beta-glucan, fiber, protein, starch).
GWAS = EMMAX (aIBS kinship), BLUP phenotypes + 3 PCs, Bonferroni alpha=0.10 (-log10p = 6.7712).
Step 03 (`03_USED_candidate_genes_around_leading_snps`) produced 108 candidate genes in 20 loci
(11 significant, 9 marginal) within 200 kb of the lead/marginal SNPs. Each gene here is tagged with
its origin lead SNP and significant|marginal class.

## What this step does

For every candidate gene, run CrossHap on the gene +/- 1000 bp variants, build haplotypes across a
range of epsilon and two MGmin thresholds, and test whether haplotype groups differ in phenotype by
an omnibus Kruskal-Wallis test, with within-gene and per-trait multiple-testing correction. Produces
screening plots (CrossHap tree, per-epsilon CrossHap viz, violins, genotype heatmaps) and statistics
tables, plus a single decision sheet (`gene_shortlist.csv`).

## Inputs

- Genes + windows: `00_config/gene_windows.tsv` (built by 00 from step-03 `candidate_genes.tsv`).
  Window = gene_start-1000 .. gene_end+1000 (region cut). 108 genes (betaglucan 72, fiber 18,
  protein 1, starch 17).
- Phenotype: `{trait}_corrected_V3.pheno` (the BLUP phenotypes the GWAS used). 290 samples, Ind = IID.
- Raw VCF source: `data/inputs/morexV3_with_ids.vcf.gz` (300 samples; REF == Morex V3 genome, no flip).
- Imputed VCF (for LD only): `USED_imputation/morexV3_final_imputed.vcf` (LEA sNMF K=3); see
  `02_imputation/PROVENANCE.md`.
- 290 analysis samples = 300 minus the 10 in `samples_to_remove_V3.txt` (HS0401-HS0411).

## Directory layout

```
00_config/   config.yaml, gene_windows.tsv, samples_keep_290.txt, README.md, final_genes.tsv.template
01_scripts/  00_build_gene_windows.R  01_make_raw_per_gene_vcfs.sh  02_make_imputed_per_gene_vcfs.sh
             03_run_crosshap_pipeline.R  03_run_pipeline_logged.sh  04_make_shortlist.R
             05_publication_plots.R  05_publication_plots.sh   R/ (utils, run_crosshap, plot_*)
02_imputation/   imputation script copies + PROVENANCE.md (imputation is reused, not re-run)
03_per_gene_vcfs/   raw_1000bp/{trait}/{gene}.vcf.gz   imputed_1000bp/{trait}/{gene}.vcf.gz   manifests
04_runs/<run_id>/   Cache/ CombinedPDF/ Heatmaps/ Logs/ Stats/ tmp/
05_results/   tables/   publication_plots/
```

## Run order

```bash
# (once) build inputs
Rscript 01_scripts/00_build_gene_windows.R
bash    01_scripts/01_make_raw_per_gene_vcfs.sh
bash    01_scripts/02_make_imputed_per_gene_vcfs.sh
# main analysis (run inside GNU screen; ~1 h for 108 genes x 2 MGmin)
bash    01_scripts/03_run_pipeline_logged.sh
# decision sheet
Rscript 01_scripts/04_make_shortlist.R
# publication plots for the curated final genes (after biological filtering)
bash    01_scripts/05_publication_plots.sh  00_config/config.yaml  00_config/final_genes.tsv
```

## CrossHap parameters (config.yaml)

MGmin = 2 and 3; epsilon = 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 0.85; minHap = 9; hetmiss_as = "allele";
keep_outliers = false. RAW per-gene VCF drives haplotyping + genotype display; IMPUTED per-gene VCF
drives the PLINK LD (r2). PLINK = /usr/local/bin/plink with --keep-allele-order.

## Statistics (the discovery criterion)

1. For each gene x MGmin x epsilon: omnibus Kruskal-Wallis of phenotype across haplotype groups,
   after dropping haplotype 0, dropping samples with missing phenotype, and requiring > 1 group.
2. Within a gene, collapse duplicate tests (same raw p + n_groups + sorted group sizes), keeping the
   smallest epsilon, then larger MGmin.
3. Within-gene Holm correction across the unique tests; the representative test = smallest Holm p
   (deterministic tie-break: smallest raw p, then smallest epsilon, then larger MGmin).
4. Per trait (not pooled across traits), across-gene BH (FDR) + Bonferroni on each gene's best
   within-gene Holm p; significance threshold alpha = 0.05.

Violin plots are post-hoc/visual: pairwise Wilcoxon of every group vs the LARGEST haplotype group,
Holm-corrected; not the discovery statistic.

## Outputs (04_runs/<run_id>/Stats/)

- `gene_summary.csv`  one row per gene with >= 1 valid test: best MGmin/epsilon, raw p, within-gene
  Holm p, trait FDR + Bonferroni p, significance flags, and the gene's origin (lead_SNP, class, ...).
- `per_test_stats.csv`  one row per unique valid test (audit trail).
- `failures_notes.csv`  all invalid/failed tests (e.g. genes with too few SNPs, epsilons with no
  marker groups).
- `gene_shortlist.csv`  THE decision sheet: one row per candidate gene (tested or not) with origin
  lead SNP + class, n_snps_window, status, best MGmin/epsilon, all corrected p-values, significance
  flags, description, and paths to the best combined PDF + heatmap dir. Sorted by trait then FDR p.
- Screening plots: `CombinedPDF/<trait>/<gene>/MGmin_<M>/...pdf` (tree + summary pages + per-epsilon
  CrossHap viz + violin) and `Heatmaps/<trait>/<gene>/MGmin_<M>/...pdf` (one page per haplotype group).

## How to pick genes / epsilon / MGmin

Open `gene_shortlist.csv`, sort within trait by `trait_fdr_p`, inspect `significant_fdr` (and the
stricter `significant_bonferroni`). For a candidate, open its best combined PDF (tree to judge
epsilon stability, violins to see group separation) and its heatmap (haplotype structure). Record the
chosen genes + MGmin + epsilon in a `final_genes.tsv` (see template) and run script 05.

## Bug fixes vs the previous tool

- PLINK: absolute `/usr/local/bin/plink` (bare `plink` is broken v0.76) + `--keep-allele-order`.
- Per-gene VCFs cut with bcftools (preserve REF/ALT, no flip) instead of plink; single-IID samples;
  290-sample subset; raw/imputed sample-join asserted in the runner.
- CrossHap tree uses `crosshap::clustree_viz` (old code called `clustree::clustree_viz`, wrong
  namespace -> tree page failed).
- Robust VCF read (strip `##` lines so data.table::fread cannot mis-detect the separator).
- Violin reference group = largest haplotype group (was farthest-from-global-mean).
- TMPDIR pinned to /mnt/data/shahar/.tmp; temp under the run's tmp/; no /tmp or $HOME writes.

## What is NOT here

- No allele harmonization (all-wild analysis; harmonization is deferred to the later wild+elite step,
  which reuses these clean wild per-gene VCFs).
- CrossHap is the statistical filter only; biological relevance is judged by the user afterwards.
- 16 of 108 genes have < 2 SNPs in their +/- 1000 bp window (incl. the single protein gene) and
  cannot form haplotypes; they are logged to failures_notes / shown as no_valid_test in the shortlist.
- Functional descriptions exist for only ~12 of the 108 genes (Morex V3 projected annotation; no GO).

## Result summary (run candidate_genes_1000bp_V1)

108 candidate genes; 76 produced a valid haplotype test (32 had no valid test, e.g. too few SNPs
in the +/- 1000 bp window, including the single protein gene which has 0 SNPs). FDR 0.05 is the
primary gene filter; Bonferroni 0.05 reported alongside.

| trait      | genes | tested | no valid test | sig FDR 0.05 | sig Bonferroni 0.05 |
|------------|-------|--------|---------------|--------------|---------------------|
| betaglucan | 72    | 52     | 20            | 30           | 11                  |
| fiber      | 18    | 13     | 5             | 9            | 7                   |
| protein    | 1     | 0      | 1             | 0            | 0                   |
| starch     | 17    | 11     | 6             | 9            | 7                   |
| total      | 108   | 76     | 32            | 48           | 25                  |

Failures (370 total): 252 epsilon-level no_valid_haplotype_grouping (e.g. eps=0.05 with no marker
groups), 86 crosshap_failed (SNP-poor gene x MGmin), 32 genes with no valid test.

Next step is the user's biological filtering of the FDR-significant genes, using
`05_results/tables/gene_annotation_review.csv` (per trait: passed -> did-not-pass -> no-test, each
sorted by significance, with curated/GFF annotations joined and NEED TO ADD ANNOTATION flags), then
script 05 for publication plots of the final keepers.
