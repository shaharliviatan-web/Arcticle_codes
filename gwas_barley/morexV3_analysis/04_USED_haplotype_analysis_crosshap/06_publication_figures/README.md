# 06_publication_figures

Publication-grade haplotype figures for the 4 final candidate genes, plus the
stats tables needed to write the results chapter. **Wild-only**: 290 wild-barley
accessions; NO elite/cultivated lines are read, plotted, or tested anywhere here.

Based on the legacy presentation tool
(`USED_Haplotype_Analysis/Wild_And_Elite/.../make_one_gene_presentation_png_NO_EliteInViolinPlot.R`)
but rebuilt clean on the step-04 results, with the elite rows / elite-guided
representative selection removed and TAG-oriented styling.

## What each figure is

A stacked two-panel figure per gene:

- **Top (violin):** BLUP phenotype distribution per haplotype group (crosshap
  group `0` = ungrouped, excluded), boxplot inside each violin, `n` per group on
  the x-axis, the omnibus **Kruskal-Wallis P** as subtitle, and pairwise
  **Wilcoxon rank-sum (Holm-adjusted)** significance brackets comparing every
  group to the **largest** group (the reference). Brackets show stars
  (`* <=0.05, ** <=0.01, *** <=0.001, **** <=1e-4, ns`); the exact p-values live
  in the CSVs (see below).
- **Bottom (heatmap):** one **representative wild accession per haplotype group**
  (fewest missing genotypes, then alphabetical), genotype shown across all SNPs
  in the gene window (REF = tan, ALT = slate, Missing = grey). SNP positions are
  not labeled (kept clean); the SNP count is in the CSV.

Colours are the **Okabe-Ito colourblind-safe** palette.

## Selected genes / parameters

Driven entirely by `config/figure_genes.tsv` (edit titles, y-axis labels, MGmin,
epsilon there). The selected MGmin x epsilon per gene matches each gene's
representative (best) grouping from step 04, except PHT where epsilon 0.85 was the
user-selected grouping (still FDR & Bonferroni significant).

| panel | symbol | gene | trait | MGmin | epsilon |
|-------|--------|------|-------|-------|---------|
| a | Pho  | 3HG0301750 | starch     | 2 | 0.40 |
| b | PHT  | 3HG0301710 | starch     | 2 | 0.85 |
| c | AP2/ERF | 3HG0299440 | betaglucan | 2 | 0.50 |
| d | BAHD | 7HG0642350 | fiber      | 2 | 0.20 |

> NOTE on panel c (AP2/ERF, 3HG0299440): an earlier draft of the curated final
> list labelled this gene "PME / pectin methylesterase". That was a mistake and
> has been corrected. The Morex V3 projected annotation ("APETALA2-like protein
> 5"), the step-05 Swiss-Prot best hit (Ethylene-responsive transcription factor
> ERF110) and the step-06 InterProScan evidence (Pfam PF00847 AP2 domain;
> InterPro AP2/ERF domain; GO transcription factor / ethylene signalling) all
> agree this is an **AP2/ERF ethylene-responsive transcription factor**, NOT a
> pectin methylesterase. The title is driven by `config/figure_genes.tsv`.

## Outputs

`results/figures/`
- `<serial>_<symbol>_<gene>.png` (600 dpi) and `.pdf` (vector, cairo)  -- per gene
- `composite_4panel.png` / `.pdf`  -- 2x2 composite with (a)-(d) tags

`results/tables/` (everything needed to write numbers into the results text)
- `<serial>_<symbol>_<gene>_summary.csv`  -- one row per gene: window (chr/start/
  end, n SNPs), MGmin/epsilon, n_total, n_groups, reference group, Kruskal-Wallis
  statistic / df / p, epsilon-squared effect size, lead SNP + significant/marginal
  class + distance to lead, within-gene Holm p, and the gene-level trait FDR &
  Bonferroni p with their significance flags (from step 04).
- `<serial>_<symbol>_<gene>_groups.csv`  -- one row per haplotype group: n, mean,
  sd, median, q25, q75, min, max; is_reference; Wilcoxon-vs-reference raw & Holm
  p + stars; and the heatmap representative accession with its ALT/REF/missing
  counts.
- `all_genes_summary.csv`, `all_genes_groups.csv`  -- the four genes combined.

## Reproduce

```
cd 06_publication_figures
TMPDIR=/mnt/data/shahar/.tmp /usr/bin/Rscript scripts/00_make_figures.R
```

Reads only the step-04 cache (`04_runs/candidate_genes_1000bp_V1/Cache/<trait>/
<gene>/MGmin_<M>/HapObject.rds`, which carries the raw per-gene VCF path and the
exact haplotyped variant IDs) and `Stats/gene_shortlist.csv`. The genotype matrix
is reconstructed from the raw VCF and asserted identical to the cached variant IDs
before plotting, so the heatmap SNP set matches what crosshap haplotyped on.

## Notes
- TMPDIR pinned to `/mnt/data/shahar/.tmp`; no `/tmp` or `$HOME` writes; no emojis.
- ggsave emits PNG via the `ragg` device (it overrides `type="cairo"` -- expected);
  PDFs use `grDevices::cairo_pdf` for embedded fonts / true vector output.
