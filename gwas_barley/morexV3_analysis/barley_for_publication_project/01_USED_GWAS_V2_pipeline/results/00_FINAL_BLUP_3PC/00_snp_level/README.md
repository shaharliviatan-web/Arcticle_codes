# 00_snp_level — SNP-level results

Produced by `scripts/11_chosen_config_tables_v3.R`. Independent of any locus
definition: these are the per-SNP GWAS results for the locked configuration.

| file | contents |
|---|---|
| `significant_snps.tsv` | all 52 SNPs above −log10p 6.0454, with A1/A2, MAF, beta, SE, p. **`beta` = effect of A1 (minor allele)**, flipped from EMMAX 2026-09-22 (see `../README.md`); same for `top_beta` in `per_trait_summary.tsv` |
| `per_trait_summary.tsv` | per-trait counts, lambda_GC, top signal |
| `analysis_parameters.tsv` | every constant used |
| `top15_per_chr__<trait>.tsv` | top 15 SNPs per chromosome. **`beta` here is the raw EMMAX sign = effect of A2 (major allele)**; no allele columns |

**Significant SNPs only** — there is no marginal/suggestive class.
