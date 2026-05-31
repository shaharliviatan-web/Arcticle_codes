# Publication folder — BLUE × 3 PCs (Bonferroni only)

## Chosen configuration

- **Phenotype value**: BLUE (Best Linear Unbiased Estimator)
- **PCs as EMMAX covariates**: 3
- **Significance threshold**: Bonferroni α = 0.10 on the LD-pruned 590,462-SNP set → −log10p = **6.7712**
- This folder shows **only** the Bonferroni threshold line. `tables/lead_snps.tsv` lists only SNPs above this level.

## Reference numbers

- Total SNPs tested: **7,110,996**
- LD-pruned SNPs used for structure correction: **590,462**
- Sample size: **290** wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*), Southern Levant
- λ_GC across the 4 traits (BLUE × 3 PCs): **0.9686–0.9847** (see `tables/per_trait_summary.tsv`)
- Bonferroni-passing SNPs (all traits combined): **17**
- Reference genome: Morex V3

## Why 3 PCs

3 PCs were chosen as a parsimonious correction for population structure. The leading components capture the bulk of the resolvable structure (PC1–PC3 ≈ 14.3% of total variance; the per-PC contribution flattens after the first few components — see `pc_selection/`). Three PCs already bring λ_GC close to 1 for all four traits (0.9686–0.9847), indicating that residual stratification is adequately controlled; adding more PCs changed λ_GC only marginally (see `pc_selection/Table_S4_lambda_vs_npcs_wide.tsv`). The full 24-cell sensitivity analysis (4 traits × {BLUP,BLUE} × {3,5,10} PCs) is in `supp/` so robustness to this choice can be assessed.

## Directory tree

```
publication_BonfOnly_BLUE_3PC/
├── main_figure/
│   └── Figure_main_manhattan_qq.{pdf,png}   — 4 traits × (Manhattan + QQ), TAG full page, no titles
├── individual_panels/
│   ├── manhattan_BLUE_pc3__{trait}.{pdf,png}  (clean, no burned-in labels)
│   └── qq_BLUE_pc3__{trait}.{pdf,png}
├── pc_selection/   — scree, λ-vs-N_PCs, PCA scatters, variance tables
├── tables/
│   ├── lead_snps.tsv             — Bonferroni-passing SNPs (MAF, alleles, β/SE, ±188 kb LD window)
│   ├── per_trait_summary.tsv     — 4 rows: trait, n_tested, λ, n_above_Bonf, top SNP
│   └── top15_per_chr__{trait}.tsv
├── supp/           — 24-run sensitivity grids + summary table
├── README.md
└── cross_trait_interpretation.md
```

## Methods (short, for the manuscript)

- PCA and aIBS kinship were computed on an LD-pruned subset (--indep-pairwise 50 5 0.2) of 590,462 SNPs.
- GWAS used EMMAX (mixed model) on all 7,110,996 SNPs, with the top 3 PCs as fixed-effect covariates.
- A 24-cell sensitivity analysis (4 traits × {BLUP,BLUE} × {3,5,10} PCs) informed the configuration; BLUE × 3 PCs was adopted for parsimony with adequate genomic control (λ_GC near 1 for all traits).
- Significance threshold: Bonferroni at α = 0.10 corrected for the number of independent (LD-pruned) tests (n = 590,462) → −log10p = 6.77.

## Citations for methods

- Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006). *Principal components analysis corrects for stratification in genome-wide association studies.* Nature Genetics 38, 904–909.
- Patterson N, Price AL, Reich D (2006). *Population structure and eigenanalysis.* PLoS Genetics 2(12): e190.
- Kang HM, Sul JH, Service SK, Zaitlen NA, Kong S, Freimer NB, Sabatti C, Eskin E (2010). *Variance component model to account for sample structure in genome-wide association studies.* Nature Genetics 42, 348–354. (EMMAX)

## Not included (downstream)

1. **Candidate gene count per ±188 kb window** around each lead SNP — needs the Morex V3 GFF + `bedtools intersect`.
2. Final manuscript prose — see `cross_trait_interpretation.md` for a suggestive starting draft.
