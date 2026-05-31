# Publication folder — BLUE × 10 PCs

## Chosen configuration

- **Phenotype correction**: BLUE (Best Linear Unbiased Estimator)
- **Number of PCs used as covariates in EMMAX**: 10
- **Significance threshold**: Bonferroni α = 0.10 on the LD-pruned 590,462 SNP set → −log10p = **6.7712**
- This set shows **only** the Bonferroni threshold line at −log10p = 6.7712. The lead-SNP table lists only SNPs above this level.

## Reference numbers

- Total SNPs tested: **7,110,996**
- LD-pruned SNPs used for population-structure correction: **590,462**
- Sample size: **290** wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*) from the Southern Levant
- λ_GC range across the 4 traits (BLUE × 10 PCs): see `tables/per_trait_summary.tsv`
- Reference genome: Morex V3

## Directory tree

```
publication_BonfOnly/
├── main_figure/
│   └── Figure_main_manhattan_qq.{pdf,png}    — composed 4 traits × (Manhattan + QQ), TAG ~174 mm
├── individual_panels/
│   ├── manhattan_BLUE_pc10__{trait}.{pdf,png}  (4 traits × 2 formats = 8 files)
│   └── qq_BLUE_pc10__{trait}.{pdf,png}         (4 traits × 2 formats = 8 files)
├── pc_selection/
│   ├── Figure_S1_PC_scree.{pdf,png}              — 20-PC scree
│   ├── Figure_S2_lambda_vs_npcs.{pdf,png}        — λ_GC vs n_PCs (PC-count justification)
│   ├── Figure_S3..S5_PCA_scatter_*.{pdf,png}     — PC1×2, PC1×3, PC2×3
│   ├── Table_S2_pc_variance.tsv                  — 20 rows: PC, eigenval, %% variance, cumulative %%
│   ├── Table_S3_lambda_vs_npcs_long.tsv          — 24 rows (BLUP+BLUE × pc3/5/10 × 4 traits)
│   └── Table_S4_lambda_vs_npcs_wide.tsv          — 8 rows × pc3/pc5/pc10 columns
├── tables/
│   ├── lead_snps.tsv                  — significant SNPs (MAF, alleles, β/SE, LD window)
│   ├── per_trait_summary.tsv          — 4 rows: trait, n_tested, λ, n_above_thr, top SNP
│   └── top15_per_chr__{trait}.tsv     — top 15 per chromosome (15 × 7 = 105 rows each)
├── supp/
│   ├── Supp_Figure_S1_sensitivity_QQ_grid.{pdf,png}
│   ├── Supp_Figure_S2_sensitivity_Manhattan_grid.{pdf,png}
│   └── Supp_Table_S1_24run_summary.tsv  — full 24-run sensitivity defends BLUE × 10 PCs
├── README.md
└── cross_trait_interpretation.md
```

## Methods (short, for the manuscript)

- PCA and aIBS kinship matrix were computed on an LD-pruned subset (--indep-pairwise 50 5 0.2) of 590,462 SNPs.
- GWAS was performed using EMMAX (mixed model) on all 7,110,996 SNPs, including the top 10 PCs as fixed-effect covariates.
- Sensitivity analysis across 24 cells (4 traits × {BLUP,BLUE} × {3,5,10} PCs) was used to select the configuration; BLUE × 10 PCs minimised the mean |λ − 1| across the four traits without any trait falling below the 0.95 deflation floor.
- Significance threshold: Bonferroni at α = 0.10, corrected for the number of independent (LD-pruned) tests (n = 590,462), yielding −log10p = 6.77.

## Citations to use in methods

- Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006). *Principal components analysis corrects for stratification in genome-wide association studies.* Nature Genetics 38, 904–909.
- Patterson N, Price AL, Reich D (2006). *Population structure and eigenanalysis.* PLoS Genetics 2(12): e190.
- Kang HM, Sul JH, Service SK, Zaitlen NA, Kong S, Freimer NB, Sabatti C, Eskin E (2010). *Variance component model to account for sample structure in genome-wide association studies.* Nature Genetics 42, 348–354. (EMMAX)

## What is NOT here (downstream work)

1. **Candidate gene count per ±188 kb window** around each lead SNP — requires the Morex V3 GFF annotation + `bedtools intersect`. Plan to add as `tables/candidate_genes_per_lead.tsv` in a follow-up step.
2. Final manuscript prose and biological interpretation — see `cross_trait_interpretation.md` for a suggestive starting draft.
