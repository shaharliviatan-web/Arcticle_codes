# Materials and Methods

## Plant material, field experiment, and phenotyping

A wild barley (*H. spontaneum*) collection comprised of 300 accessions sampled at 30 different sites along Israel following a unique sampling design which minimized the confounding effect of population structure (Potapenko et al. 2026) was used in this study. Genomic data was generated and processed for this entire collection, thus a total of 7,110,996 variants were called (Potapenko et al. 2026).

## Phenotypic analysis

Phenotyping of the entire collection was conducted in three common-garden experiments over consecutive seasons (2019--2022) in a net-house at Mataim Farm, Hula Valley, in a randomized complete block design (Pintel 2022). Grains from five spikes per plant was harvested at each season and four nutritional traits were quantified using a near-infrared spectroscopy (NIR; DA 7250, Perten Instruments) under a wild-barley-adapted calibration anchored to laboratory references (Pintel 2022). Ten accessions (site number 4 -- Hula) were missing in the phenotyping experiments and thus were excluded from all analyses. In addition, one complete block in season 2019 was strongly impacted by experimental variation as indicated from a poor correlation with all other blocks in the same season, therefore this block was removed from all downstream analyses. In each experiment, the four nutritional traits were scored and six morphological and phenological traits (flowering time, tillers, plant height, grain weight, grain number, and spike length). Prior to model fitting, each trait was mean-centered within season by subtracting the corresponding seasonal mean. Observations with missing values were excluded, and values exceeding ±3 SD from the pooled distribution of centered observations were removed as outliers.

To remove environmental effect and extract the genetic value in each trait the genotype best linear unbiased prediction scores (BLUPs) were calculated using a linear mixed model in lme4 v1.1-37 with genotype and season-by-block as crossed random effects.

Broad-sense heritability (H²) and variance partitioning were derived from a phenotypic variance partitioning which included the genetic variance (*V*~G~), season-by-block (*V*~E~), genotype-by-season as random effect (*V*~G×E~), and residual (*V*~R~) components. The broad-sense heritability was computed following (Holland et al., 2003):

$$H^{2} = \frac{V_{G}}{V_{G} + V_{E} + V_{G \times E} + V_{R}}$$

This model was performed for each trait, and the per-genotype reaction norms were generated.

Environmental associations were assessed at the collection-site level by averaging genotype BLUPs within each of the 29 sites. After excluding silt and elevation because of strong collinearity (\|r\| \> 0.8) with sand and March temperature, respectively, two-sided Pearson correlations were calculated between four nutritional traits and the eight retained environmental variables; Spearman correlations were used as a sensitivity analysis.

## Genome-wide association and identification of candidate QTLs

A GWAS was performed for each nutritional trait with EMMAX (Kang et al., 2010) using the BLUP scores as phenotypes. To control []{dir="rtl"}population structure in the model, a PCA was computed from a LD-pruned subset of 590,462 SNPs. Pruning was performed in PLINK v1.90 using a 50-SNP sliding window, a 5-SNP step size, and an LD threshold of r² = 0.2. The first three principal components were included as fixed covariates in the model. In addition, a kinship matrix was calculated from the SNP dataset for all accessions and was included as random effect.

Genome-wide significance threshold was corrected with Bonferroni over the 590,462 independent tests (pruned SNPs) to avoid overcorrection of the results. We used a permissive α = 0.10 to compensate for the small panel (n = 290) compared to the large SNP dataset (7M). For each trait analysis, Manhattan and QQ plots were generated and the genomic inflation factor (λ~GC~) was computed to evaluate the effectiveness of controlling the false positive rate.

To determine the size of a window to search for potential candidate genes we calculated the genome-wide LD decay within 2.5 Mb windows (1-kb bins, pair-count-weighted LOESS) indicating a distance of 200Kb to reach linkage equilibrium (*r*² = 0.2). Candidate genes were searched within the LD distance by intersected the significantly associated regions with the Morex V3 annotation (Yates et al., 2022) using bedtools v2.30.0.

## Haplotype analysis

To explore the haplotype structure around candidate genes we first defined their size using local LD analysis. To perform this, the SNPs dataset was imputed with the sNMF algorithm in LEA v3.6.0 defining the number of populations (K = 3) in accordance with the number ecotypes for wild barley in Israel (Hübner et al., 2009, 2013). Ten runs were performed and the lowest-cross-entropy run retained. Haplotype analysis was performed with crosshap v1.4.0 (Marsh et al., 2023) using the raw VCF under two minimum marker-group sizes (MGmin = 2 and 3) and seven ε thresholds (0.05--0.85; minHap = 9). The unresolved haplotype group was excluded and association of each haplotype and the corresponding trait was tested using a Kruskal--Wallis. The results were corrected within each gene, by collapsing duplicated solutions following the Holm's step-down procedure across the remaining Kruskal--Wallis tests and retaining the smallest P value as representative gene-level signal. The representative P values were then corrected separately across genes using the FDR and Bonferroni. Significant genes were annotated with DIAMOND BLASTP against UniProt Swiss-Prot, with InterProScan domains and BLASTP against NCBI nr.
