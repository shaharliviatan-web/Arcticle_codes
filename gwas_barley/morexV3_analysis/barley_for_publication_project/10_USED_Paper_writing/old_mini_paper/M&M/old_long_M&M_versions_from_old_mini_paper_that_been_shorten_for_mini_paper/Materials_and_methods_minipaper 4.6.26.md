<!-- Readable md of "Materials_and_methods_minipaper 4.6.26.docx" (Word metadata: created and modified 2026-06-04). Converted 2026-09-29 with pandoc 3.11 for easier reading in the M&M blueprint session. The docx has no comments or tracked changes. Read-only reference: the .docx is the original. -->

Materials and Methods

Plant material and genotyping

A wild barley (*Hordeum spontaneum*) germplasm collection of 300 accessions, established in 2017, was used across the three growing seasons of this study. The accessions were sampled from 30 sites spanning the environmental and climatic gradients of Israel — soil type, elevation, precipitation, and temperature — selected from previous population-genetic studies of wild barley (Hübner et al., 2009, 2012, 2013; Hübner & Kantar, 2021). The sites span the three recognized ecotypes (northern, coastal, and desert) and the three intermediate ecotones. Ten individuals were collected per site at a minimum spacing of ten meters, geographical coordinates were recorded, and all accessions underwent two rounds of single-seed descent (SSD) and selfing.

All 300 accessions were whole-genome sequenced at an average coverage of 5× on an Illumina NovaSeq 6000 platform (Shermeister, 2021). Raw reads were quality-checked with FastQC v0.11.5 (Andrews, 2017) and trimmed of adapters and low-quality bases with fastp v0.20.0 (Chen et al., 2018). Cleaned reads were aligned to the *Hordeum vulgare* cv. Morex reference genome assembly v3 (Mascher et al., 2021) using BWA-MEM2 v2.2.1 (Vasimuddin et al., 2019) under default parameters, and duplicates were marked with Picard MarkDuplicates v3.4.0 (Broad Institute, 2019). Variants were then called following the GATK4 best-practices workflow (GATK v4.6.2.0; Poplin et al., 2017; Heldenbrand et al., 2019): per-sample calling in GVCF mode with HaplotypeCaller (ploidy = 2), reblocking with ReblockGVCF, cross-accession consolidation with GenomicsDBImport, and joint genotyping with GenotypeGVCFs (maximum two alternate alleles per site).

Post-genotyping filtering was performed with bcftools v1.13 (Danecek et al., 2021) and VCFtools v0.1.15 (Danecek et al., 2011), retaining only biallelic SNPs and removing INDELs. Variants were excluded using the following hard-filter thresholds: Quality by Depth (QD) \< 5, Mapping Quality (MQ) \< 45, Fisher Strand bias (FS) \> 60, Strand Odds Ratio (SOR) \> 3, Mapping Quality Rank Sum Test (MQRankSum) \< −2.5, Quality (QUAL) \< 140, or total site depth (DP) outside 900–1,800. Residual heterozygous calls were set to missing, as were individual genotypes supported by fewer than three reads. The dataset was then filtered to a minimum minor allele frequency (MAF) of 0.05 and a maximum per-site missingness of 30%. Finally, ten accessions from the Hachola site (Site 04) were removed for poor sample quality, yielding a final dataset of 290 accessions and 7,110,996 high-quality biallelic SNPs that served as input for all subsequent genomic analyses.

Field experiment and phenotyping

Three common-garden field experiments were conducted over three consecutive growing seasons (2019–2020, 2020–2021, 2021–2022) in a net-house at Mataim Farm, Hula Valley, Israel (33°09′08.4″N, 35°37′15.8″E). After dormancy-breaking germination, seedlings were transplanted into 5-L pots filled with a commercial substrate. In each season the 300 accessions were arranged in a randomized complete block design (one plant per accession per block): four blocks in each of the first two seasons (1,200 plants) and five blocks in the third (1,500 plants).

At the end of each season, five representative spikes per plant were harvested, oven-dried (35 °C, 30 days), threshed, and cleaned, and the grain was stored for nutritional analysis.

Determination of seed nutrients

Grain content of the four nutritional traits — protein, starch, β-glucan, and dietary fiber — was quantified non-destructively across all samples (\~3,900) with a near-infrared (NIR) analyzer (DA 7250, Perten Instruments). A wild-barley-adapted calibration was built from the manufacturer’s cultivated-barley package using laboratory reference measurements on 30 accessions (one randomly chosen per site): total starch and β-glucan by Megazyme enzymatic assay kits (K-TSTA and Mixed-Linkage β-Glucan; Megazyme, Ireland) and protein by the Kjeldahl method (Beljkaš et al., 2010), while fiber was retained from the cultivated package as no reference data were available. The single calibration was applied uniformly across all three seasons and validated against the laboratory references by Pearson correlation, with cross-season consistency confirmed by pairwise correlations between seasons.

Quantitative-genetic and environmental analyses

Phenotypic data preparation

All phenotypic data preparation, mixed-model fitting, and variance-component estimation were performed in R v4.4.3 (R Core Team, 2025). To match the 290-accession wild barley panel used for the genomic analyses, the three-season phenotypic dataset was filtered to remove the two cultivated control lines (Morex and Clipper) and all records of the ten Hachola (Site 04) accessions. Block 4 of the 2019–2020 season was additionally excluded from all downstream analyses, as its grain-composition values showed substantially lower inter-block correlations than Blocks 1–3. Ten traits were analyzed: the four nutritional traits (protein, starch, β-glucan, fiber) and six agromorphological traits (flowering time, number of tillers, plant height, grain weight, grain number, spike length). Each trait was mean-centered per season by subtracting the seasonal population mean from each observation, and trait-specific quality control was applied at each model-fitting step: observations missing for the trait under analysis were excluded, and values deviating more than three standard deviations from the trait mean were removed as outliers before fitting.

BLUP estimation

For each of the ten traits, genotype Best Linear Unbiased Predictions (BLUPs) were estimated from a linear mixed model fitted with the lme4 R package v1.1.37 (Bates et al., 2015), with genotype and season-by-block fitted as crossed random effects (the season-by-block term indexing each unique combination of growing season and block). Each genotype’s BLUP was taken as the sum of the model intercept and its genotype-specific random-effect deviation. These genotype BLUPs underlie the heritability, variance-partitioning, and trait-correlation analyses; the BLUPs of the four nutritional traits additionally served as the phenotypic input to the genome-wide association analysis.

Heritability and genotype × environment partitioning

Broad-sense heritability (H²) and variance partitioning were derived from a second linear mixed model that additionally included a genotype-by-season interaction as a random effect. The phenotypic variance was partitioned into genetic (V~G~), season-by-block (V~E~), genotype-by-season (V~G×E~), and residual (V~R~) components, extracted from the fitted model, and broad-sense heritability was computed as

H^2^ = V~G~ / (V~G~ + V~E~ + V~G×E~ + V~R~)

This interaction model provided the reported heritability values for the ten traits, the variance partitioning (for an eight-trait subset), and the per-trait reaction norms across the three seasons.

Trait variation, correlations, and environmental drivers

Geographic variation in the four nutritional traits was assessed by summarizing the centered genotype BLUPs by sampling site and ecotype (northern, coastal, desert, and the three ecotones).

Pairwise Pearson correlations were computed among the BLUPs of the ten traits, and p-values were corrected by the Benjamini–Hochberg false discovery rate (Benjamini & Hochberg, 1995) separately within the nutritional–nutritional and nutritional–morphological pair sets.

For each of the four nutritional traits, a linear mixed model was fitted with the trait BLUP as response, environmental predictors as fixed effects, and sampling site as a random intercept. Variables came from the provided 30-site environmental dataset; from ten a priori predictors, collinear variables (pairwise \|r\| \> 0.8) were removed, the remainder were standardized, and fixed effects were reduced by backward elimination (lmerTest). Standardized coefficients of the retained predictors were reported, with Benjamini–Hochberg correction across coefficients.

Genome-wide association and candidate-gene discovery

Association model

A genome-wide association study was performed for each of the four nutritional traits (protein, starch, β-glucan, dietary fiber) with EMMAX (Efficient Mixed-Model Association eXpedited; Kang et al., 2010), using the trait BLUPs as phenotypes and the full filtered set of 7,110,996 SNPs without further variant-level filtering. Population structure and relatedness were controlled using covariates computed once on an LD-pruned subset of 590,462 SNPs, obtained with PLINK v1.90b6.4 (Purcell et al., 2007) by pruning in 50-SNP windows advanced in 5-SNP steps at an r² threshold of 0.2: the first three principal components from a PLINK PCA were included as fixed covariates, and an aIBS (average identity-by-state) kinship matrix computed with emmax-kin was included as a random effect. Downstream analyses were performed in R v4.1.2 (R Core Team, 2021).

Significance and lead SNPs

Genome-wide significance was defined by a Bonferroni threshold at α = 0.10 over the 590,462 independent (LD-pruned) tests, corresponding to −log₁₀(P) = 6.7712. Fifteen SNPs exceeded this threshold; in addition, nine SNPs falling just below it were manually retained as marginal associations. Per-trait Manhattan and quantile–quantile plots were generated in R, and the genomic inflation factor (λ~GC~) was computed for each trait.

LD decay and candidate-gene extraction

Genome-wide LD decay was estimated across the 290 accessions from pairwise r² values computed in PLINK for all SNP pairs within 2.5 Mb on the seven chromosomes (1H–7H). Mean r² was binned at 1-kb resolution and smoothed with a pair-count-weighted LOESS curve (span 0.10), and the decay distance was taken where the smoothed curve crossed r² = 0.2. LD decayed to r² = 0.2 at \~188 kb, rounded to a ±200 kb window.

The significant and marginal SNPs were grouped into loci by single-linkage clustering at 200 kb (per chromosome and trait), with the most significant SNP of each locus taken as its lead. Each locus interval — the clustered span extended by 200 kb on each side — was intersected with the Morex V3 gene annotation (Ensembl Plants release 62; Yates et al., 2022) using bedtools intersect v2.30.0 (Quinlan & Hall, 2010) to obtain the candidate genes in the window. For each candidate gene, a gene-level VCF (gene coordinates ±1 kb) was extracted from the filtered dataset with bcftools v1.13 (Danecek et al., 2021), serving as input to the haplotype analysis.

Genotype imputation

For the haplotype analysis, missing genotypes were imputed with the sparse non-negative matrix factorization (sNMF) algorithm of the LEA R package v3.6.0 (Frichot & François, 2015). The dataset was converted to LFMM format, ten sNMF runs were performed at K = 3, and the run minimizing the cross-entropy criterion was retained; genotypes were then imputed from the mode of the ancestral genotype probabilities, yielding a fully homozygous matrix. The imputed data were used solely to estimate local linkage disequilibrium in the haplotype analysis.

Haplotype-based association testing

Haplotype-based association testing was performed for the candidate genes with the crosshap R package v1.4.0 (Marsh et al., 2023), across the four nutritional traits. For each gene, a square pairwise r² matrix was computed in PLINK from the imputed gene-level VCF, and haplotypes were constructed from the corresponding raw (non-imputed) gene-level VCF using that LD matrix. Haplotyping was run under two minimum marker-group sizes (MGmin = 2 and 3) and seven clustering thresholds (epsilon = 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 0.85), with a minimum haplotype size of nine accessions (minHap = 9); accessions in the unresolved haplotype-zero group were excluded from the statistical tests. For each gene × MGmin × epsilon combination, the association between haplotype grouping and trait value was tested with the Kruskal–Wallis rank-sum test (Kruskal & Wallis, 1952).

Multiple-testing correction

The up-to-fourteen tests per gene (two MGmin × seven epsilon values) were corrected in two stages. Within each gene, tests with identical results (same raw p-value, number of groups, and sorted group sizes) were collapsed to a single representative, the Holm step-down procedure (Holm, 1979) was applied across the remaining tests, and the smallest Holm-adjusted p-value was taken as the gene’s representative association. Then, separately within each trait, the Benjamini–Hochberg false discovery rate (Benjamini & Hochberg, 1995) was applied across genes to these representative p-values — with a Bonferroni correction computed in parallel for reference — and genes with FDR \< 0.05 were retained as significant.

Candidate-gene annotation and curation

The FDR-significant genes were functionally annotated by DIAMOND BLASTP (Buchfink et al., 2021) against UniProt Swiss-Prot, complemented by InterProScan domain assignments (Pfam/InterPro and GO; Jones et al., 2014) and, for genes without a confident hit, BLASTP against the NCBI non-redundant protein database (Altschul et al., 1997); the evidence was reconciled into a single functional call per gene. A manual biological-relevance filter then retained only genes with a plausible role in the pathway underlying the corresponding grain component, defining the final presented set.

Final figures

For each final gene, a single MGmin–epsilon combination from among the significant groupings was used to render a two-panel figure with ggplot2 v4.0.1 (Wickham, 2016), ggpubr v0.6.3 (Kassambara, 2023), and patchwork v1.3.2 (Pedersen, 2024): a violin plot of the per-haplotype-group trait-BLUP distribution (290 wild accessions; the haplotype-zero group excluded), annotated with the omnibus Kruskal–Wallis p-value and post-hoc pairwise Wilcoxon comparisons (Holm-corrected) against the largest group, above a genotype heatmap of one representative wild accession per group (reference, alternate, and missing genotypes across the gene’s SNPs). The final genes were assembled into a single composite figure.
