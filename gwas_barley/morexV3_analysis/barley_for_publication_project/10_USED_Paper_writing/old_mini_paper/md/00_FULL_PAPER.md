**Genome-wide association and haplotype analysis identify candidate genes for grain nutritional quality in wild barley (*Hordeum vulgare* ssp. *spontaneum*)**

Student: Shahar Liviatan

Supervisor: Sariel Hübner

# Abstract

Wild barley (*Hordeum vulgare* ssp. *spontaneum*), the progenitor of cultivated barley, retains broad genetic diversity that has been reduced in cultivated germplasm through domestication and intensive breeding. Identifying the genes underlying its nutritional variation can provide precise targets for improving grain quality in cultivated barley while minimizing the transfer of undesirable linked variation.

Here, 290 wild-barley accessions sampled across Israel's environmental gradients were phenotyped over three growing seasons for four grain nutritional traits: protein, starch, β-glucan, and dietary fiber. Mixed-model analyses revealed trait-specific genetic architecture: starch showed the highest genetic contribution, whereas protein, β-glucan, and fiber had lower broad-sense heritability and larger residual components; β-glucan also showed the largest genotype-by-environment component. Starch and fiber displayed opposing geographic, environmental, and phenotypic patterns, suggesting a carbon-allocation trade-off between storage starch and cell-wall components. Genome-wide association mapping identified 20 loci across the four traits. Haplotype-based analysis refined these to four candidate genes: an α-glucan phosphorylase and a phosphate/anion transporter for starch, an AP2/ERF transcription factor for β-glucan, and a BAHD acyltransferase for fiber; protein yielded no significant candidate genes. Wild barley harbors genetically tractable variation in grain composition and provides candidate loci and haplotypes for future evaluation in cultivated germplasm.

**Keywords:** wild barley; crop wild relatives; grain nutritional quality; GWAS;

# Introduction

Climate change and population growth are intensifying pressure on the global food supply, making food security an urgent priority (Godfray et al., 2010).

For decades, crop breeding has targeted yield, often at the expense of grain nutritional quality and of the diversity available for further improvement (Dempewolf et al., 2017). Higher grain yield is also frequently associated with lower grain protein concentration, partly because carbohydrate and dry-matter accumulation can outpace nitrogen accumulation during grain filling (Simmonds, 1995; Bogard et al., 2010). Dietary risks remain major contributors to the global burden of disease (GBD 2019 Risk Factors Collaborators, 2020). Improving the nutritional quality of cereal grains is therefore an important breeding target.

Crop wild relatives (CWR) are a rich reservoir of beneficial alleles to improve disease resistance, abiotic-stress tolerance, and yield in major crops, yet their potential for nutritional value remains largely untapped (Dempewolf et al., 2017). A principal obstacle is genetic drag, whereby an introgressed wild allele carries linked deleterious variants that degrade agronomic performance (Hübner and Kantar 2021; Huang et al. 2023). High-resolution mapping mitigates this by pinpointing quantitative trait loci (QTLs) at fine resolution to be introgressed while minimizing genetic drag (Hübner & Kantar, 2021), thus identification of beneficial alleles at high-resolution is essential for exploiting CWR in breeding.

Barley (*Hordeum vulgare*) was domesticated in the Fertile Crescent 10,000 years ago and is now the fourth most important cereal worldwide (Badr et al., 2000; FAO 2025). Its wild progenitor, *H. vulgare* ssp. *spontaneum* is distributed along a wide range of environments in the Fertile Crescent and adjacent regions (Harlan & Zohary, 1966). Wild barley is characterized with a broad genetic and phenotypic variation and attempts to use it in breeding of disease, drought, and salt tolerance has been successful (Nevo & Chen, 2010). The grain is dominated by starch, with protein, dietary fiber, and the soluble fiber β-glucan as the other major components (Farag et al., 2022). Crucially, wild barley is nutritionally superior to cultivated barley, with on average \~50% more protein, 38% less starch, and 65% more fiber (Friedman & Atsmon, 1988), indicating a potential genetic resource for improving grain quality in cultivated barley.

Advances in genomics now allow to identify the genetic basis of target traits at high resolution. Genome-wide association studies (GWAS) is widely used to link genetic to phenotypic variation and have identified trait-associated loci and candidate genes in barley, and many other crops (Alqudah et al., 2020). To dissect the genetic basis of grain nutritional traits in wild barley, I analyzed a collection of wild barley sampled across a wide range of environments and characterized genetically and phenotypically, thus allowing to identify candidate genes for key nutritional traits.

## Hypothesis and objectives

The nutritional value of wild barley grain is attributed to ecological constraints, thus heritable genetic variation is expected. Screening wild-barley populations that represent contrasting ecological niches can highlight the genetic factors underlying these traits. To test this hypothesis, I will address the following objectives: (i) characterize the variation in nutritional traits among wild-barley populations representing different ecotypes in Israel; (ii) examine the correlation between grain nutritional value traits, environmental gradients, and plant-fitness traits; and (iii) identify QTLs and candidate genes associated with nutritional traits.

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

# Results

## Characterization of wild barley grain nutritional value traits

The four nutritional traits were quantified in 290 *H. spontaneum* accessions across three growing seasons, and for each trait genotype BLUPs were used to estimate accession-level genotypic values while accounting for season and block effects. Partitioning the phenotypic variance into genetic, season-and-block, genotype-by-season (G×E), and residual components (Figure 1A) showed that starch had the highest broad-sense heritability (H² = 0.474; Table 1), roughly twice as all other traits. The genotype-by-season interaction was small to moderate for all four traits (4.5% - 15.2%), where β-glucan combined the largest G×E component (15.2%) with the smallest genetic fraction indicating a strong genotype-specific response across the three growing seasons. Among morphological traits flowering time (H² = 0.747) had the highest heritability and tiller number the lowest (H^2^ = 0.031) indicating the respective genetic control for these traits. To further explore these trends, reaction norms were generated per genotype across seasons (Figure 1B). Overall, results were consistent with the heritability with substantial variation between specific genotypes indicating that the genotype X environment signature is notable for all traits.

  ---------------------------------------
  **Trait**               **H²**
  ----------------------- ---------------
  *Starch*                0.474

  *Fiber*                 0.248

  *Protein*               0.242

  *β-glucan*              0.234

  *Flowering time*        0.747

  *Spike length*          0.344

  *Grain weight*          0.223

  *Tillers*               0.031
  ---------------------------------------

![](media_00_FULL_PAPER/media/image1.png){width="6.427777777777778in" height="2.5708333333333333in"}**\
Table 1.** Broad-sense heritability (H²) of grain nutritional and morphological traits, estimated from the genotype × environment model as the genetic fraction of total phenotypic variance, for the four nutritional and four morphological traits, ordered by decreasing H².

**Figure 1.** Genetic architecture of grain nutritional and reference traits. (A) Variance partitioning for four nutritional and four morphological traits into genetic, season + block, genotype × environment, and residual components. (B) Per-genotype reaction norms across the three seasons for the four nutritional traits; grey lines, all accessions; colored lines, four representative trajectories.

To examine whether variation in grain nutritional traits followed the environmental gradient of the collection sites, trait BLUPs were first visualized across the 29 sites ordered by ecological region (Figure 2A). This descriptive pattern suggested that northern-derived genotypes generally had higher starch, whereas desert-derived genotypes tended to have lower starch and higher fiber and β-glucan; protein showed no clear directional trend. Among the nutritional traits, starch was strongly negatively correlated with fiber (r = −0.78; Figure 2B) and β-glucan (r = −0.53), while fiber and β-glucan were positively correlated (r = 0.35; all P \< 0.001). Thus, genotypes with higher starch tended to have lower fiber and β-glucan. Protein was poorly correlated with other traits, and significant negative correlation was obtained only for protein-starch (r = -0.12, P \< 0.05). Among the morphological and phenological traits, starch was strongly correlated with grain weight, plant height, and flowering time, while fiber showed the exact opposite pattern. Protein was also weakly negatively correlated with grain weight (r = −0.13, P \< 0.05). Overall, β-glucan and protein displayed weaker, mostly negative correlations, linking grain composition directly to whole-plant resource allocation and phenology (Figure 2C). To identify which environmental variables of the sites of origin were associated with each trait, site-mean nutritional-trait BLUPs were correlated with eight environmental variables (Figure 2D). Starch was positively correlated with clay, precipitation, and organic carbon, and negatively correlated with sand and pH. Fiber showed the opposite pattern for clay, sand, and precipitation, while β-glucan was positively correlated with sand and electrical conductivity and negatively correlated with clay, precipitation, and organic carbon. Protein was not significantly correlated with any of the measured environmental variables. Together, these associations linked the contrasting patterns of starch versus fiber and β-glucan to gradients in soil texture and precipitation, consistent with the broader geographic pattern observed across the collection sites.

## Identification of nutritional QTLs and candidate genes 

To identify the genomic regions that are associated with nutritional traits, we used the calculated BLUPs as the phenotypic input in the model, the first three principal components as fixed covariates to correct for population structure, and a kinship matrix as random effect to correct for relatedness (Figure 3B). The genomic inflation factor was calculated for each model indicating adequate control of genome-wide false positive inflation (λ~GC~ \> 0.9; Figure 3A).

Twenty-four SNPs (15 significant, 9 marginal) resolved into 20 independent QTLs. Interestingly, highly heritable starch (H² = 0.4[74]{dir="rtl"}) yielded no significant signals, while mildly heritable β-glucan (H² = 0.234) proved most mappable, spanning 11 loci. These included the strongest peak on chromosome 2H (2H:41,996,626) and cluster collapsing four significant SNPs on chromosome 4H (4H:34,825,152 - 34,971,609). Fiber was

![](media_00_FULL_PAPER/media/image2.png){width="5.954166666666667in" height="5.305555555555555in"}**Figure 2. Ecological architecture of grain nutritional traits.** (A) Centered trait BLUPs of the four nutritional traits across 29 sampling sites, ordered and coloured by ecological region. (B) Pearson correlations among nutritional-trait BLUPs. (C) Pearson correlations between nutritional- and morphological-trait BLUPs. (D) Pearson correlations between site-mean nutritional-trait BLUPs and eight environmental variables; cells outlined in black are significant at P \< 0.05. Significance (B--D): \*P \< 0.05, \*\*P \< 0.01, \*\*\*P \< 0.001.

associated with four QTLs (1H:344,520,079; 4H:24,707,376; 7H:14,817,657; 7H:573,606,306), whereas protein yielded a single locus on chromosome 3H (3H:106,623,911), and starch produced four marginal loci (3H:546,433,616; 6H:525,776,080; 7H:151,110,354; 7H:573,606,460). Notably, on the distal arm of chromosome 7H, the significant fiber SNP (chr7H:573,606,306) and a marginal starch SNP (Chr7H:573,606,460) co-localized within 154 bp, defining a nine-gene window consistent with shared or linked genetic control driving the negative starch--fiber relationship (Figure 3).

To search for potential candidate genes within the identified QTLs, a region of 200Kb was defined up and downstream of significant SNPs in accordance with the profile of LD decay (r² \< 0.2; Figure 3C). Overall 108 candidate genes were identified within the 20 loci: 72 for β-glucan, 18 for fiber, 17 for starch, and one for protein.

![](media_00_FULL_PAPER/media/image3.png){width="6.034722222222222in" height="4.186805555555556in"}**Figure 3.** Association mapping, structure correction, and LD diagnostics. (A) Manhattan (left) and QQ (right) plots for the four traits (290 accessions; 7,110,996 SNPs); the line marks the Bonferroni threshold (−log₁₀P = 6.7712), λ~GC~ the genomic inflation factor. (B) PC scree plot (PC1 omitted); bars, % genomic variance per component; line, cumulative %. (C) Genome-wide LD decay (1H--7H); grey points, mean r² in 1-kb bins (0--2.5 Mb); red curve, LOESS fit; dashed line, r² = 0.2, crossed at \~188 kb.

## Haplotype analysis of candidate genes for grain nutritional traits

Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accession carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. We retained 45 of the 108 candidates genes for biological screening where significant phenotypic difference was observed between contrasting alleles (FDR \< 0.05). The 45 genes were functionally annotated and a final biological-relevance screen kept only genes with a documented or plausible role in the biosynthesis, transport, or regulation of the corresponding grain component. This process narrowed the genome-wide signal of the four nutritional traits to five biologically plausible candidate genes (Figure 4). The fifth, a 2H pectin methylesterase associated with β-glucan showed only a small haplotype effect and is therefore not presented as a separate panel. Within this shared 7H window, three genes passed the trait-wise FDR filter for both starch and fiber: a probable high-affinity nitrate transporter (7HG0729020), a pentatricopeptide-repeat protein (7HG0729090), and a chloroplastic NifU-like protein (7HG0729100). For each gene, the starch and fiber analyses yielded identical genotype-defined haplotype assignments, and in all three cases the haplotype group with the highest mean fiber had the lowest mean starch. This reciprocal pattern provides locus-level support for shared or tightly linked genetic control of the starch--fiber relationship.

The two starch candidates lie within a single locus on the long arm of 3H (3H:546,433,616), the α-glucan phosphorylase *Pho* (3HG0301750) resolved into three haplotype groups that differed strongly in grain starch (Kruskal-Wallis P = 6.8 × 10⁻⁹; ε² = 0.18), with one small group carrying distinctly lower starch than the two larger groups. Pho encodes a plastidial α-1,4-glucan phosphorylase that reversibly transfers glucosyl units between glucose-1-phosphate and α-1,4-glucan chains. In developing barley endosperm, HvPho1 is active from the onset of grain development and can generate linear glucans proposed to serve as primers for amylopectin and starch biosynthesis []{dir="rtl"}(Cuesta-Seijo et al., 2017). The phosphate/anion transporter *PHT* (3HG0301710) resolved into four groups and showed the strongest haplotype separation among the genes examined (P = 1.7 × 10⁻¹⁰; ε² = 0.23), with additive effect between haplotypes. PHT4 transporters mediate inorganic-phosphate transport across plastid membranes and maintain plastid phosphate homeostasis, which is closely coupled to carbon partitioning and starch synthesis. Consistent with this role, disruption of the related plastidial transporter PHT4;2 markedly reduces starch accumulation, providing a functional link between PHT4-mediated phosphate transport and starch content (Guo et al., 2008; Irigoyen et al., 2011). Although the shared lead SNP was only marginal at the genome-wide level, the gene-level Kruskal--Wallis associations for both candidates remained significant after trait-wise Bonferroni correction across the tested starch genes.

For the β-glucan and fiber candidates smaller haplotype effects were noticed. Their gene-level associations remained significant across the tested genes. The β-glucan candidate, the AP2/ERF transcription factor (3HG0299440) resolved into two haplotype groups that differed modestly in grain β-glucan (Kruskal--Wallis P = 7.8 × 10⁻³; ε² = 0.03), with the smaller group showing higher β-glucan values than the larger group. AP2/ERF-binding motifs are over-represented in the proximal promoter of HvCslF6, a key synthase of grain mixed-linkage β-glucan, linking this candidate to the transcriptional regulation of β-glucan accumulation (Garcia-Gimenez et al., 2022). []{dir="rtl"}The fiber candidate, the BAHD acyltransferase (7HG0642350) resolved into five haplotype groups with a similarly modest separation in grain fiber (P = 4.0 × 10⁻³; ε² = 0.06); among the post-hoc comparisons against the largest haplotype group (A), only the A-versus-B contrast was significant, with group A showing higher fiber values than group B (Holm-adjusted P = 0.040). BAHD acyltransferases modify grass cell walls by transferring hydroxycinnamoyl groups to arabinoxylan, a major fiber polysaccharide, thereby affecting wall structure and digestibility. Experimental suppression of a BAHD acyltransferase reduced arabinoxylan-bound p-coumarate and increased cell-wall digestibility, linking this enzyme family to fiber architecture (Mota et al., 2021).

![](media_00_FULL_PAPER/media/image4.png){width="5.29375in" height="4.854861111111111in"}

**Figure 4.** Haplotype structure of the four final candidate genes. Violin/boxplots show trait-BLUP distributions by haplotype group, with genotype strips below: (a) Pho and (b) PHT for starch, (c) AP2/ERF for β-glucan, and (d) BAHD for fiber. The omnibus Kruskal--Wallis P value is shown above each panel. Brackets indicate post-hoc Wilcoxon comparisons of each haplotype group against the largest group, with P values adjusted by Holm independently of the within-gene Holm correction used to select the representative crosshap solution (\* P ≤ 0.05, \*\* P ≤ 0.01, \*\*\* P ≤ 0.001, \*\*\*\* P ≤ 10⁻⁴; ns, not significant).

# Discussion

## Ecological differentiation of grain composition

Previous studies of Israeli wild barley showed that population structure, flowering time, and growth traits are organized along gradients of temperature and precipitation (Hübner et al., 2009, 2013). The present study extends this ecological differentiation to grain nutritional traits. Because the comparisons were based on genotype BLUPs from common-garden experiments, the regional patterns reflect heritable differences among source populations rather than immediate plastic responses to the environments in which the accessions were collected. Northern genotypes tended to produce starch-rich grains, whereas desert genotypes produced less starch and more fiber and β-glucan (Figure 2). Starch and fiber showed opposing associations with sand, clay, and precipitation, while β-glucan generally followed the fiber pattern, linking the carbohydrate axis to environmental variation among the source sites (Figure 2D).

Physiological studies provide a plausible explanation for the direction of this pattern. Heat and drought during grain filling shorten the filling period, reduce the activity of starch-synthetic enzymes, and decrease starch deposition and final grain weight (Savin & Nicolas, 1996). In barley, soluble starch synthase is sensitive to high temperature, and post-anthesis heat or drought alters starch accumulation and endosperm development (Savin & Nicolas, 1996). Moreover, desert genotypes had higher β-glucan, but with strong G×E interaction. Previous studies in barley have shown that temperature can affect β-glucan content, solubility, viscosity, and molecular weight, with the direction and magnitude depending on genotype, stress intensity, and timing (Anker-Nilssen et al., 2008). Higher growth temperatures can increase soluble β-glucan in some genotype--environment combinations, whereas severe short-term heat stress reduced β-glucan and promoted its degradation (Anker-Nilssen et al., 2008). Thus our results on the large G×E corroborate with previous studies. Taken []{dir="rtl"}together, the nutritional traits were not distributed randomly across the genotype collection but were tightly associated with environmental conditions.

## A shared carbon-allocation axis for starch, fiber, and β-glucan

Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, consistent with previous evidence that grain protein concentration reflects the balance between carbon deposition and nitrogen acquisition, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and cell-wall polysaccharides (fiber and β-glucan) (Corke et al., 1989). This trade-off is consistently observed across genotypes, a geographic north to desert gradient (Figure 2), and a shared 7H locus, indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm specific overexpression of HvCslF6 increases β-glucan accumulation while reducing starch content (Lim et al., 2020), demonstrating that these nutritional components are governed by interconnected physiological and genetic systems. This metabolic axis integrates with a broader whole plant resource allocation strategy linked to plant architecture and phenology. Larger, later flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster developing plants distribute resources toward higher grain numbers and structural cell wall material. While domestication and breeding have historically shifted cultivated barley toward the starch rich end of this axis (i.e. higher starch and lower protein and fiber varieties), wild barley accessions remain a potential source of genetic variation for optimizing the fiber and β-glucan rich properties of the crop (Friedman & Atsmon, 1988). In contrast, no significant correlations were observed between protein and the measured environmental variables (Figure 2D). Additionally, a weak negative correlation between protein and grain weight was observed (Figure 2C). This pattern suggests that grain protein is governed by processes that differ from those organizing carbohydrate composition. This observation is largely supported in the literature highlighting a negative relationship between grain yield or grain weight and grain protein concentration (Simmonds, 1995). This negative relationship can potentially be explained by a dilution effect, where carbohydrates accumulate faster than nitrogen during grain filling, causing protein to become a smaller proportion of total grain mass (Bogard et al., 2010). Differences between wild and cultivated barley in dry-matter accumulation, nitrogen allocation, and nitrogen harvest index further demonstrated that grain protein depends on the balance between carbon deposition and nitrogen acquisition rather than on direct competition with starch alone (Corke et al., 1989). Total soil nitrogen and organic carbon were also tested, but neither was significantly correlated with grain protein (Figure 2D). This separation has practical implementation as strong starch--fiber trade-off may not impose the same constraint on protein, and increasing fiber or β-glucan need not cause a reduction in protein.

## Heritability and GWAS reveal different aspects of genetic architecture

Trait heritability poorly predicted genome-wide association mapping power. Despite exhibiting the highest broad-sense heritability, starch produced a sparse mapping profile, whereas β-glucan had lower broad-sense heritability than starch and higher environmental sensitivity, yet yielded a robust, distributed signal comprising 11 significant loci across all seven chromosomes (Figure 3A). This contrast emphasizes that mapping success depends on complex genetic architecture components like number of alleles and their frequencies, effect sizes, and linkage disequilibrium rather than on total genetic variation alone. For example, the starch locus on chromosome 3H illustrates how marginal signal in the GWAS analysis produced a strong separation between haplotypes once the entire haplotype which groups several markers together was analyzed.

The linked starch candidates on 3H offer complementary mechanisms for a single locus. The glucan phosphorylase (Pho) directly participates in []{dir="rtl"}metabolism and endosperm glucan synthesis (Cuesta-Seijo et al., 2017). The phosphate transporter (PHT) suggests that phosphate homeostasis contributes to the metabolic conditions required for starch accumulation, consistent with the impaired grain filling reported in rice phosphate-transporter mutants (Ma et al., 2021). []{dir="rtl"}For β-glucan, the AP2/ERF candidate is consistent with transcriptional regulation of HvCslF6. AP2/ERF-binding motifs are over-represented in the proximal promoter of HvCslF6, a key synthase of grain mixed-linkage β-glucan, placing this candidate within its broader regulatory network (Garcia-Gimenez et al., 2022). Additionally, the BAHD candidate highlights a mechanism for fiber organization. BAHD enzymes mediate cell wall feruloylation and polymer cross linking, directly impacting structural rigidity and digestibility (de Souza et al., 2018). These four genes therefore identify distinct biological pathways for future investigation: direct starch metabolism, phosphate associated grain filling, glucan transcriptional regulation, and cell wall fiber structural modification.

# Conclusions

Grain composition in wild barley is organized around carbon allocation involving starch, fiber, and β glucan, alongside a distinct protein axis governed by carbon and nitrogen balance. By applying genome wide association studies and haplotype analysis to ecologically diverse populations, this research identified heritable variation and specific genetic architectures for these traits. Trait associations varied significantly; β glucan was distributed across multiple loci, whereas starch required haplotype level resolution to identify. Key regions, including a shared 7H locus and the 3H Pho/PHT locus, successfully link grain composition to glucan metabolism, phosphate transport, transcriptional regulation, and cell wall organization.

These defined loci provide precise targets for improving cultivated barley nutrition while minimizing the linkage drag commonly associated with wild germplasm introgression. Previous research successfully transferred beneficial wild alleles to increase β amylase activity and yield, yet the simultaneous transfer of undesirable traits remains a persistent challenge. Refining broad genetic associations into a focused set of biologically plausible targets simplifies the monitoring of transferred genomic segments. Future breeding efforts must validate these haplotypes in elite cultivars, replicate their effects in independent populations, and conduct functional tests to assess their overall agronomic impact.

# References

Alqudah AM, Sallam A, Baenziger PS, Borner A (2020) GWAS: fast-forwarding gene identification and characterization in temperate cereals - lessons from barley, a review. Journal of Advanced Research 22:119-135.

Anker-Nilssen K, Sahlstrom S, Knutsen SH, Holtekjolen AK, Uhlen AK (2008) Influence of growth temperature on content, viscosity and relative molecular weight of water-soluble beta-glucans in barley (Hordeum vulgare L.). Journal of Cereal Science 48:670-677.

Badr A, Muller K, Schafer-Pregl R, El Rabey H, Effgen S, Ibrahim HH, Pozzi C, Rohde W, Salamini F (2000) On the origin and domestication history of barley (Hordeum vulgare). Molecular Biology and Evolution 17:499-510.

Bogard M, Allard V, Brancourt-Hulmel M, Heumez E, Machet JM, Jeuffroy MH, Gate P, Martre P, Le Gouis J (2010) Deviation from the grain protein concentration-grain yield negative relationship is highly correlated to post-anthesis N uptake in winter wheat. Journal of Experimental Botany 61:4303-4312.

Corke H, Avivi N, Atsmon D (1989) Pre- and post-anthesis accumulation of dry matter and nitrogen in wild barley (Hordeum spontaneum) and in barley cultivars (H. vulgare) differing in final grain size and protein content. Euphytica 40:127-134.

Cuesta-Seijo JA, Ruzanski C, Krucewicz K, Meier S, Hagglund P, Svensson B, Palcic MM (2017) Functional and structural characterization of plastidic starch phosphorylase during barley endosperm development. PLoS ONE 12:e0175488.

Dempewolf H, Baute G, Anderson J, Kilian B, Smith C, Guarino L (2017) Past and future use of wild relatives in crop breeding. Crop Science 57:1070-1082.

de Souza WR, Martins PK, Freeman J, Pellny TK, Michaelson LV, et al. (2018) Suppression of a single BAHD gene in Setaria viridis causes large, stable decreases in cell wall feruloylation and increases biomass digestibility. New Phytologist 218:81-93.

FAO (2025) FAOSTAT statistical database. Food and Agriculture Organization of the United Nations, Rome.

Farag MA, Xiao J, Abdallah HM (2022) Nutritional value of barley cereal and better opportunities for its processing as a value-added food: a comprehensive review. Critical Reviews in Food Science and Nutrition 62:1092-1104.

Friedman M, Atsmon D (1988) Comparison of grain composition and nutritional quality in wild barley (Hordeum spontaneum) and in a standard cultivar. Journal of Agricultural and Food Chemistry 36:1167-1172.

Garcia-Gimenez G, et al. (2022) Identification of candidate MYB transcription factors that influence CslF6 expression in barley grain. Frontiers in Plant Science 13:883139.

GBD 2019 Risk Factors Collaborators (2020) Global burden of 87 risk factors in 204 countries and territories, 1990-2019: a systematic analysis for the Global Burden of Disease Study 2019. The Lancet 396:1223-1249.

Godfray HCJ, Beddington JR, Crute IR, Haddad L, Lawrence D, Muir JF, Pretty J, Robinson S, Thomas SM, Toulmin C (2010) Food security: the challenge of feeding 9 billion people. Science 327:812-818.

Guo B, Jin Y, Wussler C, Blancaflor EB, Motes CM, Versaw WK (2008) Functional analysis of the Arabidopsis PHT4 family of intracellular phosphate transporters. New Phytologist 177:889-898.

Harlan JR, Zohary D (1966) Distribution of wild wheats and barley. Science 153:1074-1080.

Holland JB, Nyquist WE, Cervantes-Martinez CT (2003) Estimating and interpreting heritability for plant breeding: an update. Plant Breeding Reviews 22:9-112.

Huang K, et al. (2023) The genomics of linkage drag in inbred lines of sunflower. Proceedings of the National Academy of Sciences of the United States of America 120:e2205783119.

Hubner S, Hoffken M, Oren E, Haseneyer G, Stein N, Graner A, Schmid K, Fridman E (2009) Strong correlation of wild barley (Hordeum spontaneum) population structure with temperature and precipitation variation. Molecular Ecology 18:1523-1536.

Hubner S, Bdolach E, Ein-Gedi S, Schmid KJ, Korol A, Fridman E (2013) Phenotypic landscapes: phenological patterns in wild and cultivated barley. Journal of Evolutionary Biology 26:163-174.

Hubner S, Kantar MB (2021) Tapping diversity from the wild: from sampling to implementation. Frontiers in Plant Science 12:626565.

Irigoyen S, Karlsson PM, Kuruvilla J, Spetea C, Versaw WK (2011) The sink-specific plastidic phosphate transporter PHT4;2 influences starch accumulation and leaf size in Arabidopsis. Plant Physiology 157:1765-1777.

Kang HM, Sul JH, Service SK, Zaitlen NA, Kong SY, Freimer NB, Sabatti C, Eskin E (2010) Variance component model to account for sample structure in genome-wide association studies. Nature Genetics 42:348-354.

Lim WL, Collins HM, Byrt CS, Lahnstein J, Shirley NJ, Aubert MK, Tucker MR, Peukert M, Matros A, Burton RA (2020) Overexpression of HvCslF6 in barley grain alters carbohydrate partitioning plus transfer tissue and endosperm development. Journal of Experimental Botany 71:138-153.

Ma B, Zhang L, Gao Q, et al. (2021) A plasma membrane transporter coordinates phosphate reallocation and grain filling in cereals. Nature Genetics 53:906-915.

Marsh JI, Petereit J, Johnston BA, Bayer PE, Tay Fernandez CG, Al-Mamun HA, Batley J, Edwards D (2023) crosshap: R package for local haplotype visualization for trait association analysis. Bioinformatics 39:btad518.

Mota TR, de Souza WR, Oliveira DM, Martins PK, Sampaio BL, et al. (2021) Suppression of a BAHD acyltransferase decreases p-coumaroyl on arabinoxylan and improves biomass digestibility in the model grass Setaria viridis. The Plant Journal 105:136-150.

Nevo E, Chen G (2010) Drought and salt tolerances in wild relatives for wheat and barley improvement. Plant, Cell and Environment 33:670-685.

Pintel N (2022) The Genomic Basis of Nutritional Traits in Wild Barley. MSc thesis, Tel-Hai College, Israel.

Potapenko E, Shermeister B, Mandel T, Hubner S (2026) Intraspecific genome size variation is attributed to adaptive silencing of transposable elements in Hordeum species. Molecular Biology and Evolution 43:msag051.

Savin R, Nicolas ME (1996) Effects of short periods of drought and high temperature on grain growth and starch accumulation of two malting barley cultivars. Australian Journal of Plant Physiology 23:201-210.

Simmonds NW (1995) The relation between yield and protein in cereal grain. Journal of the Science of Food and Agriculture 67:309-315.

Yates AD, Allen J, Amode RM, Azov AG, Barba M, Becerra A, Bhai J, Campbell LI, et al. (2022) Ensembl Genomes 2022: an expanding genome resource for non-vertebrates. Nucleic Acids Research 50:D996-D1003.
