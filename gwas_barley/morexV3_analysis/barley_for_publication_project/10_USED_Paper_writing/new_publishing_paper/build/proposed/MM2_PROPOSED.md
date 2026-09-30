# Review item MM2: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM2 of the M&M voice review** (item 2 in the chat): standard procedures are named, not explained, as in Evgeny's papers. The same paragraphs also get the clarity fixes (B12, B13, A6 from the chat list) and the voice fixes ("we", aim first).
>
> **Removed:** the λ~GC~ formula · "crossed" · "the season-by-block term indexed…" · "BLUP = intercept + genotype effect" · the EMMAX build string ("beta; Intel 64-bit build of February 2012") · the LEA detail ("the sNMF algorithm", "ancestral populations", "assigned its most likely value") · repeated "of the 290 accessions". All of them stay in the hidden source notes.
>
> **Clarity fixes:** crosshap now gets one clause saying what it does (it groups SNPs by LD, then defines haplotypes), so "LD matrix" and "imputed vs observed genotypes" make sense on first reading. It no longer says crosshap "tested" the genes; the Kruskal–Wallis test does that. "Following Potapenko et al. (2026a)" now clearly refers to the choice of three PCs; before, the citation seemed to be for our PCs themselves.
>
> **Checked:** Mol Ecol uses the first three PCs and imputes with the LEA impute function at K = 3 (three genetic clusters). The crosshap wording matches the package code: dbscan on LD gives the marker groups, then the haplotypes are built from them. The LEA call is snmf K = 3, 10 runs, the lowest cross-entropy run, impute with the mode.
>
> Length: 4 paragraphs, 409 → 345 words. The crosshap clause adds about 25 words back, for clarity.
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** **EMMAX version.** The text now says only "EMMAX (Kang et al. 2010)", as Evgeny writes it. The binary prints no version, and the build (February 2012) is kept in the hidden note. OK to leave the version out of the text?
>

---

## Change 1 of 4 · M&M · Phenotypic data analysis ¶2 (BLUP model)

*Why:* Standard lme4 detail removed; aim first, 'we'.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Genotype best linear unbiased predictions (BLUPs) were estimated</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To estimate the genetic value of each accession, we fitted a linear mixed model</ins> for each of the ten traits <del style="color:#c0392b;background:#fdecea">with a linear mixed model fitted</del> in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as <del style="color:#c0392b;background:#fdecea">crossed </del>random effects<del style="color:#c0392b;background:#fdecea">; the season-by-block term indexed each block within each season. The BLUP of each accession was the model intercept plus its genotype random effect. These BLUPs were</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. The genotype best linear unbiased predictions (BLUPs) of this model were</ins> used as the phenotypes in all downstream analyses.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> To estimate the genetic value of each accession, we fitted a linear mixed model for each of the ten traits in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as random effects. The genotype best linear unbiased predictions (BLUPs) of this model were used as the phenotypes in all downstream analyses.

---

## Change 2 of 4 · M&M · Genome-wide association ¶1

*Why:* Aim first, 'we'; build string removed; the pruned set introduced in one sentence; the Potapenko citation placed on the choice it supports.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Genome-wide association was tested</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To identify genomic regions associated with the nutritional traits, we performed a genome-wide association study (GWAS)</ins> for each <del style="color:#c0392b;background:#fdecea">nutritional </del>trait with EMMAX (<del style="color:#c0392b;background:#fdecea">beta; Intel 64-bit build of February 2012; </del>Kang et al. 2010), using the genotype BLUPs as phenotypes and all 7,110,996 SNPs<del style="color:#c0392b;background:#fdecea"> of the 290 accessions</del>. <del style="color:#c0392b;background:#fdecea">A single linkage disequilibrium (LD)-pruned SNP set was used to correct for population structure and relatedness and to set the significance threshold. It was obtained</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Population structure, relatedness and the significance threshold were all derived from one linkage disequilibrium (LD)-pruned set of 111,017 SNPs, obtained</ins> with PLINK v1.90b6.4 (Purcell et al. 2007) with a 1,000-kb window, a step of one SNP and an r² threshold of 0.2<del style="color:#c0392b;background:#fdecea">, and comprised 111,017 SNPs</del>. <del style="color:#c0392b;background:#fdecea">The</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Following Potapenko et al. (2026a), the</ins> first three principal components of this set<del style="color:#c0392b;background:#fdecea">, computed in PLINK,</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> (computed in PLINK)</ins> were included as fixed covariates<del style="color:#c0392b;background:#fdecea"> (Potapenko et al. 2026a)</del>, and an identity-by-state kinship matrix <del style="color:#c0392b;background:#fdecea">computed </del>from the same set <del style="color:#c0392b;background:#fdecea">with emmax-kin was included</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">(emmax-kin)</ins> as a random effect.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> To identify genomic regions associated with the nutritional traits, we performed a genome-wide association study (GWAS) for each trait with EMMAX (Kang et al. 2010), using the genotype BLUPs as phenotypes and all 7,110,996 SNPs. Population structure, relatedness and the significance threshold were all derived from one linkage disequilibrium (LD)-pruned set of 111,017 SNPs, obtained with PLINK v1.90b6.4 (Purcell et al. 2007) with a 1,000-kb window, a step of one SNP and an r² threshold of 0.2. Following Potapenko et al. (2026a), the first three principal components of this set (computed in PLINK) were included as fixed covariates, and an identity-by-state kinship matrix from the same set (emmax-kin) as a random effect.

---

## Change 3 of 4 · M&M · Genome-wide association ¶2

*Why:* λGC formula removed (standard); threshold sentence untangled.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Genome-wide significance was</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The genome-wide significance threshold was</ins> set by a Bonferroni correction at α = 0.10<del style="color:#c0392b;background:#fdecea"> over the 111,017 LD-pruned SNPs, taken</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, with the 111,017 pruned SNPs taken</ins> as the number of independent tests (*P* < 9.008 × 10⁻⁷; −log~10~(*P*) = 6.0454). For each trait, <del style="color:#c0392b;background:#fdecea">Manhattan and quantile–quantile plots were drawn, and the genomic inflation factor (λ~GC~) was computed as the median of the χ² values (1 df) corresponding to the SNP *P* values, divided by the median of the χ² distribution with 1 df (0.455)</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">we drew Manhattan and quantile–quantile plots and calculated the genomic inflation factor (λ~GC~)</ins> (Fig. 3).<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins> Allelic effects (β) are given for the minor allele, in the units of the trait (%).

**Reads after approval:**

> The genome-wide significance threshold was set by a Bonferroni correction at α = 0.10, with the 111,017 pruned SNPs taken as the number of independent tests (*P* < 9.008 × 10⁻⁷; −log~10~(*P*) = 6.0454). For each trait, we drew Manhattan and quantile–quantile plots and calculated the genomic inflation factor (λ~GC~) (Fig. 3). Allelic effects (β) are given for the minor allele, in the units of the trait (%).

---

## Change 4 of 4 · M&M · Haplotype analysis ¶1

*Why:* Say what crosshap does before using its terms; crosshap defines haplotypes, it does not test; LEA detail shortened.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Candidate genes were tested for haplotype–trait associations with crosshap v1.4.0 (Marsh et al. 2023).</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To test haplotype–trait associations, we defined the haplotypes of each candidate gene with crosshap v1.4.0 (Marsh et al. 2023), which groups SNPs into marker groups by their LD and defines haplotypes by the marker-group alleles that each accession carries.</ins> For each gene, the SNPs within the gene and 1 kb on either side were extracted <del style="color:#c0392b;background:#fdecea">from the genotypes of the 290 accessions </del>with bcftools v1.13 (Danecek et al. 2021). <del style="color:#c0392b;background:#fdecea">Haplotypes were built from these observed genotypes, whereas the LD matrix used by crosshap to group SNPs was computed in PLINK (r²; Purcell et al. 2007) from imputed genotypes at the same positions.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The LD matrix (r², PLINK; Purcell et al. 2007) was computed from imputed genotypes, whereas haplotypes were assigned from the observed genotypes.</ins> Missing genotypes (18.2%<del style="color:#c0392b;background:#fdecea"> of all genotypes</del>) were imputed with <del style="color:#c0392b;background:#fdecea">the sNMF algorithm of </del>the LEA R package v3.6.0 (Frichot and François 2015)<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> from sNMF ancestry estimates</ins> with K = 3<del style="color:#c0392b;background:#fdecea"> ancestral populations</del>, the number of genetic clusters in this collection (Potapenko et al. 2026a)<del style="color:#c0392b;background:#fdecea">. Of ten runs, the run with the lowest cross-entropy was used, and each missing genotype was assigned its most likely value.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, using the run with the lowest cross-entropy of ten.</ins><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> To test haplotype–trait associations, we defined the haplotypes of each candidate gene with crosshap v1.4.0 (Marsh et al. 2023), which groups SNPs into marker groups by their LD and defines haplotypes by the marker-group alleles that each accession carries. For each gene, the SNPs within the gene and 1 kb on either side were extracted with bcftools v1.13 (Danecek et al. 2021). The LD matrix (r², PLINK; Purcell et al. 2007) was computed from imputed genotypes, whereas haplotypes were assigned from the observed genotypes. Missing genotypes (18.2%) were imputed with the LEA R package v3.6.0 (Frichot and François 2015) from sNMF ancestry estimates with K = 3, the number of genetic clusters in this collection (Potapenko et al. 2026a), using the run with the lowest cross-entropy of ten.

---

**Your decision:** approve · discard · or tell me what to change.
