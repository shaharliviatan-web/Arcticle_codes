# Review item S3: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item S3 — renumbering Figs. 3, 4, 5 → 2, 3, 4** (old Fig. 2 is dissolved). References and captions only; applied together with S4 and S5 so the file never has two "Fig. 2" captions.
>
> Your Word comment 311 ("See figure 3 start of chrom 4") also names the Manhattan figure, so its number changes; nothing else in it.
>

## Change 1 of 25 · M&M · Genome-wide association ¶2

*Why:* Renumbering: old Fig. 3 is now Fig. 2.

**With the changes marked:**

The genome-wide significance threshold was set by a Bonferroni correction at α = 0.10, with the 111,017 pruned SNPs taken as the number of independent tests (*P* < 9.008 × 10⁻⁷; −log~10~(*P*) = 6.0454). For each trait, we drew Manhattan and quantile–quantile plots and calculated the genomic inflation factor (λ~GC~) (Fig. <del style="color:#c0392b;background:#fdecea">3</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2</ins>). Allelic effects (β) are given for the minor allele, in the units of the trait (%).

**Reads after approval:**

> The genome-wide significance threshold was set by a Bonferroni correction at α = 0.10, with the 111,017 pruned SNPs taken as the number of independent tests (*P* < 9.008 × 10⁻⁷; −log~10~(*P*) = 6.0454). For each trait, we drew Manhattan and quantile–quantile plots and calculated the genomic inflation factor (λ~GC~) (Fig. 2). Allelic effects (β) are given for the minor allele, in the units of the trait (%).

---

## Change 2 of 25 · M&M · Haplotype analysis ¶3

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource N)<sup>(TODO 224)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014)<sup>(TODO 225)</sup>, and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. <del style="color:#c0392b;background:#fdecea">4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3</ins>).

**Reads after approval:**

> For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource N)<sup>(TODO 224)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014)<sup>(TODO 225)</sup>, and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 3).

---

## Change 3 of 25 · M&M · Functional annotation ¶3 (7H gene)

*Why:* Renumbering: old Fig. 5 is now Fig. 4.

**With the changes marked:**

To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes (Online Resource N)<sup>(TODO 229)</sup>. As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. <del style="color:#c0392b;background:#fdecea">5</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4</ins>).

**Reads after approval:**

> To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes (Online Resource N)<sup>(TODO 229)</sup>. As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. 4).

---

## Change 4 of 25 · M&M · Comparison with elite cultivars ¶2

*Why:* Renumbering: old Fig. 4 and 5 is now Fig. 3 and 4.

**With the changes marked:**

For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. <del style="color:#c0392b;background:#fdecea">4d–f</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3d–f</ins>). The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups. For the 7H gene (Fig. 5c), the same five cultivars and procedure were used without a new call-rate screen.

**Reads after approval:**

> For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. 3d–f). The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups. For the 7H gene (Fig. 5c), the same five cultivars and procedure were used without a new call-rate screen.

---

## Change 5 of 25 · M&M · Comparison with elite cultivars ¶2 (7H)

*Why:* Renumbering: old Fig. 5 is now Fig. 4.

**With the changes marked:**

For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. 4d–f). The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups. For the 7H gene (Fig. <del style="color:#c0392b;background:#fdecea">5c</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4c</ins>), the same five cultivars and procedure were used without a new call-rate screen.

**Reads after approval:**

> For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. 4d–f). The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups. For the 7H gene (Fig. 4c), the same five cultivars and procedure were used without a new call-rate screen.

---

## Change 6 of 25 · M&M · Geographic origin

*Why:* Renumbering: old Fig. 4, 5 is now Fig. 3, 4.

**With the changes marked:**

To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (Figs. <del style="color:#c0392b;background:#fdecea">4, 5</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3, 4</ins>) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure (Online Resource N)<sup>(TODO 237)</sup>. To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

**Reads after approval:**

> To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (Figs. 3, 4) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure (Online Resource N)<sup>(TODO 237)</sup>. To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

---

## Change 7 of 25 · Results Ch. 2 ¶1 (1/2)

*Why:* Renumbering: old Fig. 3 is now Fig. 2.

**With the changes marked:**

To identify the genomic regions associated with grain nutritional traits, genome-wide association was performed for each trait using the genotype BLUPs as the phenotypic input, the first three principal components as fixed covariates to correct for population structure, and a kinship matrix as a random effect to correct for relatedness among the 290 accessions, testing 7,110,996 SNPs (Fig. <del style="color:#c0392b;background:#fdecea">3</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2</ins>). The genomic inflation factor was close to unity for all four traits (λ~GC~ = 0.975–1.013; Fig. 3), indicating adequate control of population structure and of genome-wide false-positive inflation.

**Reads after approval:**

> To identify the genomic regions associated with grain nutritional traits, genome-wide association was performed for each trait using the genotype BLUPs as the phenotypic input, the first three principal components as fixed covariates to correct for population structure, and a kinship matrix as a random effect to correct for relatedness among the 290 accessions, testing 7,110,996 SNPs (Fig. 2). The genomic inflation factor was close to unity for all four traits (λ~GC~ = 0.975–1.013; Fig. 3), indicating adequate control of population structure and of genome-wide false-positive inflation.

---

## Change 8 of 25 · Results Ch. 2 ¶1 (2/2)

*Why:* Renumbering: old Fig. 3 is now Fig. 2.

**With the changes marked:**

To identify the genomic regions associated with grain nutritional traits, genome-wide association was performed for each trait using the genotype BLUPs as the phenotypic input, the first three principal components as fixed covariates to correct for population structure, and a kinship matrix as a random effect to correct for relatedness among the 290 accessions, testing 7,110,996 SNPs (Fig. 3). The genomic inflation factor was close to unity for all four traits (λ~GC~ = 0.975–1.013; Fig. <del style="color:#c0392b;background:#fdecea">3</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2</ins>), indicating adequate control of population structure and of genome-wide false-positive inflation.

**Reads after approval:**

> To identify the genomic regions associated with grain nutritional traits, genome-wide association was performed for each trait using the genotype BLUPs as the phenotypic input, the first three principal components as fixed covariates to correct for population structure, and a kinship matrix as a random effect to correct for relatedness among the 290 accessions, testing 7,110,996 SNPs (Fig. 3). The genomic inflation factor was close to unity for all four traits (λ~GC~ = 0.975–1.013; Fig. 2), indicating adequate control of population structure and of genome-wide false-positive inflation.

---

## Change 9 of 25 · Results Ch. 2 ¶2

*Why:* Renumbering: old Fig. 3 is now Fig. 2.

**With the changes marked:**

Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. <del style="color:#c0392b;background:#fdecea">3</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2</ins>; Online Resource 1<sup>(TODO 301)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource N<sup>(TODO 322)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many. The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which resolved into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(TODO 309)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci. Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 1<sup>(TODO 301)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource N<sup>(TODO 322)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many. The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which resolved into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(TODO 309)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci. Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 10 of 25 · Results Ch. 2 ¶2 (your comment 311)

*Why:* Your comment names the Manhattan figure; renumbered only.

**With the changes marked:**

Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 3; Online Resource 1<sup>(TODO 301)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource N<sup>(TODO 322)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many. The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which resolved into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(TODO 309)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci. Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 3; Online Resource 1<sup>(TODO 301)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource N<sup>(TODO 322)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many. The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which resolved into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(TODO 309)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci. Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 11 of 25 · Results Ch. 2, Manhattan caption

*Why:* Renumbering: old Fig. 3 is now Fig. 2.

**With the changes marked:**

**Fig. <del style="color:#c0392b;background:#fdecea">3</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2</ins>** Genome-wide association mapping of grain nutritional traits in 290 *H. spontaneum* accessions. Manhattan plots (left) and quantile–quantile plots (right) for **a** β-glucan, **b** fiber, **c** protein and **d** starch (7,110,996 SNPs). In the Manhattan plots, the SNPs of each of the 36 loci are colored by locus over the alternating gray chromosomes, with adjacent loci in different colors; the horizontal line marks the Bonferroni threshold (−log~10~(*P*) = 6.0454). λ~GC~, genomic inflation factor

**Reads after approval:**

> **Fig. 2** Genome-wide association mapping of grain nutritional traits in 290 *H. spontaneum* accessions. Manhattan plots (left) and quantile–quantile plots (right) for **a** β-glucan, **b** fiber, **c** protein and **d** starch (7,110,996 SNPs). In the Manhattan plots, the SNPs of each of the 36 loci are colored by locus over the alternating gray chromosomes, with adjacent loci in different colors; the horizontal line marks the Bonferroni threshold (−log~10~(*P*) = 6.0454). λ~GC~, genomic inflation factor

---

## Change 12 of 25 · Results Ch. 3 ¶3 (GPAT6)

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

The first of the two fiber genes, the glycerol-3-phosphate 2-*O*-acyltransferase *GPAT6* (*HORVU.MOREX.r3.3HG0301300*), resolved into three haplotype groups that differed in grain fiber (*q* = 5.3 × 10⁻⁶, η² = 0.208; Online Resource N<sup>(TODO 323)</sup>; Fig. <del style="color:#c0392b;background:#fdecea">4a</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3a</ins>). It belongs to the land-plant GPAT4/6/8 clade, whose members are bifunctional acyltransferase–phosphatases producing *sn*-2 monoacylglycerol, the committed precursor of cutin; cutin and suberin of the pericarp and testa are recovered in the insoluble fiber fraction (Yang et al. 2010; Petit et al. 2016)<sup>(TODO 304)</sup>.

**Reads after approval:**

> The first of the two fiber genes, the glycerol-3-phosphate 2-*O*-acyltransferase *GPAT6* (*HORVU.MOREX.r3.3HG0301300*), resolved into three haplotype groups that differed in grain fiber (*q* = 5.3 × 10⁻⁶, η² = 0.208; Online Resource N<sup>(TODO 323)</sup>; Fig. 3a). It belongs to the land-plant GPAT4/6/8 clade, whose members are bifunctional acyltransferase–phosphatases producing *sn*-2 monoacylglycerol, the committed precursor of cutin; cutin and suberin of the pericarp and testa are recovered in the insoluble fiber fraction (Yang et al. 2010; Petit et al. 2016)<sup>(TODO 304)</sup>.

---

## Change 13 of 25 · Results Ch. 3 ¶4 (GH17)

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

The second, the glucan endo-1,3-β-glucosidase *GH17* (*HORVU.MOREX.r3.5HG0487060*), resolved into five haplotype groups that differed in grain fiber (*q* = 7.5 × 10⁻⁵, η² = 0.222; Online Resource N<sup>(TODO 324)</sup>; Fig. <del style="color:#c0392b;background:#fdecea">4b</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3b</ins>). It carries the glycoside hydrolase family 17 catalytic domain together with an X8 carbohydrate-binding module. Barley’s own (1,3;1,4)-β-D-glucan endohydrolases EI and EII belong to this family, and mixed-linkage β-glucan is the dominant soluble fiber of the grain (Hrmova and Fincher 2024)<sup>(TODO 305)</sup>. This gene, on 5H, is an uncharacterized member of the family rather than either endohydrolase, which map to 1H and 7H.

**Reads after approval:**

> The second, the glucan endo-1,3-β-glucosidase *GH17* (*HORVU.MOREX.r3.5HG0487060*), resolved into five haplotype groups that differed in grain fiber (*q* = 7.5 × 10⁻⁵, η² = 0.222; Online Resource N<sup>(TODO 324)</sup>; Fig. 3b). It carries the glycoside hydrolase family 17 catalytic domain together with an X8 carbohydrate-binding module. Barley’s own (1,3;1,4)-β-D-glucan endohydrolases EI and EII belong to this family, and mixed-linkage β-glucan is the dominant soluble fiber of the grain (Hrmova and Fincher 2024)<sup>(TODO 305)</sup>. This gene, on 5H, is an uncharacterized member of the family rather than either endohydrolase, which map to 1H and 7H.

---

## Change 14 of 25 · Results Ch. 3 ¶5 (PHT4;3)

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

For starch, the plastid-localized anion/phosphate transporter *PHT4;3* (*HORVU.MOREX.r3.3HG0301710*), an ortholog of rice PHT4;3, resolved into four haplotype groups that differed in grain starch (*q* = 3.8 × 10⁻⁴, η² = 0.078; Online Resource N<sup>(TODO 325)</sup>; Fig. <del style="color:#c0392b;background:#fdecea">4c</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3c</ins>). Phosphate is the allosteric inhibitor of ADP-glucose pyrophosphorylase, the committed step of starch synthesis, so the rate at which phosphate crosses the plastid envelope constrains starch accumulation. Consistent with this role, loss of the related plastidic transporter PHT4;2 alters starch accumulation in *Arabidopsis* (Guo et al. 2008; Irigoyen et al. 2011)<sup>(TODO 307)</sup>.

**Reads after approval:**

> For starch, the plastid-localized anion/phosphate transporter *PHT4;3* (*HORVU.MOREX.r3.3HG0301710*), an ortholog of rice PHT4;3, resolved into four haplotype groups that differed in grain starch (*q* = 3.8 × 10⁻⁴, η² = 0.078; Online Resource N<sup>(TODO 325)</sup>; Fig. 3c). Phosphate is the allosteric inhibitor of ADP-glucose pyrophosphorylase, the committed step of starch synthesis, so the rate at which phosphate crosses the plastid envelope constrains starch accumulation. Consistent with this role, loss of the related plastidic transporter PHT4;2 alters starch accumulation in *Arabidopsis* (Guo et al. 2008; Irigoyen et al. 2011)<sup>(TODO 307)</sup>.

---

## Change 15 of 25 · Results Ch. 3 ¶6 (elite)

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

To ask which of these haplotypes are present in cultivated barley, five European spring malting cultivars released between 2012 and 2018 were genotyped across the same three windows from an independent MorexV3 call set. Their genotypes were drawn beneath the wild haplotype groups (Fig. <del style="color:#c0392b;background:#fdecea">4d–f</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3d–f</ins>; Online Resource 3<sup>(comment 319)</sup><sup>(TODO 306)</sup>). At the acyltransferase (*GPAT6*) the five cultivars matched the two low-fiber haplotype groups closely, at 95–98% and 88–91% of called sites, and the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. Like the two low-fiber groups, the cultivars carried the reference allele at only 12–14% of these sites, against 65% in the high-fiber group. At the glucanase (*GH17*) and the phosphate transporter (*PHT4;3*), by contrast, all five cultivars carried the reference allele at every called site.

**Reads after approval:**

> To ask which of these haplotypes are present in cultivated barley, five European spring malting cultivars released between 2012 and 2018 were genotyped across the same three windows from an independent MorexV3 call set. Their genotypes were drawn beneath the wild haplotype groups (Fig. 3d–f; Online Resource 3<sup>(comment 319)</sup><sup>(TODO 306)</sup>). At the acyltransferase (*GPAT6*) the five cultivars matched the two low-fiber haplotype groups closely, at 95–98% and 88–91% of called sites, and the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. Like the two low-fiber groups, the cultivars carried the reference allele at only 12–14% of these sites, against 65% in the high-fiber group. At the glucanase (*GH17*) and the phosphate transporter (*PHT4;3*), by contrast, all five cultivars carried the reference allele at every called site.

---

## Change 16 of 25 · Results Ch. 3, haplotype caption

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

**Fig. <del style="color:#c0392b;background:#fdecea">4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3</ins>** Haplotype structure of the three candidate genes carried forward, and the haplotypes of five elite malting cultivars. **a**–**c** Grain-trait BLUPs by wild haplotype group for **a** *HORVU.MOREX.r3.3HG0301300* (GPAT6, fiber), **b** *HORVU.MOREX.r3.5HG0487060* (GH17, fiber) and **c** *HORVU.MOREX.r3.3HG0301710* (PHT4;3, starch); violins with boxplots, group size below each. Above each panel, the Benjamini–Hochberg *q* and η² of the gene-level Kruskal–Wallis test; brackets, Wilcoxon comparisons of each group against the largest, Holm-adjusted (\*\*\*\**P* ≤ 10⁻⁴, \*\*\**P* ≤ 10⁻³, \*\**P* ≤ 0.01, \**P* ≤ 0.05; ns, not significant). **d**–**f** Genotypes at the same genes: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

**Reads after approval:**

> **Fig. 3** Haplotype structure of the three candidate genes carried forward, and the haplotypes of five elite malting cultivars. **a**–**c** Grain-trait BLUPs by wild haplotype group for **a** *HORVU.MOREX.r3.3HG0301300* (GPAT6, fiber), **b** *HORVU.MOREX.r3.5HG0487060* (GH17, fiber) and **c** *HORVU.MOREX.r3.3HG0301710* (PHT4;3, starch); violins with boxplots, group size below each. Above each panel, the Benjamini–Hochberg *q* and η² of the gene-level Kruskal–Wallis test; brackets, Wilcoxon comparisons of each group against the largest, Holm-adjusted (\*\*\*\**P* ≤ 10⁻⁴, \*\*\**P* ≤ 10⁻³, \*\**P* ≤ 0.01, \**P* ≤ 0.05; ns, not significant). **d**–**f** Genotypes at the same genes: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

---

## Change 17 of 25 · Results Ch. 4 ¶1

*Why:* Renumbering: old Fig. 3 is now Fig. 2.

**With the changes marked:**

To search for a candidate gene at the tightest of the three fiber–starch co-localizations, the region around the two 7H loci was examined (Fig. <del style="color:#c0392b;background:#fdecea">3b, d</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2b, d</ins>). Both loci were too narrow to contain a gene, spanning 185 bp (fiber) and 1.56 kb (starch). A single gene, *HORVU.MOREX.r3.7HG0729030*, lay immediately beside them, within 732 bp of both lead SNPs and 554 bp beyond the edge of the starch locus, and the nearest other annotated gene was 112 kb away. The gene encodes a GDSL esterase/lipase predicted to be secreted, identified at the level of the family rather than as a characterized ortholog. Members of this family act on ester bonds in the plant cell wall and cuticle. In rice, the GDSL esterase DARX1 removes acetyl groups from arabinoxylan, a major fiber polysaccharide of the cereal grain, and in tomato, GDSL1/CD1 polymerizes cutin (Zhang et al. 2019; Girard et al. 2012)<sup>(TODO 308)</sup>.

**Reads after approval:**

> To search for a candidate gene at the tightest of the three fiber–starch co-localizations, the region around the two 7H loci was examined (Fig. 2b, d). Both loci were too narrow to contain a gene, spanning 185 bp (fiber) and 1.56 kb (starch). A single gene, *HORVU.MOREX.r3.7HG0729030*, lay immediately beside them, within 732 bp of both lead SNPs and 554 bp beyond the edge of the starch locus, and the nearest other annotated gene was 112 kb away. The gene encodes a GDSL esterase/lipase predicted to be secreted, identified at the level of the family rather than as a characterized ortholog. Members of this family act on ester bonds in the plant cell wall and cuticle. In rice, the GDSL esterase DARX1 removes acetyl groups from arabinoxylan, a major fiber polysaccharide of the cereal grain, and in tomato, GDSL1/CD1 polymerizes cutin (Zhang et al. 2019; Girard et al. 2012)<sup>(TODO 308)</sup>.

---

## Change 18 of 25 · Results Ch. 4 ¶3 (1/4)

*Why:* Renumbering: old Fig. 5 is now Fig. 4.

**With the changes marked:**

In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. <del style="color:#c0392b;background:#fdecea">5</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4</ins>). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 5a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 5b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 5c).

**Reads after approval:**

> In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 4). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 5a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 5b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 5c).

---

## Change 19 of 25 · Results Ch. 4 ¶3 (2/4)

*Why:* Renumbering: old Fig. 5 is now Fig. 4.

**With the changes marked:**

In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 5). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. <del style="color:#c0392b;background:#fdecea">5a</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4a</ins>) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 5b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 5c).

**Reads after approval:**

> In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 5). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 4a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 5b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 5c).

---

## Change 20 of 25 · Results Ch. 4 ¶3 (3/4)

*Why:* Renumbering: old Fig. 5 is now Fig. 4.

**With the changes marked:**

In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 5). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 5a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. <del style="color:#c0392b;background:#fdecea">5b</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4b</ins>) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 5c).

**Reads after approval:**

> In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 5). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 5a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 4b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 5c).

---

## Change 21 of 25 · Results Ch. 4 ¶3 (4/4)

*Why:* Renumbering: old Fig. 5 is now Fig. 4.

**With the changes marked:**

In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 5). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 5a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 5b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. <del style="color:#c0392b;background:#fdecea">5c</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4c</ins>).

**Reads after approval:**

> In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 5). The haplotypes were defined by a single group of four SNPs, the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 5a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 5b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 2b). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 4c).

---

## Change 22 of 25 · Results Ch. 4, GDSL caption

*Why:* Renumbering: old Fig. 5 is now Fig. 4.

**With the changes marked:**

**Fig. <del style="color:#c0392b;background:#fdecea">5</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4</ins>** Haplotype structure of the GDSL esterase/lipase *HORVU.MOREX.r3.7HG0729030*, and the haplotypes of five elite malting cultivars. **a**, **b** Grain-trait BLUPs by wild haplotype group for **a** fiber and **b** starch; violins with boxplots, group size below each. Above each panel, the Kruskal–Wallis *P* and η²; brackets, Wilcoxon comparison of the two groups (\*\**P* ≤ 0.01). **c** Genotypes at the same gene: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

**Reads after approval:**

> **Fig. 4** Haplotype structure of the GDSL esterase/lipase *HORVU.MOREX.r3.7HG0729030*, and the haplotypes of five elite malting cultivars. **a**, **b** Grain-trait BLUPs by wild haplotype group for **a** fiber and **b** starch; violins with boxplots, group size below each. Above each panel, the Kruskal–Wallis *P* and η²; brackets, Wilcoxon comparison of the two groups (\*\**P* ≤ 0.01). **c** Genotypes at the same gene: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

---

## Change 23 of 25 · Discussion §3 ¶3 (β-glucan)

*Why:* Renumbering: old Fig. 3 is now Fig. 2.

**With the changes marked:**

β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP<sup>(comment 414)</sup><sup>(TODO 410)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. <del style="color:#c0392b;background:#fdecea">3a</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2a</ins>). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Burton et al. 2008)<sup>(TODO 409)</sup>, is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors attributed this conservation to the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014)<sup>(TODO 408)</sup>, indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

**Reads after approval:**

> β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP<sup>(comment 414)</sup><sup>(TODO 410)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Burton et al. 2008)<sup>(TODO 409)</sup>, is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors attributed this conservation to the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014)<sup>(TODO 408)</sup>, indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

---

## Change 24 of 25 · Discussion §5 (1/2)

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

The five elite cultivars shared a single haplotype at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. <del style="color:#c0392b;background:#fdecea">4f</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3f</ins>), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 4e). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resource X<sup>(TODO 420)</sup>). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

**Reads after approval:**

> The five elite cultivars shared a single haplotype at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 3f), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 4e). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resource X<sup>(TODO 420)</sup>). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

---

## Change 25 of 25 · Discussion §5 (2/2)

*Why:* Renumbering: old Fig. 4 is now Fig. 3.

**With the changes marked:**

The five elite cultivars shared a single haplotype at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 4f), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. <del style="color:#c0392b;background:#fdecea">4e</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3e</ins>). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resource X<sup>(TODO 420)</sup>). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

**Reads after approval:**

> The five elite cultivars shared a single haplotype at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 4f), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 3e). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resource X<sup>(TODO 420)</sup>). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

---

**Your decision:** approve · discard · or tell me what to change.
