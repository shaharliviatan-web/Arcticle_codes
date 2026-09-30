# Review item S5: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item S5 — the new Fig. 1 and the renamed figure files** (Q1 triangle → lower triangle; Q2 A; Q3 approved). Applied together with S3 and S4.
>
> Fig. 1 caption: a and b word for word; c reworded for the lower triangle; d and the significance line from the old Fig. 2 caption, stating that one test (raw two-sided P, no multiple-testing correction) gives all stars in c and d.
>
> Image links and src notes follow the renamed files (Figs. 2–4 re-rendered pixel-identical). A hidden note at the top of the file records the old → new numbers.
>

## Change 1 of 11 · File top (hidden note)

*Why:* Record the old → new numbering once, for the hidden notes that keep the old numbers.

**With the changes marked:**

<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">

</ins>

**Reads after approval:**

> *(paragraph removed; nothing is left of it in Word)*

---

## Change 2 of 11 · M&M · Phenotypic data analysis ¶1 (hidden src note)

*Why:* Scripts renumbered.

**With the changes marked:**

All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup>. Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The nutritional-trait data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

**Reads after approval:**

> All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup>. Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The nutritional-trait data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

---

## Change 3 of 11 · Results Ch. 1, Fig. 1 caption

*Why:* New Fig. 1 (Q1, Q2 approved): title merges the two old figure titles; a, b unchanged; c reworded for the lower triangle; d and the significance line from the old Fig. 2 caption, stating one test for all stars.

**With the changes marked:**

**Fig. 1** Genetic <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">and ecological </ins>architecture of grain nutritional <del style="color:#c0392b;background:#fdecea">and reference </del>traits. **a** Variance partitioning for four nutritional and four morphological traits into genetic, season + block, genotype × environment, and residual components. **b** Per-genotype reaction norms across the three seasons for the four nutritional traits; gray lines, all accessions; colored lines, four highlighted accessions per trait with values in all three seasons: the most stable, the strongest increase and the strongest decrease from 2019–20 to 2021–22, and the strongest crossover between seasons<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. **c** Pearson correlations of the nutritional-trait BLUPs (rows) with one another (lower triangle) and with four morphological-trait BLUPs. **d** Pearson correlations between site-mean nutritional-trait BLUPs and eight environmental variables; cells outlined in black are significant at *P* \< 0.05. Significance (**c**, **d**), from the two-sided test of each correlation, without correction for multiple testing: \**P* \< 0.05, \*\**P* \< 0.01, \*\*\**P* \< 0.001
</ins>

**Reads after approval:**

> **Fig. 1** Genetic and ecological architecture of grain nutritional traits. **a** Variance partitioning for four nutritional and four morphological traits into genetic, season + block, genotype × environment, and residual components. **b** Per-genotype reaction norms across the three seasons for the four nutritional traits; gray lines, all accessions; colored lines, four highlighted accessions per trait with values in all three seasons: the most stable, the strongest increase and the strongest decrease from 2019–20 to 2021–22, and the strongest crossover between seasons. **c** Pearson correlations of the nutritional-trait BLUPs (rows) with one another (lower triangle) and with four morphological-trait BLUPs. **d** Pearson correlations between site-mean nutritional-trait BLUPs and eight environmental variables; cells outlined in black are significant at *P* \< 0.05. Significance (**c**, **d**), from the two-sided test of each correlation, without correction for multiple testing: \**P* \< 0.05, \*\**P* \< 0.01, \*\*\**P* \< 0.001

---

## Change 4 of 11 · Results Ch. 1, old Fig. 2 image

*Why:* Old Fig. 2 dissolved (a → Online Resource, b + c → Fig. 1c, d → Fig. 1d).

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">![](../../08_USED_creating_figures/Figure_2/Fig2.png){width="6.85in"}</del>

**Reads after approval:**

> *(paragraph removed; nothing is left of it in Word)*

---

## Change 5 of 11 · Results Ch. 1, old Fig. 2 caption

*Why:* Removed with the figure; the former caption is kept in a hidden note.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">**Fig. 2** Ecological architecture of grain nutritional traits. **a** Centered trait BLUPs of the four nutritional traits across 29 sampling sites, ordered and colored by ecological region. **b** Pearson correlations among nutritional-trait BLUPs. **c** Pearson correlations between nutritional- and morphological-trait BLUPs. **d** Pearson correlations between site-mean nutritional-trait BLUPs and eight environmental variables; cells outlined in black are significant at *P* \< 0.05. Significance (**b**–**d**): \**P* \< 0.05, \*\**P* \< 0.01, \*\*\**P* \< 0.001</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"></ins>

**Reads after approval:**

> *(paragraph removed; nothing is left of it in Word)*

---

## Change 6 of 11 · Results, image link of the old Fig. 3

*Why:* Renamed output: Figure_3/Fig3.png is now Figure_2/Fig2.png (pixel-identical).

**With the changes marked:**

![](../../08_USED_creating_figures/Figure_<del style="color:#c0392b;background:#fdecea">3</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2</ins>/Fig<del style="color:#c0392b;background:#fdecea">3</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2</ins>.png){width="6.85in"}

**Reads after approval:**

> ![](../../08_USED_creating_figures/Figure_2/Fig2.png){width="6.85in"}

---

## Change 7 of 11 · Results, image link of the old Fig. 4

*Why:* Renamed output: Figure_4/Fig4.png is now Figure_3/Fig3.png (pixel-identical).

**With the changes marked:**

![](../../08_USED_creating_figures/Figure_<del style="color:#c0392b;background:#fdecea">4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3</ins>/Fig<del style="color:#c0392b;background:#fdecea">4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3</ins>.png){width="6.85in"}

**Reads after approval:**

> ![](../../08_USED_creating_figures/Figure_3/Fig3.png){width="6.85in"}

---

## Change 8 of 11 · Results, image link of the old Fig. 5

*Why:* Renamed output: Figure_5/Fig5.png is now Figure_4/Fig4.png (pixel-identical).

**With the changes marked:**

![](../../08_USED_creating_figures/Figure_<del style="color:#c0392b;background:#fdecea">5</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4</ins>/Fig<del style="color:#c0392b;background:#fdecea">5</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">4</ins>.png){width="6.85in"}

**Reads after approval:**

> ![](../../08_USED_creating_figures/Figure_4/Fig4.png){width="6.85in"}

---

## Change 9 of 11 · Results, hidden src note under Fig. 2

*Why:* Script and file renamed.

**With the changes marked:**

**Fig. 2** Genome-wide association mapping of grain nutritional traits in 290 *H. spontaneum* accessions. Manhattan plots (left) and quantile–quantile plots (right) for **a** β-glucan, **b** fiber, **c** protein and **d** starch (7,110,996 SNPs). In the Manhattan plots, the SNPs of each of the 36 loci are colored by locus over the alternating gray chromosomes, with adjacent loci in different colors; the horizontal line marks the Bonferroni threshold (−log~10~(*P*) = 6.0454). λ~GC~, genomic inflation factor

**Reads after approval:**

> **Fig. 2** Genome-wide association mapping of grain nutritional traits in 290 *H. spontaneum* accessions. Manhattan plots (left) and quantile–quantile plots (right) for **a** β-glucan, **b** fiber, **c** protein and **d** starch (7,110,996 SNPs). In the Manhattan plots, the SNPs of each of the 36 loci are colored by locus over the alternating gray chromosomes, with adjacent loci in different colors; the horizontal line marks the Bonferroni threshold (−log~10~(*P*) = 6.0454). λ~GC~, genomic inflation factor

---

## Change 10 of 11 · Results, hidden src note under Fig. 3

*Why:* Script and file renamed.

**With the changes marked:**

**Fig. 3** Haplotype structure of the three candidate genes carried forward, and the haplotypes of five elite malting cultivars. **a**–**c** Grain-trait BLUPs by wild haplotype group for **a** *HORVU.MOREX.r3.3HG0301300* (GPAT6, fiber), **b** *HORVU.MOREX.r3.5HG0487060* (GH17, fiber) and **c** *HORVU.MOREX.r3.3HG0301710* (PHT4;3, starch); violins with boxplots, group size below each. Above each panel, the Benjamini–Hochberg *q* and η² of the gene-level Kruskal–Wallis test; brackets, Wilcoxon comparisons of each group against the largest, Holm-adjusted (\*\*\*\**P* ≤ 10⁻⁴, \*\*\**P* ≤ 10⁻³, \*\**P* ≤ 0.01, \**P* ≤ 0.05; ns, not significant). **d**–**f** Genotypes at the same genes: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

**Reads after approval:**

> **Fig. 3** Haplotype structure of the three candidate genes carried forward, and the haplotypes of five elite malting cultivars. **a**–**c** Grain-trait BLUPs by wild haplotype group for **a** *HORVU.MOREX.r3.3HG0301300* (GPAT6, fiber), **b** *HORVU.MOREX.r3.5HG0487060* (GH17, fiber) and **c** *HORVU.MOREX.r3.3HG0301710* (PHT4;3, starch); violins with boxplots, group size below each. Above each panel, the Benjamini–Hochberg *q* and η² of the gene-level Kruskal–Wallis test; brackets, Wilcoxon comparisons of each group against the largest, Holm-adjusted (\*\*\*\**P* ≤ 10⁻⁴, \*\*\**P* ≤ 10⁻³, \*\**P* ≤ 0.01, \**P* ≤ 0.05; ns, not significant). **d**–**f** Genotypes at the same genes: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

---

## Change 11 of 11 · Results, hidden src note under Fig. 4

*Why:* Script and file renamed.

**With the changes marked:**

**Fig. 4** Haplotype structure of the GDSL esterase/lipase *HORVU.MOREX.r3.7HG0729030*, and the haplotypes of five elite malting cultivars. **a**, **b** Grain-trait BLUPs by wild haplotype group for **a** fiber and **b** starch; violins with boxplots, group size below each. Above each panel, the Kruskal–Wallis *P* and η²; brackets, Wilcoxon comparison of the two groups (\*\**P* ≤ 0.01). **c** Genotypes at the same gene: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

**Reads after approval:**

> **Fig. 4** Haplotype structure of the GDSL esterase/lipase *HORVU.MOREX.r3.7HG0729030*, and the haplotypes of five elite malting cultivars. **a**, **b** Grain-trait BLUPs by wild haplotype group for **a** fiber and **b** starch; violins with boxplots, group size below each. Above each panel, the Kruskal–Wallis *P* and η²; brackets, Wilcoxon comparison of the two groups (\*\**P* ≤ 0.01). **c** Genotypes at the same gene: one row per haplotype-group consensus, then one row per cultivar; columns are the SNPs shared by the two call sets, in genomic order

---

**Your decision:** approve · discard · or tell me what to change.
