# Review item MM9: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM9 of the M&M voice review** (item 6, voice and clarity, part 2 of 3): Phenotypic data analysis ¶1 and ¶4, and Loci and candidate genes ¶1–2. No facts change.
>
> **Phenotypic ¶1 (B15):** the outlier sentence ("for each trait" twice; "the mean of the centered values", which is about 0) is split, and it now says what matters: the mean and SD are pooled over the three seasons, as the code does. "Against" becomes "compared with".
>
> **Phenotypic ¶4:** opens with its aim and uses "we"; "with two-sided tests" is no longer tacked on at the end.
>
> **Loci ¶1 (A8):** "limited to the contiguous run of member SNPs…" now describes the walk: moving outward from the lead, the locus ends at the first gap of more than 50 kb between consecutive member SNPs. "Lead a locus" becomes "serve as lead SNPs". The repeat rule is stated positively (assigned SNPs may join but not lead). The paragraph opens with its aim and uses "we".
>
> **Loci ¶2:** opens with its aim and uses "we"; "annotation … of Ensembl Plants" becomes "from Ensembl Plants".
>
> Length: 283 → 291 words (a clarity pass).
>

## Change 1 of 4 · M&M · Phenotypic data analysis ¶1

*Why:* Outlier sentence split and made exact; 'against' → 'compared with'.

**With the changes marked:**

All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup>. Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded<del style="color:#c0392b;background:#fdecea">,</del> because its β-glucan values <del style="color:#c0392b;background:#fdecea">were poorly correlated</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">correlated poorly</ins> with those of the other blocks (Pearson r = 0.27–0.34, <del style="color:#c0392b;background:#fdecea">against</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">compared with</ins> 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The grain-composition data thus comprised 3,164 plants of the 290 accessions. Each trait was <del style="color:#c0392b;background:#fdecea">mean-centered</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">centered</ins> within season by subtracting the season mean<del style="color:#c0392b;background:#fdecea">, and for each trait, missing</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. Missing</ins> values were <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">then </ins>excluded<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">,</ins> and values <del style="color:#c0392b;background:#fdecea">deviating by </del>more than three standard deviations from the mean of the centered values<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, pooled over the three seasons,</ins> were removed as outliers.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup>. Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The grain-composition data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

---

## Change 2 of 4 · M&M · Phenotypic data analysis ¶4 (trait correlations)

*Why:* Aim first + 'we'.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Pairwise Pearson correlations</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To examine the relationships among traits, we calculated two-sided Pearson correlations</ins> between accession BLUPs (n = 290) <del style="color:#c0392b;background:#fdecea">were calculated </del>among the four nutritional traits (Fig. 2b) and between the nutritional traits and flowering time, plant height, grain weight and grain number (Fig. 2c)<del style="color:#c0392b;background:#fdecea">, with two-sided tests</del>.

**Reads after approval:**

> To examine the relationships among traits, we calculated two-sided Pearson correlations between accession BLUPs (n = 290) among the four nutritional traits (Fig. 2b) and between the nutritional traits and flowering time, plant height, grain weight and grain number (Fig. 2c).

---

## Change 3 of 4 · M&M · Loci and candidate genes ¶1

*Why:* Clumping and gap rule in plain words; aim first + 'we'.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">For each trait, the significant SNPs were grouped into loci</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To define loci, we grouped the significant SNPs of each trait</ins> by LD clumping in PLINK (Purcell et al. 2007)<del style="color:#c0392b;background:#fdecea">, using the genotypes of the 290 accessions</del>. Only <del style="color:#c0392b;background:#fdecea">SNPs above the genome-wide threshold could lead a locus</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">significant SNPs could serve as lead SNPs</ins>, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then <del style="color:#c0392b;background:#fdecea">limited to the contiguous run of member SNPs around its lead SNP, ending at the first gap longer than 50 kb on either side</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs</ins>. Clumping was repeated<del style="color:#c0392b;background:#fdecea">, with SNPs already assigned to a locus excluded as leads,</del> until every significant SNP lay within a locus<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, with SNPs already assigned to a locus allowed to join but not to lead</ins>, giving 36 loci (Online Resource N)<sup>(TODO 220)</sup>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> To define loci, we grouped the significant SNPs of each trait by LD clumping in PLINK (Purcell et al. 2007). Only significant SNPs could serve as lead SNPs, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs. Clumping was repeated until every significant SNP lay within a locus, with SNPs already assigned to a locus allowed to join but not to lead, giving 36 loci (Online Resource N)<sup>(TODO 220)</sup>.

---

## Change 4 of 4 · M&M · Loci and candidate genes ¶2

*Why:* Aim first + 'we'; 'of Ensembl' → 'from Ensembl'.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">The interval of each locus, from its first to its last member SNP, was intersected</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To identify candidate genes, we intersected the interval of each locus, from its first to its last member SNP,</ins> with the MorexV3 high-confidence gene annotation (Mascher et al. 2021) <del style="color:#c0392b;background:#fdecea">of</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">from</ins> Ensembl Plants release 62 (Yates et al. 2025)<sup>(TODO 221)</sup> <del style="color:#c0392b;background:#fdecea">with</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">using</ins> bedtools intersect v2.30.0 (Quinlan and Hall 2010). Every gene overlapping a locus interval was taken as a candidate gene (Online Resource N)<sup>(TODO 222)</sup>.

**Reads after approval:**

> To identify candidate genes, we intersected the interval of each locus, from its first to its last member SNP, with the MorexV3 high-confidence gene annotation (Mascher et al. 2021) from Ensembl Plants release 62 (Yates et al. 2025)<sup>(TODO 221)</sup> using bedtools intersect v2.30.0 (Quinlan and Hall 2010). Every gene overlapping a locus interval was taken as a candidate gene (Online Resource N)<sup>(TODO 222)</sup>.

---

**Your decision:** approve · discard · or tell me what to change.
