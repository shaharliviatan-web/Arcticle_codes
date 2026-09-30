# Review item MM11: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM11, the minor round** (the last item of the M&M voice review).
>
> **One name for the trait set:** the Results, the Discussion and most of the M&M say "nutritional traits", but three places in the M&M still said "grain-composition". "Grain composition" stays where it names the concept, including the heading "Grain composition by near-infrared spectroscopy".
>
> **"Ten traits"** is spelled out once, so the reader does not have to count 4 + 6.
>
> **One stray double space**, left over from an earlier cut, is removed. Pandoc would collapse it anyway.
>
> **Checked, no change needed:** every abbreviation is defined at its first use (SNP, MAF, CTAB, NIR, BLUP, H², GWAS, LD, MGmin, minHap). The SD wording now matches Table 2, since MM10b. The "we" count is 20, close to Evgeny's 17 in his Mol Biol Evol M&M.
>
> **Left as is (outside the M&M):** the Results call flowering time a "morphological" trait, while the M&M says "morphological and phenological". That is S. Hübner's Results wording; tell me if you want it aligned.
>

## Change 1 of 4 · M&M · Common garden ¶2 (last sentence)

*Why:* Trait-set name as in the rest of the paper.

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026a): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four <del style="color:#c0392b;background:#fdecea">grain-composition</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">nutritional</ins> traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026a): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 2 of 4 · M&M · Common garden ¶3 (post-harvest)

*Why:* Same.

**With the changes marked:**

At the end of each season, five spikes per plant were bagged and oven-dried at 35 °C for 30 days. The dried spikes were threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22<sup>(TODO 211)</sup>, and the grain was cleaned with a South Dakota seed blower (Seedburo Equipment Company) and stored for <del style="color:#c0392b;background:#fdecea">grain-composition</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">nutritional</ins> analysis.

**Reads after approval:**

> At the end of each season, five spikes per plant were bagged and oven-dried at 35 °C for 30 days. The dried spikes were threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22<sup>(TODO 211)</sup>, and the grain was cleaned with a South Dakota seed blower (Seedburo Equipment Company) and stored for nutritional analysis.

---

## Change 3 of 4 · M&M · Phenotypic data analysis ¶1

*Why:* Same.

**With the changes marked:**

All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup>. Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The <del style="color:#c0392b;background:#fdecea">grain-composition</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">nutritional-trait</ins> data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

**Reads after approval:**

> All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup>. Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The nutritional-trait data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

---

## Change 4 of 4 · M&M · Phenotypic data analysis ¶2 (BLUP model)

*Why:* 'Ten traits' spelled out once; a stray double space removed.

**With the changes marked:**

To estimate the genetic value of each accession, we fitted a linear mixed model for each of the ten traits <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">(the four nutritional and six morphological and phenological traits) </ins><del style="color:#c0392b;background:#fdecea"> </del>in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as random effects. The genotype best linear unbiased predictions (BLUPs) of this model were used as the phenotypes in all downstream analyses.

**Reads after approval:**

> To estimate the genetic value of each accession, we fitted a linear mixed model for each of the ten traits (the four nutritional and six morphological and phenological traits) in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as random effects. The genotype best linear unbiased predictions (BLUPs) of this model were used as the phenotypes in all downstream analyses.

---

**Your decision:** approve · discard · or tell me what to change.
