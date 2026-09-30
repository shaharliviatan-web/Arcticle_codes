# Review item MM5: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM5 of the M&M voice review** (item 5, closing paragraphs; approved in principle in the chat). Today the last three paragraphs (graphics, data, LLM statement) have no heading, so in Word they read as part of "Geographic origin of allele and haplotype carriers". The two other parts of item 5 (the 7H section, and shortening Geographic origin) come next as MM6 and MM7.
>
> **New layout of the M&M ending:** … Geographic origin → **## Data and material availability** (the ENA accession number + the material-sharing placeholder, TODO 235) → **## Use of generative artificial intelligence** (the LLM statement, unchanged, TODO 236).
>
> **Moved to the Declarations (not written yet):** the sentence on the Zenodo VCF, the BLUPs Online Resource and the GitHub code. On approval it is saved word for word in `build/DECLARATIONS_HANDOVER.md` for the Declarations session. TODO 234 keeps its full list for Sariel and says where that sentence went.
>
> **Graphics sentence:** moves into Phenotypic data analysis ¶1, where R is first named, with its TODO 233 (the citations to check). "Figures were drawn in R with ggplot2 v4.0.1 …"
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** **Heading names:** OK with "Data and material availability" and "Use of generative artificial intelligence"? I avoided the plain "Data availability", because the Declarations will carry a section with that exact title.
>
> **Q2.** **Graphics sentence placement:** in Phenotypic data analysis ¶1, right after R is named (proposed)? The alternative is to keep it at the end of the M&M, as the first sentence under the new "Data and material availability" heading. It fits there less well, but it would not interrupt the phenotypic paragraph.
>

---

## Change 1 of 4 · M&M · Phenotypic data analysis ¶1

*Why:* The graphics-software sentence joins the first mention of R, as Evgeny names plotting packages inside the method they belong to; TODO 233 moves with it.

**With the changes marked:**

All phenotypic analyses were performed in R v4.1.2 <del style="color:#c0392b;background:#fdecea">(R Core Team 2021)</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">(R Core Team 2021)<sup>(TODO 233)</sup></ins>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. </ins> Block 4 of the 2019–20 season was excluded, because its β-glucan values were poorly correlated with those of the other blocks (Pearson r = 0.27–0.34, against 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The grain-composition data thus comprised 3,164 plants of the 290 accessions. Each trait was mean-centered within season by subtracting the season mean, and for each trait, missing values were excluded and values deviating by more than three standard deviations from the mean of the centered values were removed as outliers.

**Reads after approval:**

> All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup>. Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded, because its β-glucan values were poorly correlated with those of the other blocks (Pearson r = 0.27–0.34, against 0.61–0.64 among blocks 1–3; Online Resource N<sup>(TODO 215)</sup>). The grain-composition data thus comprised 3,164 plants of the 290 accessions. Each trait was mean-centered within season by subtracting the season mean, and for each trait, missing values were excluded and values deviating by more than three standard deviations from the mean of the centered values were removed as outliers.

---

## Change 2 of 4 · M&M · closing paragraph 1 (graphics)

*Why:* The sentence moves to Phenotypic ¶1; in its place, the heading of the closing subsection.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Figures were drawn in R v4.1.2 (R Core Team 2021)<sup>(TODO 233)</sup> with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. </del><span style="color:#1f4e9c;font-weight:bold;font-style:italic"> (→ Phenotypic data analysis ¶1) </span><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">## Data and material availability</ins>

**Reads after approval:**

> ## Data and material availability

---

## Change 3 of 4 · M&M · closing paragraph 2 (data)

*Why:* The M&M keeps the accession number (TAG: at the end of M&M); Zenodo, GitHub and BLUPs go to the Declarations' Data availability, where both of Evgeny's papers put them. TODO 234 keeps its full list for Sariel and says where the sentence went.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">[TODO: Sariel needs to sort this out and approve it. To fill in and confirm: (1) the MorexV3 SNP calls on Zenodo with a DOI (user choice 2026-09-30; /mnt/data/shahar/gwas_barley/data/inputs/morexV3_with_ids.vcf.gz, 6.9 GB, 300 accessions x 7,110,996 SNPs - decide whether to deposit the 300- or the 290-accession file); (2) the GitHub repository under https://github.com/hubner-lab (as in Potapenko et al. 2026 Mol Ecol and Mol Biol Evol) and a Zenodo archive of its release for a DOI; (3) phenotypes as in Potapenko et al.: the genotype BLUPs as an Online Resource (00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables/BLUP_10_traits.csv), no separate deposit; (4) the ENA sample alias of each accession (VCF HS0103 = ENA alias 01_03) in the sampling-site Online Resource (TODO 218)</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026a); the sample alias of each accession is listed in Online Resource N.<del style="color:#c0392b;background:#fdecea"> The MorexV3 SNP calls are available from Zenodo (https://doi.org/XXXX), the genotype BLUPs are given in Online Resource N, and the analysis code is available at https://github.com/hubner-lab/XXXX (archived at https://doi.org/XXXX).</del><span style="color:#1f4e9c;font-weight:bold;font-style:italic"> (→ Declarations, Data availability) </span><sup>(TODO 234)</sup> material-sharing statement<sup>(TODO 235)</sup>

**Reads after approval:**

> The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026a); the sample alias of each accession is listed in Online Resource N.<sup>(TODO 234)</sup> material-sharing statement<sup>(TODO 235)</sup>

---

## Change 4 of 4 · M&M · closing paragraph 3 (LLM statement)

*Why:* Own heading, so it no longer reads as part of 'Geographic origin…' or of the data paragraph. Text unchanged.

**With the changes marked:**

<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">## Use of generative artificial intelligence

</ins>Generative artificial intelligence tools, Claude Opus and Claude Sonnet (Anthropic; used through Claude Code) and ChatGPT with GPT-5 and GPT-5.5 (OpenAI), were used between February and October 2026 to (i) discuss analysis strategies and the interpretation of results; (ii) write, debug and document the R, Python and shell scripts of the analysis; (iii) check the consistency of the manuscript with the code and its outputs; and (iv) draft and revise the manuscript text. All analytical decisions and interpretations were made by the authors, all code and its outputs were checked by the authors, and all text was reviewed and edited by the authors, who take full responsibility for the content of this publication.<sup>(TODO 236)</sup>

**Reads after approval:**

> ## Use of generative artificial intelligence

Generative artificial intelligence tools, Claude Opus and Claude Sonnet (Anthropic; used through Claude Code) and ChatGPT with GPT-5 and GPT-5.5 (OpenAI), were used between February and October 2026 to (i) discuss analysis strategies and the interpretation of results; (ii) write, debug and document the R, Python and shell scripts of the analysis; (iii) check the consistency of the manuscript with the code and its outputs; and (iv) draft and revise the manuscript text. All analytical decisions and interpretations were made by the authors, all code and its outputs were checked by the authors, and all text was reviewed and edited by the authors, who take full responsibility for the content of this publication.<sup>(TODO 236)</sup>

---

**Your decision:** approve · discard · or tell me what to change.
