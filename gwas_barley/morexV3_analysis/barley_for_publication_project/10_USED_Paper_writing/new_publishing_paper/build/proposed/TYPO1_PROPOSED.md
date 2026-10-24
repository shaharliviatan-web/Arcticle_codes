# Review item TYPO1: proposed changes

Proposal only. `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Applied 2026-10-17 as decided by the user:** changes 1 and 2 only (the open-items line; protein–starch with the article). Not applied: changes 3–5 (inside Word comments; user: do not touch comments or markup) and change 6 (numerals in a series: not a TAG requirement; the TAG guidelines set no rule on number style or US/British spelling, only correct grammar, spelling and punctuation). The reading-copy generator fix (Q4) was approved.
>
> **Item TYPO1 — typos, spelling and spacing** (user request 2026-10-17: a proofreading pass over the whole manuscript). Checks run: a spell check of all visible text and of every Word comment (pyspellchecker; every flagged word reviewed by hand: the rest are names, gene and tool names, Latin and technical terms, all correct); double spaces, spaces before punctuation, repeated words, unbalanced brackets, hyphens in number ranges, British spellings; and a full read for wrong-but-real words and grammar. **No misspelled word was found in the manuscript text.** What was found is below (changes 1–6).
>
> **Changes 3–5 are inside Word comments.** The preview hides comment text, so for these three the exact before → after is given in the *Why* line; the paragraph shown has no visible change.
>
> **Double spaces.** In the main file, the only real double spaces are the two in your reply to Sariel (change 5). The other double spaces you see in the reading copy are where a hidden note or Word comment was removed by the reading-copy generator; Word and the markdown preview never show them. **Proposed together with this item (Q4):** one line added to `build/make_reading_copy.py` so that it collapses those leftover runs of spaces into one; the main file is not affected.
>
> **Noticed, not changed (tell me if you want an item for them):** (a) three abbreviations are defined twice: SNP and GWAS in the Introduction and again in the Materials and methods, MAF in the Materials and methods and again in the Results (TAG: define at first mention, then use consistently); (b) Sariel's Results ch. 1 sentences use the serial comma ("clay, precipitation, and organic carbon"), the rest of the text does not.
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** Change 2 is in Sariel's sentence: besides the en dash ("protein–starch", as in "starch–fiber" elsewhere), it adds the missing article ("and **a** significant negative correlation was obtained"). Approve both, or only the dash?
>
> **Q2.** Change 6 is optional style: in a series that includes a number above nine, Springer style uses numerals throughout ("2, 11 and 18 loci" instead of "two, 11 and 18 loci"). Approve or discard?
>
> **Q3.** Your reply to Sariel (comment on "unresolved in this panel") says "…or just from the **embryonate** or ??". Did you mean "embryo"? Change 5 fixes only the two double spaces; say if you want "embryo" too.
>
> **Q4.** Approve the one-line change to `build/make_reading_copy.py` (collapse the double spaces left where hidden notes are removed; reading copy only)?
>

---

## Change 1 of 2 · Open-items list, item 9 (figure accessibility)

*Why:* US spelling, as in the rest of the text ("colored", "gray"); and the manuscript has had three figures since the merge, so the color-only figures are now Figs. 2 and 3.

**With the changes marked:**

1. **Title page (not written yet):** the author names and their affiliations (institution, department, city, country); the corresponding author's e-mail and ORCID (ORCID is mandatory); and the Acknowledgments, which TAG places on the title page, with the funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, affiliation and e-mail also go into the title block of every Online Resource file, where they are placeholders now; the files are then rebuilt.
2. **Author contribution statement (not written yet):** one short line per author, using initials. TAG prints it before the Acknowledgments, and it is also entered in the submission system.
3. **Statements and Declarations (not written yet; they go after the References):** TAG returns a manuscript that lacks any of them. Needed: Funding (agencies and grant numbers); Competing interests; Author contributions; Data availability; Ethics approval (not applicable to this study). The Data availability text is already drafted and names Online Resource 1 for the genotype BLUPs; it still needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub link (with its Zenodo archive) of the analysis code (comment 234).
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary): done.** The 13 files are built, numbered in the order they are first cited, cited in the text, and captioned in the last section ("Supplementary Information"). Left to do: put the author names, affiliation and e-mail into the files (item 1); at submission, upload the 13 files and enter each one's caption in the submission system; and decide whether the extra, uncited supplement of raw-phenotype checks should join the supplement (it would then need a number and a citation in the text).
6. **Other open comments** (each marked by a Word comment).
7. **Cover letter (not written yet):** give the accession number of the raw reads (ENA PRJEB79623); state any restriction on sharing the seeds of the accessions (the same statement belongs at the end of the Materials and methods); and suggest reviewers who are independent of the authors, from different countries, with institutional e-mail addresses.
8. **Final Word file:** 10-pt Times Roman; automatic page numbering; no field functions. The three figures are embedded in the text and are also uploaded as separate TIFF files, Fig1 to Fig3 (600 dpi, 174 mm wide). Before submitting, delete this list and every remaining Word comment.
9. *(Optional)* **Figure accessibility:** TAG recommends patterns or symbols in addition to <del style="color:#c0392b;background:#fdecea">colour</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">color</ins>, but Figs. 2<del style="color:#c0392b;background:#fdecea">–4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> and 3</ins> tell loci and haplotype groups apart by <del style="color:#c0392b;background:#fdecea">colour</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">color</ins> alone.

**Reads after approval:**

> 1. **Title page (not written yet):** the author names and their affiliations (institution, department, city, country); the corresponding author's e-mail and ORCID (ORCID is mandatory); and the Acknowledgments, which TAG places on the title page, with the funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, affiliation and e-mail also go into the title block of every Online Resource file, where they are placeholders now; the files are then rebuilt.
2. **Author contribution statement (not written yet):** one short line per author, using initials. TAG prints it before the Acknowledgments, and it is also entered in the submission system.
3. **Statements and Declarations (not written yet; they go after the References):** TAG returns a manuscript that lacks any of them. Needed: Funding (agencies and grant numbers); Competing interests; Author contributions; Data availability; Ethics approval (not applicable to this study). The Data availability text is already drafted and names Online Resource 1 for the genotype BLUPs; it still needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub link (with its Zenodo archive) of the analysis code (comment 234).
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary): done.** The 13 files are built, numbered in the order they are first cited, cited in the text, and captioned in the last section ("Supplementary Information"). Left to do: put the author names, affiliation and e-mail into the files (item 1); at submission, upload the 13 files and enter each one's caption in the submission system; and decide whether the extra, uncited supplement of raw-phenotype checks should join the supplement (it would then need a number and a citation in the text).
6. **Other open comments** (each marked by a Word comment).
7. **Cover letter (not written yet):** give the accession number of the raw reads (ENA PRJEB79623); state any restriction on sharing the seeds of the accessions (the same statement belongs at the end of the Materials and methods); and suggest reviewers who are independent of the authors, from different countries, with institutional e-mail addresses.
8. **Final Word file:** 10-pt Times Roman; automatic page numbering; no field functions. The three figures are embedded in the text and are also uploaded as separate TIFF files, Fig1 to Fig3 (600 dpi, 174 mm wide). Before submitting, delete this list and every remaining Word comment.
9. *(Optional)* **Figure accessibility:** TAG recommends patterns or symbols in addition to color, but Figs. 2 and 3 tell loci and haplotype groups apart by color alone.

---

## Change 2 of 2 · Results ch. 1 ¶2 (trait correlations)

*Why:* En dash for a pair of traits (as "starch–fiber", "fiber–starch" elsewhere) and the missing article (Q1).

**With the changes marked:**

Among the nutritional traits, starch was strongly negatively correlated with fiber (r = −0.78; Fig. 1c) and β-glucan (r = −0.53), while fiber and β-glucan were positively correlated (r = 0.35; all *P* \< 0.001). Thus, genotypes with higher starch tended to have lower fiber and β-glucan. Protein was poorly correlated with other traits, and <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">a </ins>significant negative correlation was obtained only for <del style="color:#c0392b;background:#fdecea">protein-starch</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">protein–starch</ins> (r = −0.12, *P* \< 0.05). Among the morphological and phenological traits, starch was strongly correlated with grain weight, plant height, and flowering time, while fiber showed the exact opposite pattern. Protein was also weakly negatively correlated with grain weight (r = −0.13, *P* \< 0.05). Overall, β-glucan and protein displayed weaker, mostly negative correlations, linking grain composition directly to whole-plant resource allocation and phenology (Fig. 1c).

**Reads after approval:**

> Among the nutritional traits, starch was strongly negatively correlated with fiber (r = −0.78; Fig. 1c) and β-glucan (r = −0.53), while fiber and β-glucan were positively correlated (r = 0.35; all *P* \< 0.001). Thus, genotypes with higher starch tended to have lower fiber and β-glucan. Protein was poorly correlated with other traits, and a significant negative correlation was obtained only for protein–starch (r = −0.12, *P* \< 0.05). Among the morphological and phenological traits, starch was strongly correlated with grain weight, plant height, and flowering time, while fiber showed the exact opposite pattern. Protein was also weakly negatively correlated with grain weight (r = −0.13, *P* \< 0.05). Overall, β-glucan and protein displayed weaker, mostly negative correlations, linking grain composition directly to whole-plant resource allocation and phenology (Fig. 1c).

---

**Your decision:** approve · discard · or tell me what to change.
