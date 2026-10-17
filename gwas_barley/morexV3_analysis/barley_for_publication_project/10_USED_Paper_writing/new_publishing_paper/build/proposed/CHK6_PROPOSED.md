# Review item CHK6: proposed changes

Proposal only. `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **The list as it will read** (exactly; item 4 unchanged). In the changes further down, red struck-through text is the old wording that is deleted; only the green text stays.
>
> 
>
> > **⚠ BEFORE SUBMISSION: open items (TAG requirements; updated 2026-10-01). Delete this list before submitting.**
>
> >
>
> > 1. **Title page (not written yet):** the author names and their affiliations (institution, department, city, country); the corresponding author's e-mail and ORCID (ORCID is mandatory); and the Acknowledgments, which TAG places on the title page, with the funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, affiliation and e-mail also go into the title block of every Online Resource file, where they are placeholders now; the files are then rebuilt.
>
> > 2. **Author contribution statement (not written yet):** one short line per author, using initials. TAG prints it before the Acknowledgments, and it is also entered in the submission system.
>
> > 3. **Statements and Declarations (not written yet; they go after the References):** TAG returns a manuscript that lacks any of them. Needed: Funding (agencies and grant numbers); Competing interests; Author contributions; Data availability; Ethics approval (not applicable to this study). The Data availability text is already drafted and names Online Resource 1 for the genotype BLUPs; it still needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub link (with its Zenodo archive) of the analysis code (comment 234).
>
> > 4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
>
> > 5. **Online Resources (supplementary): done.** The 13 files are built, numbered in the order they are first cited, cited in the text, and captioned in the last section ("Supplementary Information"). Left to do: put the author names, affiliation and e-mail into the files (item 1); at submission, upload the 13 files and enter each one's caption in the submission system; and decide whether the extra, uncited supplement of raw-phenotype checks should join the supplement (it would then need a number and a citation in the text).
>
> > 6. **Other open comments** (each marked by a Word comment).
>
> > 7. **Cover letter (not written yet):** give the accession number of the raw reads (ENA PRJEB79623); state any restriction on sharing the seeds of the accessions (the same statement belongs at the end of the Materials and methods); and suggest reviewers who are independent of the authors, from different countries, with institutional e-mail addresses.
>
> > 8. **Final Word file:** 10-pt Times Roman; automatic page numbering; no field functions. The four figures are embedded in the text and are also uploaded as separate TIFF files, Fig1 to Fig4 (600 dpi, 174 mm wide). Before submitting, delete this list and every remaining Word comment.
>
> > 9. *(Optional)* **Figure accessibility:** TAG recommends patterns or symbols in addition to colour, but Figs. 2–4 tell loci and haplotype groups apart by colour alone.
>
> 
>
> **Item CHK6 (revised twice: no names of who should answer) — the open-items list, the resolved comments, and the open comments rewritten** (user, 2026-10-01). Only the list and the Word comments change; the manuscript text is unchanged except where a deleted comment leaves its anchor as plain text.
>
> **Checklist:** no paths in any item (the map of locations stays in the project rules; I add this rule there on apply). Item 6 is now only its heading. Item 4 is unchanged.
>
> **Deleted (resolved):** your comments 311, 312, 313, 314, 315, 316, 317, 318, 319, 414, 415, 416, and TODO 203 (the reasons are in each change below).
>
> **Left unchanged:** the reference comments of item 4 (102, 104, 216, 257, 428, 501, 502; already detailed by the References session). **Also deleted (change 5, user):** the 37 supmat comments, all resolved.
>
> **The new text of each rewritten comment** (the preview below hides comment texts, so here they are in full):
>
> **201.** TODO - number of rounds of single-seed descent. The text says three rounds, as in the paper cited for the plant material, Potapenko et al. 2026b (Mol Ecol): "This procedure (single-seed-descent and selfing) was repeated three times". The genome-size paper on the same collection, Potapenko et al. 2026a (Mol Biol Evol), says two rounds, and so do our earlier Materials and methods drafts (May and June 2026). N. Pintel's thesis does not say. Please confirm the number of rounds.
>
> **202.** TODO - how the reads were cleaned before the MorexV3 re-call. The text says the raw reads were quality-checked and trimmed as described by Potapenko et al. 2026b (Mol Ecol), that is, with FastQC v0.11.5 and fastp v0.20.0 as in that paper. This is right only if the MorexV3 re-call started from the reads already cleaned for that paper. If the reads were cleaned again for this study, please give the tools, versions and settings, and they will replace the citation.
>
> **204.** TODO - the reason the ten accessions of site 04 (Hachola) were removed. The text now says "because of an error in the sampling site", which needs to be confirmed or replaced. Earlier write-ups gave three different reasons: "missing in the phenotyping experiments" (mini paper), "poor sample quality" (old Materials and methods drafts) and "no phenotype data" (notes in our analysis scripts). None of them fits the data well: the phenotype table does hold site-04 plants in all three seasons (130 plants, 100 with grain-composition values), but in 2019-20 and 2020-21 their IDs (for example HS0403, HS0414, HS0415) do not match the ten sequenced accessions (HS0401-HS0411, without HS0403); only 2021-22 matches. Please give the reason to write here.
>
> **205.** TODO - the missing-data filter of the SNPs. The text says that SNPs with more than 30% missing data were removed, but the final SNP set contains SNPs with up to 31.7% missing data (40 SNPs above 30% over the 300 accessions; 155 over the 290 analysed). The likely reason is the order of the filters: the 30% filter was applied before genotypes supported by fewer than three reads were set to missing, which added missing calls afterwards. Please confirm the order and the threshold, so the sentence can state them correctly. The other filter, minor-allele frequency 0.05, holds exactly.
>
> **206.** TODO - the pre-sowing treatment of the seeds. The text gives no details and cites Potapenko et al. 2026b (Mol Ecol), whose protocol is: seeds dried at 30 °C for 14 days, dipped in 4% sodium hypochlorite, sown in mid-October on germination paper, kept 10 days in the dark at 4 °C, then moved to 5-L pots with a commercial soil mix. The sources disagree: Potapenko et al. 2026a (Mol Biol Evol) gives a 3% chlorine solution, sowing trays, 10 days in the dark at 4 °C and 14 days in a growth chamber at 20 °C before transplanting; N. Pintel's thesis and our old drafts give 3% sodium hypochlorite for 5 minutes, a 5-minute rinse in distilled water and germination paper at 16 °C for 18 days. Please confirm which protocol was used in this experiment; if it differs from the cited paper, the details will be written here.
>
> **207.** TODO - the definition of flowering time. The text follows the cited paper, Potapenko et al. 2026b (Mol Ecol): days from sowing to heading. Other sources define it differently: Potapenko et al. 2026a (Mol Biol Evol) counts days from transplanting to the emergence of the awns (Zadoks stage 49), and N. Pintel's thesis and the old draft count days from transplanting to the emergence of the first spike from the flag-leaf sheath. Our data fit the cited definition: in 2020-21 they range from 73 to 147 days, as reported in Mol Ecol. Please confirm.
>
> **208.** TODO - the definition of plant height. The text follows the cited paper, Potapenko et al. 2026b (Mol Ecol): the length of the tallest tiller, which is the measurement used in our analyses. N. Pintel's thesis defines plant height instead as the height from the ground to the base of the highest flag leaf, which Potapenko et al. 2026a (Mol Biol Evol) records as a separate trait (flag-leaf height). Please confirm that the length of the tallest tiller is right.
>
> **209.** TODO - the definition of grain weight. Neither cited paper defines it, so the text gives no definition. Our old draft (May 2026) says the cleaned grain of the five harvested spikes was weighed per plant, but in the data the grains were counted on three spikes per plant (two for 13 plants in 2019-20), so it is not clear how many spikes were weighed. Please give the definition: what was weighed, on how many spikes, and in which unit.
>
> **210.** TODO - the definition of grain number. The text follows the cited paper, Potapenko et al. 2026b (Mol Ecol): number of grains per spike. The values used in our analyses, however, are the total number of grains on the spikes counted per plant (three spikes; two for 13 plants in 2019-20); grains per spike is a separate column that we did not use, and the Mol Ecol figures themselves label the trait "Number of grains". Please confirm which it is. If it is the total, the text should say "number of grains on three spikes per plant".
>
> **211.** TODO - how the grain was threshed. The text says by hand in 2019-20 and with a Haldrup LT-21 thresher in 2020-21 and 2021-22. N. Pintel's thesis, which covers only the first two seasons, says by hand in the first season and with the Haldrup LT-21 in the second. For 2021-22 the only source is our old draft (May 2026). Please confirm the threshing in 2021-22.
>
> **213.** TODO - was one NIR calibration used for all three seasons? The text describes one calibration: the manufacturer's calibration for cultivated barley, adjusted to the laboratory values for protein and starch, plus a new beta-glucan calibration built from the kit values. It does not say whether this calibration was applied to all three seasons; only our old draft (June 2026) says so, and N. Pintel's thesis covers the first two seasons only. N. Pintel's files show the calibration step itself: the 30 calibration accessions were scanned on 7 January 2020 (five scans each) with the factory calibrations for cultivated barley, before any adjustment and with no beta-glucan calibration. The scans of the seasons (June 2020 for 2019-20, October 2021 for 2021-22) do not record which calibration was used. Two signs that the seasons may differ: the season means of starch as recorded are 33.0, 42.9 and 20.7% (2019-20, 2020-21, 2021-22), and in 2021-22, 25 of 1,163 plants have negative starch values as recorded (down to -21.6%). The analysis scripts add 8.85 to the 2021-22 starch values, an adjustment carried over from an earlier pipeline that is not described in the text; it does not change the results, because each season is centred on its own mean. Please check whether one calibration was used for all seasons, and how it was adjusted.
>
> **234.** TODO - data and code availability. This paragraph keeps only the accession number of the raw reads, which TAG requires at the end of the Materials and methods. The full Data availability statement goes into the Declarations (not written yet) and still needs: (1) the MorexV3 SNP calls on Zenodo with a DOI (one VCF file of 6.9 GB, 300 accessions and 7,110,996 SNPs; decide whether to deposit all 300 accessions or the 290 analysed); (2) the analysis code in a GitHub repository of the lab (as for the two Potapenko et al. 2026 papers), archived on Zenodo for a DOI. Already done: the genotype BLUPs (the phenotypes) and the ENA identifiers of every accession are in Online Resource 1.
>
> **235.** TODO - statement on sharing the seeds. TAG asks that any restriction on sharing the material, here the seeds of the 300 accessions, be stated in the Materials and methods and in the cover letter; TAG expects materials to be freely available for non-commercial use within 60 days of a request. No document of the project states the lab's policy. Please give the sentence to write here, for example that the seeds are available from the corresponding author on request, or the restriction that applies.
>
> **236.** TODO - statement on the use of generative AI. TAG, following the Springer Nature policy, requires the use of large language models to be described in the Methods; only AI-assisted copy editing of text written by the authors is exempt. This statement was drafted on 30 September 2026 from the list of uses given by Shahar Liviatan. It names the tools (Claude Opus and Claude Sonnet from Anthropic, used through Claude Code; ChatGPT with GPT-5 and GPT-5.5 from OpenAI), the period (February to October 2026) and the tasks, and states that the authors checked all code, outputs and text and take responsibility for the content. Please approve it, or add the exact model versions or any other AI tool used in the project.
>
> **401.** TODO - the wild-versus-cultivated comparison rests on one old paper. The sentence says that domestication and breeding moved cultivated barley toward higher starch and lower protein and fiber, citing only Friedman and Atsmon (1988), which could not be checked (no PDF or abstract available; see the reference comment in the Introduction). Our experiment grew two cultivated checks, Morex and Clipper, but they are not in the analysed data. N. Pintel's NIR files contain one scan of each, from 17 November 2020: Morex 7.0% protein, 64.7% starch, 4.1% beta-glucan and 2.7% fiber; Clipper 12.9%, 59.7%, 5.4% and 4.2%; the wild accessions averaged about 16% protein and 33-43% starch in the first two seasons. One scan per cultivar is too little to replace the citation, but it points the same way. Decide whether to keep the citation alone, add these values, or find a second source.
>
> **402.** TODO - optional: mention how protein was measured where the Discussion explains why no protein gene was found. The protein results are the weakest: two loci and no gene that could be tested. One possible factor is the measurement: protein was predicted by NIR, calibrated against Kjeldahl values of 30 accessions (r = 0.92 on the same grain, 0.56 against the 2019-20 values of the same accessions; Online Resource 2). On 25 September 2026 it was decided to give no measurement detail here. Keep that decision and delete this comment, or add one sentence.
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** Approve as written, or tell me what to change (item and comment number).
>

---

## Change 1 of 125 · Checklist · heading

*Why:* Date of this update.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">**⚠ BEFORE SUBMISSION: open items (TAG requirements; checked 2026-10-01). Delete this list before submitting.**</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">**⚠ BEFORE SUBMISSION: open items (TAG requirements; updated 2026-10-01). Delete this list before submitting.**</ins>

**Reads after approval:**

> **⚠ BEFORE SUBMISSION: open items (TAG requirements; updated 2026-10-01). Delete this list before submitting.**

---

## Change 2 of 125 · Checklist · item 1

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">1. **Title page (not written yet):** the author names and their affiliations (institution, department, city, country); the corresponding author's e-mail and ORCID (ORCID is mandatory); and the Acknowledgments, which TAG places on the title page, with the funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, affiliation and e-mail also go into the title block of every Online Resource file, where they are placeholders now; the files are then rebuilt.</ins>
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page (not written yet):** the author names and their affiliations (institution, department, city, country); the corresponding author's e-mail and ORCID (ORCID is mandatory); and the Acknowledgments, which TAG places on the title page, with the funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, affiliation and e-mail also go into the title block of every Online Resource file, where they are placeholders now; the files are then rebuilt.
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 3 of 125 · Checklist · item 2

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
<del style="color:#c0392b;background:#fdecea">2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">2. **Author contribution statement (not written yet):** one short line per author, using initials. TAG prints it before the Acknowledgments, and it is also entered in the submission system.</ins>
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement (not written yet):** one short line per author, using initials. TAG prints it before the Acknowledgments, and it is also entered in the submission system.
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 4 of 125 · Checklist · item 3

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
<del style="color:#c0392b;background:#fdecea">3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3. **Statements and Declarations (not written yet; they go after the References):** TAG returns a manuscript that lacks any of them. Needed: Funding (agencies and grant numbers); Competing interests; Author contributions; Data availability; Ethics approval (not applicable to this study). The Data availability text is already drafted and names Online Resource 1 for the genotype BLUPs; it still needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub link (with its Zenodo archive) of the analysis code (comment 234).</ins>
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations (not written yet; they go after the References):** TAG returns a manuscript that lacks any of them. Needed: Funding (agencies and grant numbers); Competing interests; Author contributions; Data availability; Ethics approval (not applicable to this study). The Data availability text is already drafted and names Online Resource 1 for the genotype BLUPs; it still needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub link (with its Zenodo archive) of the analysis code (comment 234).
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 5 of 125 · Checklist · item 5

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
<del style="color:#c0392b;background:#fdecea">5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">5. **Online Resources (supplementary): done.** The 13 files are built, numbered in the order they are first cited, cited in the text, and captioned in the last section ("Supplementary Information"). Left to do: put the author names, affiliation and e-mail into the files (item 1); at submission, upload the 13 files and enter each one's caption in the submission system; and decide whether the extra, uncited supplement of raw-phenotype checks should join the supplement (it would then need a number and a citation in the text).</ins>
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary): done.** The 13 files are built, numbered in the order they are first cited, cited in the text, and captioned in the last section ("Supplementary Information"). Left to do: put the author names, affiliation and e-mail into the files (item 1); at submission, upload the 13 files and enter each one's caption in the submission system; and decide whether the extra, uncited supplement of raw-phenotype checks should join the supplement (it would then need a number and a citation in the text).
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 6 of 125 · Checklist · item 6

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
<del style="color:#c0392b;background:#fdecea">6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">6. **Other open comments** (each marked by a Word comment).</ins>
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Other open comments** (each marked by a Word comment).
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 7 of 125 · Checklist · item 7

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
<del style="color:#c0392b;background:#fdecea">7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">7. **Cover letter (not written yet):** give the accession number of the raw reads (ENA PRJEB79623); state any restriction on sharing the seeds of the accessions (the same statement belongs at the end of the Materials and methods); and suggest reviewers who are independent of the authors, from different countries, with institutional e-mail addresses.</ins>
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter (not written yet):** give the accession number of the raw reads (ENA PRJEB79623); state any restriction on sharing the seeds of the accessions (the same statement belongs at the end of the Materials and methods); and suggest reviewers who are independent of the authors, from different countries, with institutional e-mail addresses.
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 8 of 125 · Checklist · item 8

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
<del style="color:#c0392b;background:#fdecea">8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">8. **Final Word file:** 10-pt Times Roman; automatic page numbering; no field functions. The four figures are embedded in the text and are also uploaded as separate TIFF files, Fig1 to Fig4 (600 dpi, 174 mm wide). Before submitting, delete this list and every remaining Word comment.</ins>
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word file:** 10-pt Times Roman; automatic page numbering; no field functions. The four figures are embedded in the text and are also uploaded as separate TIFF files, Fig1 to Fig4 (600 dpi, 174 mm wide). Before submitting, delete this list and every remaining Word comment.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 9 of 125 · Checklist · item 9

*Why:* Updated to the current state, plain status line like item 4, no paths (user, 2026-10-01).

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
<del style="color:#c0392b;background:#fdecea">9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">9. *(Optional)* **Figure accessibility:** TAG recommends patterns or symbols in addition to colour, but Figs. 2–4 tell loci and haplotype groups apart by colour alone.</ins>

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG recommends patterns or symbols in addition to colour, but Figs. 2–4 tell loci and haplotype groups apart by colour alone.

---

## Change 10 of 125 · Checklist · hidden note

*Why:* Record (invisible in Word).

**With the changes marked:**



**Reads after approval:**

> *(paragraph removed; nothing is left of it in Word)*

---

## Change 11 of 125 · Delete your comment 311

*Why:* Resolved: answered in the text (Results Ch. 2: the 4H cluster is described as two adjacent loci split by the 50-kb gap rule).

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), <del style="color:#c0392b;background:#fdecea"></del>which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 12 of 125 · Delete your comment 311 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><del style="color:#c0392b;background:#fdecea"><sup>(comment 311)</sup></del><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 13 of 125 · Delete your comment 312

*Why:* Resolved: same as 311.

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), <del style="color:#c0392b;background:#fdecea"></del>which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 14 of 125 · Delete your comment 312 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<del style="color:#c0392b;background:#fdecea"><sup>(comment 312)</sup></del><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 15 of 125 · Delete your comment 313

*Why:* Resolved: answered in the text (Results Ch. 2: six fiber loci on 1H, five of them in two clusters).

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). <del style="color:#c0392b;background:#fdecea"></del>Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 16 of 125 · Delete your comment 313 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <del style="color:#c0392b;background:#fdecea"><sup>(comment 313)</sup></del>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb.

---

## Change 17 of 125 · Delete your comment 314

*Why:* Resolved: answered: the 52 significant SNPs with MAF and effects are in Online Resource 5, cited in this paragraph.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea"></del>Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

**Reads after approval:**

> Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

---

## Change 18 of 125 · Delete your comment 314 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <del style="color:#c0392b;background:#fdecea"><sup>(comment 314)</sup></del>

**Reads after approval:**

> Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839).

---

## Change 19 of 125 · Delete your comment 315

*Why:* Resolved: same as 314.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea"></del>Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

**Reads after approval:**

> Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

---

## Change 20 of 125 · Delete your comment 315 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <del style="color:#c0392b;background:#fdecea"><sup>(comment 315)</sup></del>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

**Reads after approval:**

> Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

---

## Change 21 of 125 · Delete your comment 316

*Why:* Resolved: answered: all locus coordinates are in Online Resource 5.

**With the changes marked:**

To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. <del style="color:#c0392b;background:#fdecea"></del>Fiber and starch loci co-localized in three genomic regions<sup>(comment 316)</sup><sup>(comment 332)</sup>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

**Reads after approval:**

> To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. Fiber and starch loci co-localized in three genomic regions<sup>(comment 316)</sup><sup>(comment 332)</sup>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

---

## Change 22 of 125 · Delete your comment 316 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. Fiber and starch loci co-localized in three genomic regions<del style="color:#c0392b;background:#fdecea"><sup>(comment 316)</sup></del><sup>(comment 332)</sup>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

**Reads after approval:**

> To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. Fiber and starch loci co-localized in three genomic regions<sup>(comment 332)</sup>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

---

## Change 23 of 125 · Delete your comment 317

*Why:* Resolved: answered: all 55 genes are in Online Resource 6, cited in this paragraph.

**With the changes marked:**

Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. <del style="color:#c0392b;background:#fdecea"></del>Of the 55 candidate genes, 30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6)<sup>(comment 333)</sup>. The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested. <sup>(comment 317)</sup>

**Reads after approval:**

> Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. Of the 55 candidate genes, 30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6)<sup>(comment 333)</sup>. The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested. <sup>(comment 317)</sup>

---

## Change 24 of 125 · Delete your comment 317 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. Of the 55 candidate genes, 30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6)<sup>(comment 333)</sup>. The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested. <del style="color:#c0392b;background:#fdecea"><sup>(comment 317)</sup></del>

**Reads after approval:**

> Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. Of the 55 candidate genes, 30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6)<sup>(comment 333)</sup>. The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested.

---

## Change 25 of 125 · Delete your comment 318

*Why:* Resolved: answered: Online Resource 6.

**With the changes marked:**

A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (<del style="color:#c0392b;background:#fdecea"></del>Online Resource 6<sup>(comment 318)</sup><sup>(comment 334)</sup>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

**Reads after approval:**

> A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (Online Resource 6<sup>(comment 318)</sup><sup>(comment 334)</sup>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

---

## Change 26 of 125 · Delete your comment 318 (end mark)

*Why:* Its end mark.

**With the changes marked:**

A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (Online Resource 6<del style="color:#c0392b;background:#fdecea"><sup>(comment 318)</sup></del><sup>(comment 334)</sup>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

**Reads after approval:**

> A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (Online Resource 6<sup>(comment 334)</sup>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

---

## Change 27 of 125 · Delete your comment 319

*Why:* Resolved: answered: the compared and excluded SNP sites are in Online Resource 10.

**With the changes marked:**

To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 10<sup>(comment 319)</sup><sup>(comment 338)</sup>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

**Reads after approval:**

> To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; Online Resource 10<sup>(comment 319)</sup><sup>(comment 338)</sup>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

---

## Change 28 of 125 · Delete your comment 319 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; Online Resource 10<del style="color:#c0392b;background:#fdecea"><sup>(comment 319)</sup></del><sup>(comment 338)</sup>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

**Reads after approval:**

> To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; Online Resource 10<sup>(comment 338)</sup>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

---

## Change 29 of 125 · Delete your comment 414

*Why:* Resolved: answered: the distances are in Online Resource 8, cited in this sentence.

**With the changes marked:**

β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, <del style="color:#c0392b;background:#fdecea"></del>lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<sup>(comment 414)</sup><sup>(comment 426)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

**Reads after approval:**

> β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<sup>(comment 414)</sup><sup>(comment 426)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

---

## Change 30 of 125 · Delete your comment 414 (end mark)

*Why:* Its end mark.

**With the changes marked:**

β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<del style="color:#c0392b;background:#fdecea"><sup>(comment 414)</sup></del><sup>(comment 426)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

**Reads after approval:**

> β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<sup>(comment 426)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

---

## Change 31 of 125 · Delete your comment 415

*Why:* Resolved: answered: the sentence refers to the Materials and methods.

**With the changes marked:**

Ranking genes by the significance of their haplotype test was a poor guide to biological candidacy. The genes with the strongest haplotype associations were often the least interpretable, unannotated or derived from mobile elements, whereas the genes with a clear route to the trait ranked lower. Because genes within a locus share one linkage-disequilibrium block, their significance largely reflects linkage with the lead signal rather than their own function. Candidate genes were therefore selected in three steps: the genome-wide association identified the locus, the haplotype test retained the genes whose variation was associated with the trait, and biological relevance chose among these genes. The first two steps can also miss genuine candidates. The definition of loci delimits each association by linkage disequilibrium with its lead SNP, so a gene lying just outside a narrow locus is never considered. The GDSL esterase/lipase at 7H illustrates this: it lay beside two of the strongest signals of the study and was found only by examining the region directly. The haplotype analysis requires sufficient variation within a gene to define haplotype groups, so a gene that could not be tested is not a gene without effect. Likewise, one fixed grouping setting for all genes can leave the associated SNPs of a sparse gene ungrouped, so a real association can be missed even in a gene that was tested. <del style="color:#c0392b;background:#fdecea"></del>This was also the case at the same gene, whose associated SNPs were grouped only at a different setting (see Materials and methods).<sup>(comment 415)</sup> The three genes carried forward are therefore those that passed both filters, not a complete list of the genes underlying these traits.

**Reads after approval:**

> Ranking genes by the significance of their haplotype test was a poor guide to biological candidacy. The genes with the strongest haplotype associations were often the least interpretable, unannotated or derived from mobile elements, whereas the genes with a clear route to the trait ranked lower. Because genes within a locus share one linkage-disequilibrium block, their significance largely reflects linkage with the lead signal rather than their own function. Candidate genes were therefore selected in three steps: the genome-wide association identified the locus, the haplotype test retained the genes whose variation was associated with the trait, and biological relevance chose among these genes. The first two steps can also miss genuine candidates. The definition of loci delimits each association by linkage disequilibrium with its lead SNP, so a gene lying just outside a narrow locus is never considered. The GDSL esterase/lipase at 7H illustrates this: it lay beside two of the strongest signals of the study and was found only by examining the region directly. The haplotype analysis requires sufficient variation within a gene to define haplotype groups, so a gene that could not be tested is not a gene without effect. Likewise, one fixed grouping setting for all genes can leave the associated SNPs of a sparse gene ungrouped, so a real association can be missed even in a gene that was tested. This was also the case at the same gene, whose associated SNPs were grouped only at a different setting (see Materials and methods).<sup>(comment 415)</sup> The three genes carried forward are therefore those that passed both filters, not a complete list of the genes underlying these traits.

---

## Change 32 of 125 · Delete your comment 415 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Ranking genes by the significance of their haplotype test was a poor guide to biological candidacy. The genes with the strongest haplotype associations were often the least interpretable, unannotated or derived from mobile elements, whereas the genes with a clear route to the trait ranked lower. Because genes within a locus share one linkage-disequilibrium block, their significance largely reflects linkage with the lead signal rather than their own function. Candidate genes were therefore selected in three steps: the genome-wide association identified the locus, the haplotype test retained the genes whose variation was associated with the trait, and biological relevance chose among these genes. The first two steps can also miss genuine candidates. The definition of loci delimits each association by linkage disequilibrium with its lead SNP, so a gene lying just outside a narrow locus is never considered. The GDSL esterase/lipase at 7H illustrates this: it lay beside two of the strongest signals of the study and was found only by examining the region directly. The haplotype analysis requires sufficient variation within a gene to define haplotype groups, so a gene that could not be tested is not a gene without effect. Likewise, one fixed grouping setting for all genes can leave the associated SNPs of a sparse gene ungrouped, so a real association can be missed even in a gene that was tested. This was also the case at the same gene, whose associated SNPs were grouped only at a different setting (see Materials and methods).<del style="color:#c0392b;background:#fdecea"><sup>(comment 415)</sup></del> The three genes carried forward are therefore those that passed both filters, not a complete list of the genes underlying these traits.

**Reads after approval:**

> Ranking genes by the significance of their haplotype test was a poor guide to biological candidacy. The genes with the strongest haplotype associations were often the least interpretable, unannotated or derived from mobile elements, whereas the genes with a clear route to the trait ranked lower. Because genes within a locus share one linkage-disequilibrium block, their significance largely reflects linkage with the lead signal rather than their own function. Candidate genes were therefore selected in three steps: the genome-wide association identified the locus, the haplotype test retained the genes whose variation was associated with the trait, and biological relevance chose among these genes. The first two steps can also miss genuine candidates. The definition of loci delimits each association by linkage disequilibrium with its lead SNP, so a gene lying just outside a narrow locus is never considered. The GDSL esterase/lipase at 7H illustrates this: it lay beside two of the strongest signals of the study and was found only by examining the region directly. The haplotype analysis requires sufficient variation within a gene to define haplotype groups, so a gene that could not be tested is not a gene without effect. Likewise, one fixed grouping setting for all genes can leave the associated SNPs of a sparse gene ungrouped, so a real association can be missed even in a gene that was tested. This was also the case at the same gene, whose associated SNPs were grouped only at a different setting (see Materials and methods). The three genes carried forward are therefore those that passed both filters, not a complete list of the genes underlying these traits.

---

## Change 33 of 125 · Delete your comment 416

*Why:* Resolved: answered: the References session checked every citation of this paragraph (all supported).

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea"></del>Each candidate gene connects its trait to a specific process in the developing grain. The acyltransferase *GPAT6* suggests that cutin deposition in the pericarp and testa contributes to the insoluble fiber fraction, consistent with the central role of GPAT6 in fruit cutin biosynthesis in tomato (Petit et al. 2016). The glucanase *GH17* suggests that the hydrolysis of cell-wall glucans contributes to grain fiber, consistent with the role of its family, which includes barley’s own (1,3;1,4)-β-glucan endohydrolases, in degrading mixed-linkage glucan (Hrmova et al. 2024). The phosphate transporter *PHT4;3* suggests that plastid phosphate homeostasis contributes to the metabolic conditions required for starch accumulation, consistent with the impaired grain filling reported in rice phosphate-transporter mutants (Ma et al. 2021). The GDSL esterase/lipase at 7H suggests that the modification of cell-wall esters contributes to grain fiber, consistent with the roles of this family in arabinoxylan deacetylation and cutin deposition (Girard et al. 2012; Zhang et al. 2019). Its associated SNPs were confined to the putative promoter, while variants within the gene body showed no association, pointing to regulatory rather than coding variation. These four genes therefore identify distinct biological pathways for future investigation: cutin deposition, glucan hydrolysis, plastid phosphate transport and cell-wall ester modification. <sup>(comment 416)</sup>

**Reads after approval:**

> Each candidate gene connects its trait to a specific process in the developing grain. The acyltransferase *GPAT6* suggests that cutin deposition in the pericarp and testa contributes to the insoluble fiber fraction, consistent with the central role of GPAT6 in fruit cutin biosynthesis in tomato (Petit et al. 2016). The glucanase *GH17* suggests that the hydrolysis of cell-wall glucans contributes to grain fiber, consistent with the role of its family, which includes barley’s own (1,3;1,4)-β-glucan endohydrolases, in degrading mixed-linkage glucan (Hrmova et al. 2024). The phosphate transporter *PHT4;3* suggests that plastid phosphate homeostasis contributes to the metabolic conditions required for starch accumulation, consistent with the impaired grain filling reported in rice phosphate-transporter mutants (Ma et al. 2021). The GDSL esterase/lipase at 7H suggests that the modification of cell-wall esters contributes to grain fiber, consistent with the roles of this family in arabinoxylan deacetylation and cutin deposition (Girard et al. 2012; Zhang et al. 2019). Its associated SNPs were confined to the putative promoter, while variants within the gene body showed no association, pointing to regulatory rather than coding variation. These four genes therefore identify distinct biological pathways for future investigation: cutin deposition, glucan hydrolysis, plastid phosphate transport and cell-wall ester modification. <sup>(comment 416)</sup>

---

## Change 34 of 125 · Delete your comment 416 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Each candidate gene connects its trait to a specific process in the developing grain. The acyltransferase *GPAT6* suggests that cutin deposition in the pericarp and testa contributes to the insoluble fiber fraction, consistent with the central role of GPAT6 in fruit cutin biosynthesis in tomato (Petit et al. 2016). The glucanase *GH17* suggests that the hydrolysis of cell-wall glucans contributes to grain fiber, consistent with the role of its family, which includes barley’s own (1,3;1,4)-β-glucan endohydrolases, in degrading mixed-linkage glucan (Hrmova et al. 2024). The phosphate transporter *PHT4;3* suggests that plastid phosphate homeostasis contributes to the metabolic conditions required for starch accumulation, consistent with the impaired grain filling reported in rice phosphate-transporter mutants (Ma et al. 2021). The GDSL esterase/lipase at 7H suggests that the modification of cell-wall esters contributes to grain fiber, consistent with the roles of this family in arabinoxylan deacetylation and cutin deposition (Girard et al. 2012; Zhang et al. 2019). Its associated SNPs were confined to the putative promoter, while variants within the gene body showed no association, pointing to regulatory rather than coding variation. These four genes therefore identify distinct biological pathways for future investigation: cutin deposition, glucan hydrolysis, plastid phosphate transport and cell-wall ester modification. <del style="color:#c0392b;background:#fdecea"><sup>(comment 416)</sup></del>

**Reads after approval:**

> Each candidate gene connects its trait to a specific process in the developing grain. The acyltransferase *GPAT6* suggests that cutin deposition in the pericarp and testa contributes to the insoluble fiber fraction, consistent with the central role of GPAT6 in fruit cutin biosynthesis in tomato (Petit et al. 2016). The glucanase *GH17* suggests that the hydrolysis of cell-wall glucans contributes to grain fiber, consistent with the role of its family, which includes barley’s own (1,3;1,4)-β-glucan endohydrolases, in degrading mixed-linkage glucan (Hrmova et al. 2024). The phosphate transporter *PHT4;3* suggests that plastid phosphate homeostasis contributes to the metabolic conditions required for starch accumulation, consistent with the impaired grain filling reported in rice phosphate-transporter mutants (Ma et al. 2021). The GDSL esterase/lipase at 7H suggests that the modification of cell-wall esters contributes to grain fiber, consistent with the roles of this family in arabinoxylan deacetylation and cutin deposition (Girard et al. 2012; Zhang et al. 2019). Its associated SNPs were confined to the putative promoter, while variants within the gene body showed no association, pointing to regulatory rather than coding variation. These four genes therefore identify distinct biological pathways for future investigation: cutin deposition, glucan hydrolysis, plastid phosphate transport and cell-wall ester modification.

---

## Change 35 of 125 · Delete TODO 203

*Why:* Resolved: the two GATK citations that replaced the preprint passed the References session's bibliographic check (tool citations).

**With the changes marked:**

Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026b). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed as described by Potapenko et al. (2026b)<sup>(TODO 202)</sup>. Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow <del style="color:#c0392b;background:#fdecea">(McKenna et al. 2010; Van der Auwera and O'Connor 2020)<sup>(TODO 203)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">(McKenna et al. 2010; Van der Auwera and O'Connor 2020)</ins>. For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, with at most two alternate alleles per site. To reduce the computational load of the large genome, calling was run separately in non-overlapping 60-Mb intervals.

**Reads after approval:**

> Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026b). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed as described by Potapenko et al. (2026b)<sup>(TODO 202)</sup>. Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow (McKenna et al. 2010; Van der Auwera and O'Connor 2020). For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, with at most two alternate alleles per site. To reduce the computational load of the large genome, calling was run separately in non-overlapping 60-Mb intervals.

---

## Change 36 of 125 · Delete supmat comment 239

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 1)<sup>(comment 239)</sup>. At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through three rounds of single-seed descent and selfing<sup>(TODO 201)</sup>.

**Reads after approval:**

> We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil (Online Resource 1)<sup>(comment 239)</sup>. At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through three rounds of single-seed descent and selfing<sup>(TODO 201)</sup>.

---

## Change 37 of 125 · Delete supmat comment 239 (end mark)

*Why:* Its end mark.

**With the changes marked:**

We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil (Online Resource 1)<del style="color:#c0392b;background:#fdecea"><sup>(comment 239)</sup></del>. At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through three rounds of single-seed descent and selfing<sup>(TODO 201)</sup>.

**Reads after approval:**

> We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil (Online Resource 1). At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through three rounds of single-seed descent and selfing<sup>(TODO 201)</sup>.

---

## Change 38 of 125 · Delete supmat comment 240

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites. Their grain was scanned five times and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), with the protocols adapted to small samples; each accession was assayed in four subsamples of about 200 mg, and their mean was used. Total protein was determined by the Kjeldahl method (nitrogen × 6.25) in three subsamples of about 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College. The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.<sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (<del style="color:#c0392b;background:#fdecea"></del>Online Resource 2<sup>(comment 240)</sup>).

**Reads after approval:**

> The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites. Their grain was scanned five times and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), with the protocols adapted to small samples; each accession was assayed in four subsamples of about 200 mg, and their mean was used. Total protein was determined by the Kjeldahl method (nitrogen × 6.25) in three subsamples of about 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College. The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.<sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (Online Resource 2<sup>(comment 240)</sup>).

---

## Change 39 of 125 · Delete supmat comment 240 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites. Their grain was scanned five times and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), with the protocols adapted to small samples; each accession was assayed in four subsamples of about 200 mg, and their mean was used. Total protein was determined by the Kjeldahl method (nitrogen × 6.25) in three subsamples of about 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College. The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.<sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (Online Resource 2<del style="color:#c0392b;background:#fdecea"><sup>(comment 240)</sup></del>).

**Reads after approval:**

> The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites. Their grain was scanned five times and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), with the protocols adapted to small samples; each accession was assayed in four subsamples of about 200 mg, and their mean was used. Total protein was determined by the Kjeldahl method (nitrogen × 6.25) in three subsamples of about 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College. The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.<sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (Online Resource 2).

---

## Change 40 of 125 · Delete supmat comment 241

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021). Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 2<sup>(comment 241)</sup>). The nutritional-trait data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

**Reads after approval:**

> All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021). Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource 2<sup>(comment 241)</sup>). The nutritional-trait data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

---

## Change 41 of 125 · Delete supmat comment 241 (end mark)

*Why:* Its end mark.

**With the changes marked:**

All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021). Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource 2<del style="color:#c0392b;background:#fdecea"><sup>(comment 241)</sup></del>). The nutritional-trait data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

**Reads after approval:**

> All phenotypic analyses were performed in R v4.1.2 (R Core Team 2021). Figures were drawn in R with ggplot2 v4.0.1 (Wickham 2016), patchwork v1.3.2 (Pedersen 2025) and base graphics. Block 4 of the 2019–20 season was excluded because its β-glucan values correlated poorly with those of the other blocks (Pearson r = 0.27–0.34, compared with 0.61–0.64 among blocks 1–3; Online Resource 2). The nutritional-trait data thus comprised 3,164 plants of the 290 accessions. Each trait was centered within season by subtracting the season mean. Missing values were then excluded, and values more than three standard deviations from the mean of the centered values, pooled over the three seasons, were removed as outliers.

---

## Change 42 of 125 · Delete supmat comment 242

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To estimate the genetic value of each accession, we fitted a linear mixed model for each of the ten traits (the four nutritional and six morphological and phenological traits) in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as random effects. The genotype best linear unbiased predictions (BLUPs) of this model <del style="color:#c0392b;background:#fdecea"></del>were used as the phenotypes in all downstream analyses (Online Resource 1)<sup>(comment 242)</sup>.

**Reads after approval:**

> To estimate the genetic value of each accession, we fitted a linear mixed model for each of the ten traits (the four nutritional and six morphological and phenological traits) in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as random effects. The genotype best linear unbiased predictions (BLUPs) of this model were used as the phenotypes in all downstream analyses (Online Resource 1)<sup>(comment 242)</sup>.

---

## Change 43 of 125 · Delete supmat comment 242 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To estimate the genetic value of each accession, we fitted a linear mixed model for each of the ten traits (the four nutritional and six morphological and phenological traits) in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as random effects. The genotype best linear unbiased predictions (BLUPs) of this model were used as the phenotypes in all downstream analyses (Online Resource 1)<del style="color:#c0392b;background:#fdecea"><sup>(comment 242)</sup></del>.

**Reads after approval:**

> To estimate the genetic value of each accession, we fitted a linear mixed model for each of the ten traits (the four nutritional and six morphological and phenological traits) in lme4 v1.1.35.3 (Bates et al. 2015), with genotype and season-by-block as random effects. The genotype best linear unbiased predictions (BLUPs) of this model were used as the phenotypes in all downstream analyses (Online Resource 1).

---

## Change 44 of 125 · Delete supmat comment 243

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

This model was fitted for the four nutritional traits and for flowering time, number of tillers, grain weight and spike length (Fig. 1a; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 3<sup>(comment 243)</sup>). To visualize genotype-by-season interaction, we drew reaction norms from the mean of the centered values of each accession in each season (Fig. 1b).

**Reads after approval:**

> This model was fitted for the four nutritional traits and for flowering time, number of tillers, grain weight and spike length (Fig. 1a; Online Resource 3<sup>(comment 243)</sup>). To visualize genotype-by-season interaction, we drew reaction norms from the mean of the centered values of each accession in each season (Fig. 1b).

---

## Change 45 of 125 · Delete supmat comment 243 (end mark)

*Why:* Its end mark.

**With the changes marked:**

This model was fitted for the four nutritional traits and for flowering time, number of tillers, grain weight and spike length (Fig. 1a; Online Resource 3<del style="color:#c0392b;background:#fdecea"><sup>(comment 243)</sup></del>). To visualize genotype-by-season interaction, we drew reaction norms from the mean of the centered values of each accession in each season (Fig. 1b).

**Reads after approval:**

> This model was fitted for the four nutritional traits and for flowering time, number of tillers, grain weight and spike length (Fig. 1a; Online Resource 3). To visualize genotype-by-season interaction, we drew reaction norms from the mean of the centered values of each accession in each season (Fig. 1b).

---

## Change 46 of 125 · Delete supmat comment 244

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To test whether the nutritional traits varied with the environment of the sites of origin, we averaged the BLUPs of the ten accessions of each of the 29 sites for each trait. Ten environmental variables were available for each site: maximum temperature in March, average annual precipitation and elevation from the Israeli Meteorological Service, and seven soil properties measured in a soil sample from each site: pH, electrical conductivity, clay, silt, sand, organic carbon and total nitrogen<sup>(TODO 216)</sup> (Potapenko et al. 2026a). Two pairs of variables were strongly correlated across sites (|r| > 0.8), elevation with March temperature and silt with sand, and elevation and silt were removed <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 1)<sup>(comment 244)</sup>. The site means were then correlated with the eight remaining variables by two-sided Pearson correlation (Fig. 1d), with Spearman correlation as a sensitivity analysis.

**Reads after approval:**

> To test whether the nutritional traits varied with the environment of the sites of origin, we averaged the BLUPs of the ten accessions of each of the 29 sites for each trait. Ten environmental variables were available for each site: maximum temperature in March, average annual precipitation and elevation from the Israeli Meteorological Service, and seven soil properties measured in a soil sample from each site: pH, electrical conductivity, clay, silt, sand, organic carbon and total nitrogen<sup>(TODO 216)</sup> (Potapenko et al. 2026a). Two pairs of variables were strongly correlated across sites (|r| > 0.8), elevation with March temperature and silt with sand, and elevation and silt were removed (Online Resource 1)<sup>(comment 244)</sup>. The site means were then correlated with the eight remaining variables by two-sided Pearson correlation (Fig. 1d), with Spearman correlation as a sensitivity analysis.

---

## Change 47 of 125 · Delete supmat comment 244 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To test whether the nutritional traits varied with the environment of the sites of origin, we averaged the BLUPs of the ten accessions of each of the 29 sites for each trait. Ten environmental variables were available for each site: maximum temperature in March, average annual precipitation and elevation from the Israeli Meteorological Service, and seven soil properties measured in a soil sample from each site: pH, electrical conductivity, clay, silt, sand, organic carbon and total nitrogen<sup>(TODO 216)</sup> (Potapenko et al. 2026a). Two pairs of variables were strongly correlated across sites (|r| > 0.8), elevation with March temperature and silt with sand, and elevation and silt were removed (Online Resource 1)<del style="color:#c0392b;background:#fdecea"><sup>(comment 244)</sup></del>. The site means were then correlated with the eight remaining variables by two-sided Pearson correlation (Fig. 1d), with Spearman correlation as a sensitivity analysis.

**Reads after approval:**

> To test whether the nutritional traits varied with the environment of the sites of origin, we averaged the BLUPs of the ten accessions of each of the 29 sites for each trait. Ten environmental variables were available for each site: maximum temperature in March, average annual precipitation and elevation from the Israeli Meteorological Service, and seven soil properties measured in a soil sample from each site: pH, electrical conductivity, clay, silt, sand, organic carbon and total nitrogen<sup>(TODO 216)</sup> (Potapenko et al. 2026a). Two pairs of variables were strongly correlated across sites (|r| > 0.8), elevation with March temperature and silt with sand, and elevation and silt were removed (Online Resource 1). The site means were then correlated with the eight remaining variables by two-sided Pearson correlation (Fig. 1d), with Spearman correlation as a sensitivity analysis.

---

## Change 48 of 125 · Delete supmat comment 245

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Genome-wide LD decay was estimated from the pairwise r² of all SNP pairs within 2.5 Mb, computed in PLINK, averaged in 1-kb distance bins and smoothed with a LOESS curve weighted by the number of pairs per bin <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 4)<sup>(comment 245)</sup>.

**Reads after approval:**

> Genome-wide LD decay was estimated from the pairwise r² of all SNP pairs within 2.5 Mb, computed in PLINK, averaged in 1-kb distance bins and smoothed with a LOESS curve weighted by the number of pairs per bin (Online Resource 4)<sup>(comment 245)</sup>.

---

## Change 49 of 125 · Delete supmat comment 245 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Genome-wide LD decay was estimated from the pairwise r² of all SNP pairs within 2.5 Mb, computed in PLINK, averaged in 1-kb distance bins and smoothed with a LOESS curve weighted by the number of pairs per bin (Online Resource 4)<del style="color:#c0392b;background:#fdecea"><sup>(comment 245)</sup></del>.

**Reads after approval:**

> Genome-wide LD decay was estimated from the pairwise r² of all SNP pairs within 2.5 Mb, computed in PLINK, averaged in 1-kb distance bins and smoothed with a LOESS curve weighted by the number of pairs per bin (Online Resource 4).

---

## Change 50 of 125 · Delete supmat comment 246

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To define loci, we grouped the significant SNPs of each trait by LD clumping in PLINK (Purcell et al. 2007). Only significant SNPs could serve as lead SNPs, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs. Significant SNPs left outside every locus then served as lead SNPs in further rounds of clumping, until every significant SNP lay within a locus, giving 36 non-overlapping loci <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 5)<sup>(comment 246)</sup>.

**Reads after approval:**

> To define loci, we grouped the significant SNPs of each trait by LD clumping in PLINK (Purcell et al. 2007). Only significant SNPs could serve as lead SNPs, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs. Significant SNPs left outside every locus then served as lead SNPs in further rounds of clumping, until every significant SNP lay within a locus, giving 36 non-overlapping loci (Online Resource 5)<sup>(comment 246)</sup>.

---

## Change 51 of 125 · Delete supmat comment 246 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To define loci, we grouped the significant SNPs of each trait by LD clumping in PLINK (Purcell et al. 2007). Only significant SNPs could serve as lead SNPs, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs. Significant SNPs left outside every locus then served as lead SNPs in further rounds of clumping, until every significant SNP lay within a locus, giving 36 non-overlapping loci (Online Resource 5)<del style="color:#c0392b;background:#fdecea"><sup>(comment 246)</sup></del>.

**Reads after approval:**

> To define loci, we grouped the significant SNPs of each trait by LD clumping in PLINK (Purcell et al. 2007). Only significant SNPs could serve as lead SNPs, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs. Significant SNPs left outside every locus then served as lead SNPs in further rounds of clumping, until every significant SNP lay within a locus, giving 36 non-overlapping loci (Online Resource 5).

---

## Change 52 of 125 · Delete supmat comment 247

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To identify candidate genes, we intersected the interval of each locus, from its first to its last member SNP, with the MorexV3 high-confidence gene annotation (Mascher et al. 2021) from Ensembl Plants release 62 (Yates et al. 2026) using bedtools intersect v2.30.0 (Quinlan and Hall 2010). Every gene overlapping a locus interval was taken as a candidate gene <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 6)<sup>(comment 247)</sup>.

**Reads after approval:**

> To identify candidate genes, we intersected the interval of each locus, from its first to its last member SNP, with the MorexV3 high-confidence gene annotation (Mascher et al. 2021) from Ensembl Plants release 62 (Yates et al. 2026) using bedtools intersect v2.30.0 (Quinlan and Hall 2010). Every gene overlapping a locus interval was taken as a candidate gene (Online Resource 6)<sup>(comment 247)</sup>.

---

## Change 53 of 125 · Delete supmat comment 247 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To identify candidate genes, we intersected the interval of each locus, from its first to its last member SNP, with the MorexV3 high-confidence gene annotation (Mascher et al. 2021) from Ensembl Plants release 62 (Yates et al. 2026) using bedtools intersect v2.30.0 (Quinlan and Hall 2010). Every gene overlapping a locus interval was taken as a candidate gene (Online Resource 6)<del style="color:#c0392b;background:#fdecea"><sup>(comment 247)</sup></del>.

**Reads after approval:**

> To identify candidate genes, we intersected the interval of each locus, from its first to its last member SNP, with the MorexV3 high-confidence gene annotation (Mascher et al. 2021) from Ensembl Plants release 62 (Yates et al. 2026) using bedtools intersect v2.30.0 (Quinlan and Hall 2010). Every gene overlapping a locus interval was taken as a candidate gene (Online Resource 6).

---

## Change 54 of 125 · Delete supmat comment 248

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

All genes were analyzed with the same crosshap settings: a DBSCAN clustering radius (ε) of 0.6, a minimum marker-group size (MGmin) of two SNPs and a minimum haplotype size (minHap) of nine accessions, with all other settings at their defaults. MGmin = 2 allowed genes with few SNPs to be tested, and ε = 0.6 was chosen without reference to the phenotypes, as the value that balanced the number of testable genes against the proportion of accessions assigned to haplotype groups <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 7)<sup>(comment 248)</sup>.

**Reads after approval:**

> All genes were analyzed with the same crosshap settings: a DBSCAN clustering radius (ε) of 0.6, a minimum marker-group size (MGmin) of two SNPs and a minimum haplotype size (minHap) of nine accessions, with all other settings at their defaults. MGmin = 2 allowed genes with few SNPs to be tested, and ε = 0.6 was chosen without reference to the phenotypes, as the value that balanced the number of testable genes against the proportion of accessions assigned to haplotype groups (Online Resource 7)<sup>(comment 248)</sup>.

---

## Change 55 of 125 · Delete supmat comment 248 (end mark)

*Why:* Its end mark.

**With the changes marked:**

All genes were analyzed with the same crosshap settings: a DBSCAN clustering radius (ε) of 0.6, a minimum marker-group size (MGmin) of two SNPs and a minimum haplotype size (minHap) of nine accessions, with all other settings at their defaults. MGmin = 2 allowed genes with few SNPs to be tested, and ε = 0.6 was chosen without reference to the phenotypes, as the value that balanced the number of testable genes against the proportion of accessions assigned to haplotype groups (Online Resource 7)<del style="color:#c0392b;background:#fdecea"><sup>(comment 248)</sup></del>.

**Reads after approval:**

> All genes were analyzed with the same crosshap settings: a DBSCAN clustering radius (ε) of 0.6, a minimum marker-group size (MGmin) of two SNPs and a minimum haplotype size (minHap) of nine accessions, with all other settings at their defaults. MGmin = 2 allowed genes with few SNPs to be tested, and ε = 0.6 was chosen without reference to the phenotypes, as the value that balanced the number of testable genes against the proportion of accessions assigned to haplotype groups (Online Resource 7).

---

## Change 56 of 125 · Delete supmat comment 249

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 6)<sup>(comment 249)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014), and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 3).

**Reads after approval:**

> For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource 6)<sup>(comment 249)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014), and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 3).

---

## Change 57 of 125 · Delete supmat comment 249 (end mark)

*Why:* Its end mark.

**With the changes marked:**

For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource 6)<del style="color:#c0392b;background:#fdecea"><sup>(comment 249)</sup></del>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014), and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 3).

**Reads after approval:**

> For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource 6). *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014), and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 3).

---

## Change 58 of 125 · Delete supmat comment 250

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Only 3 of the 55 candidate genes had a functional description in the MorexV3 annotation. We therefore annotated the 23 genes with significant haplotype associations from the protein sequence of their longest isoform in the MorexV3 high-confidence proteome (Mascher et al. 2021). Each protein was first searched against UniProtKB/Swiss-Prot release 2026_03 (UniProt Consortium 2025) by BLASTP in DIAMOND v2.0.14 (Buchfink et al. 2021), and its best hit was accepted when the e-value was ≤ 10⁻⁵, the hit covered ≥ 50% of the query at ≥ 30% identity, and the hit was not described as uncharacterized, hypothetical, predicted or of unknown function. Protein domains and Gene Ontology terms were then assigned with InterProScan v5.78-109.0 (Jones et al. 2014). Proteins without a result from either source were searched with BLASTP (BLAST+ v2.12.0; Camacho et al. 2009) against the NCBI non-redundant protein database, restricted to green plants (Viridiplantae; accessed 9 September 2026), with the same criteria. Each gene took its annotation from the first of these sources that gave one, in the order Swiss-Prot, InterPro and nr <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 6)<sup>(comment 250)</sup>.

**Reads after approval:**

> Only 3 of the 55 candidate genes had a functional description in the MorexV3 annotation. We therefore annotated the 23 genes with significant haplotype associations from the protein sequence of their longest isoform in the MorexV3 high-confidence proteome (Mascher et al. 2021). Each protein was first searched against UniProtKB/Swiss-Prot release 2026_03 (UniProt Consortium 2025) by BLASTP in DIAMOND v2.0.14 (Buchfink et al. 2021), and its best hit was accepted when the e-value was ≤ 10⁻⁵, the hit covered ≥ 50% of the query at ≥ 30% identity, and the hit was not described as uncharacterized, hypothetical, predicted or of unknown function. Protein domains and Gene Ontology terms were then assigned with InterProScan v5.78-109.0 (Jones et al. 2014). Proteins without a result from either source were searched with BLASTP (BLAST+ v2.12.0; Camacho et al. 2009) against the NCBI non-redundant protein database, restricted to green plants (Viridiplantae; accessed 9 September 2026), with the same criteria. Each gene took its annotation from the first of these sources that gave one, in the order Swiss-Prot, InterPro and nr (Online Resource 6)<sup>(comment 250)</sup>.

---

## Change 59 of 125 · Delete supmat comment 250 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Only 3 of the 55 candidate genes had a functional description in the MorexV3 annotation. We therefore annotated the 23 genes with significant haplotype associations from the protein sequence of their longest isoform in the MorexV3 high-confidence proteome (Mascher et al. 2021). Each protein was first searched against UniProtKB/Swiss-Prot release 2026_03 (UniProt Consortium 2025) by BLASTP in DIAMOND v2.0.14 (Buchfink et al. 2021), and its best hit was accepted when the e-value was ≤ 10⁻⁵, the hit covered ≥ 50% of the query at ≥ 30% identity, and the hit was not described as uncharacterized, hypothetical, predicted or of unknown function. Protein domains and Gene Ontology terms were then assigned with InterProScan v5.78-109.0 (Jones et al. 2014). Proteins without a result from either source were searched with BLASTP (BLAST+ v2.12.0; Camacho et al. 2009) against the NCBI non-redundant protein database, restricted to green plants (Viridiplantae; accessed 9 September 2026), with the same criteria. Each gene took its annotation from the first of these sources that gave one, in the order Swiss-Prot, InterPro and nr (Online Resource 6)<del style="color:#c0392b;background:#fdecea"><sup>(comment 250)</sup></del>.

**Reads after approval:**

> Only 3 of the 55 candidate genes had a functional description in the MorexV3 annotation. We therefore annotated the 23 genes with significant haplotype associations from the protein sequence of their longest isoform in the MorexV3 high-confidence proteome (Mascher et al. 2021). Each protein was first searched against UniProtKB/Swiss-Prot release 2026_03 (UniProt Consortium 2025) by BLASTP in DIAMOND v2.0.14 (Buchfink et al. 2021), and its best hit was accepted when the e-value was ≤ 10⁻⁵, the hit covered ≥ 50% of the query at ≥ 30% identity, and the hit was not described as uncharacterized, hypothetical, predicted or of unknown function. Protein domains and Gene Ontology terms were then assigned with InterProScan v5.78-109.0 (Jones et al. 2014). Proteins without a result from either source were searched with BLASTP (BLAST+ v2.12.0; Camacho et al. 2009) against the NCBI non-redundant protein database, restricted to green plants (Viridiplantae; accessed 9 September 2026), with the same criteria. Each gene took its annotation from the first of these sources that gave one, in the order Swiss-Prot, InterPro and nr (Online Resource 6).

---

## Change 60 of 125 · Delete supmat comment 251

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, we located ten canonical genes in MorexV3: the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*, and the endohydrolases *HvGlb1* and *HvGlb2*. Because their published positions refer to genetic maps or earlier genome assemblies, they were mapped from their reference protein sequences. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 8)<sup>(comment 251)</sup>.

**Reads after approval:**

> To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, we located ten canonical genes in MorexV3: the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*, and the endohydrolases *HvGlb1* and *HvGlb2*. Because their published positions refer to genetic maps or earlier genome assemblies, they were mapped from their reference protein sequences. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene (Online Resource 8)<sup>(comment 251)</sup>.

---

## Change 61 of 125 · Delete supmat comment 251 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, we located ten canonical genes in MorexV3: the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*, and the endohydrolases *HvGlb1* and *HvGlb2*. Because their published positions refer to genetic maps or earlier genome assemblies, they were mapped from their reference protein sequences. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene (Online Resource 8)<del style="color:#c0392b;background:#fdecea"><sup>(comment 251)</sup></del>.

**Reads after approval:**

> To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, we located ten canonical genes in MorexV3: the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*, and the endohydrolases *HvGlb1* and *HvGlb2*. Because their published positions refer to genetic maps or earlier genome assemblies, they were mapped from their reference protein sequences. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene (Online Resource 8).

---

## Change 62 of 125 · Delete supmat comment 252

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 9)<sup>(comment 252)</sup>. As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. 4).

**Reads after approval:**

> To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes (Online Resource 9)<sup>(comment 252)</sup>. As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. 4).

---

## Change 63 of 125 · Delete supmat comment 252 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes (Online Resource 9)<del style="color:#c0392b;background:#fdecea"><sup>(comment 252)</sup></del>. As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. 4).

**Reads after approval:**

> To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes (Online Resource 9). As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. 4).

---

## Change 64 of 125 · Delete supmat comment 253

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 10)<sup>(comment 253)</sup>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource 10)<sup>(comment 254)</sup>.

**Reads after approval:**

> To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource 10)<sup>(comment 253)</sup>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource 10)<sup>(comment 254)</sup>.

---

## Change 65 of 125 · Delete supmat comment 253 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource 10)<del style="color:#c0392b;background:#fdecea"><sup>(comment 253)</sup></del>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource 10)<sup>(comment 254)</sup>.

**Reads after approval:**

> To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource 10). Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource 10)<sup>(comment 254)</sup>.

---

## Change 66 of 125 · Delete supmat comment 254

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource 10)<sup>(comment 253)</sup>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records <del style="color:#c0392b;background:#fdecea"></del>(Online Resource 10)<sup>(comment 254)</sup>.

**Reads after approval:**

> To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource 10)<sup>(comment 253)</sup>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource 10)<sup>(comment 254)</sup>.

---

## Change 67 of 125 · Delete supmat comment 254 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource 10)<sup>(comment 253)</sup>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource 10)<del style="color:#c0392b;background:#fdecea"><sup>(comment 254)</sup></del>.

**Reads after approval:**

> To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024) in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource 10)<sup>(comment 253)</sup>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource 10).

---

## Change 68 of 125 · Delete supmat comment 255

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (Figs. 3, 4) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure <del style="color:#c0392b;background:#fdecea"></del>(Online Resources 11, 12)<sup>(comment 255)</sup>. To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

**Reads after approval:**

> To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (Figs. 3, 4) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure (Online Resources 11, 12)<sup>(comment 255)</sup>. To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

---

## Change 69 of 125 · Delete supmat comment 255 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (Figs. 3, 4) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure (Online Resources 11, 12)<del style="color:#c0392b;background:#fdecea"><sup>(comment 255)</sup></del>. To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

**Reads after approval:**

> To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (Figs. 3, 4) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure (Online Resources 11, 12). To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

---

## Change 70 of 125 · Delete supmat comment 256

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea"></del>The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup><sup>(comment 256)</sup> material-sharing statement<sup>(TODO 235)</sup>

**Reads after approval:**

> The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup><sup>(comment 256)</sup> material-sharing statement<sup>(TODO 235)</sup>

---

## Change 71 of 125 · Delete supmat comment 256 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup><del style="color:#c0392b;background:#fdecea"><sup>(comment 256)</sup></del> material-sharing statement<sup>(TODO 235)</sup>

**Reads after approval:**

> The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup> material-sharing statement<sup>(TODO 235)</sup>

---

## Change 72 of 125 · Delete supmat comment 326

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

The four nutritional traits were quantified in 290 *H. spontaneum* accessions across three growing seasons, and for each trait genotype BLUPs were used to estimate accession-level genotypic values while accounting for season and block effects. Partitioning the phenotypic variance into genetic, season-and-block, genotype-by-season (G×E), and residual components (Fig. 1a) showed that starch had the highest broad-sense heritability (H² = 0.474; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 3<sup>(comment 326)</sup>), roughly twice that of the other nutritional traits. The genotype-by-season interaction was small to moderate for all four traits (4.5–15.2%), and β-glucan combined the largest G×E component (15.2%) with the smallest genetic fraction, indicating a strong genotype-specific response across the three growing seasons. Among morphological traits, flowering time had the highest heritability (H² = 0.747) and tiller number the lowest (H² \< 0.001), indicating the respective genetic control for these traits. To further explore these trends, reaction norms were generated per genotype across seasons (Fig. 1b). They were consistent with the heritability estimates and showed substantial variation among accessions in their response to season, indicating a notable G×E component for all four traits.

**Reads after approval:**

> The four nutritional traits were quantified in 290 *H. spontaneum* accessions across three growing seasons, and for each trait genotype BLUPs were used to estimate accession-level genotypic values while accounting for season and block effects. Partitioning the phenotypic variance into genetic, season-and-block, genotype-by-season (G×E), and residual components (Fig. 1a) showed that starch had the highest broad-sense heritability (H² = 0.474; Online Resource 3<sup>(comment 326)</sup>), roughly twice that of the other nutritional traits. The genotype-by-season interaction was small to moderate for all four traits (4.5–15.2%), and β-glucan combined the largest G×E component (15.2%) with the smallest genetic fraction, indicating a strong genotype-specific response across the three growing seasons. Among morphological traits, flowering time had the highest heritability (H² = 0.747) and tiller number the lowest (H² \< 0.001), indicating the respective genetic control for these traits. To further explore these trends, reaction norms were generated per genotype across seasons (Fig. 1b). They were consistent with the heritability estimates and showed substantial variation among accessions in their response to season, indicating a notable G×E component for all four traits.

---

## Change 73 of 125 · Delete supmat comment 326 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The four nutritional traits were quantified in 290 *H. spontaneum* accessions across three growing seasons, and for each trait genotype BLUPs were used to estimate accession-level genotypic values while accounting for season and block effects. Partitioning the phenotypic variance into genetic, season-and-block, genotype-by-season (G×E), and residual components (Fig. 1a) showed that starch had the highest broad-sense heritability (H² = 0.474; Online Resource 3<del style="color:#c0392b;background:#fdecea"><sup>(comment 326)</sup></del>), roughly twice that of the other nutritional traits. The genotype-by-season interaction was small to moderate for all four traits (4.5–15.2%), and β-glucan combined the largest G×E component (15.2%) with the smallest genetic fraction, indicating a strong genotype-specific response across the three growing seasons. Among morphological traits, flowering time had the highest heritability (H² = 0.747) and tiller number the lowest (H² \< 0.001), indicating the respective genetic control for these traits. To further explore these trends, reaction norms were generated per genotype across seasons (Fig. 1b). They were consistent with the heritability estimates and showed substantial variation among accessions in their response to season, indicating a notable G×E component for all four traits.

**Reads after approval:**

> The four nutritional traits were quantified in 290 *H. spontaneum* accessions across three growing seasons, and for each trait genotype BLUPs were used to estimate accession-level genotypic values while accounting for season and block effects. Partitioning the phenotypic variance into genetic, season-and-block, genotype-by-season (G×E), and residual components (Fig. 1a) showed that starch had the highest broad-sense heritability (H² = 0.474; Online Resource 3), roughly twice that of the other nutritional traits. The genotype-by-season interaction was small to moderate for all four traits (4.5–15.2%), and β-glucan combined the largest G×E component (15.2%) with the smallest genetic fraction, indicating a strong genotype-specific response across the three growing seasons. Among morphological traits, flowering time had the highest heritability (H² = 0.747) and tiller number the lowest (H² \< 0.001), indicating the respective genetic control for these traits. To further explore these trends, reaction norms were generated per genotype across seasons (Fig. 1b). They were consistent with the heritability estimates and showed substantial variation among accessions in their response to season, indicating a notable G×E component for all four traits.

---

## Change 74 of 125 · Delete supmat comment 327

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To examine whether variation in grain nutritional traits followed the environmental gradient of the collection sites, site-mean nutritional-trait BLUPs were correlated with eight environmental variables (Fig. 1d). Starch was positively correlated with clay, precipitation, and organic carbon, and negatively correlated with sand and pH. Fiber showed the opposite pattern for clay, sand, and precipitation, while β-glucan was positively correlated with sand and electrical conductivity and negatively correlated with clay, precipitation, and organic carbon. Protein was not significantly correlated with any of the measured environmental variables. Plotting the accession BLUPs across the 29 sites ordered by ecological region showed a matching descriptive pattern (<del style="color:#c0392b;background:#fdecea"></del>Online Resource 13<sup>(comment 327)</sup>): northern-derived genotypes generally had higher starch, whereas desert-derived genotypes tended to have lower starch and higher fiber and β-glucan, and protein showed no clear directional trend. Together, these associations linked the contrasting patterns of starch versus fiber and β-glucan to gradients in soil texture and precipitation, consistent with the broader geographic pattern observed across the collection sites.

**Reads after approval:**

> To examine whether variation in grain nutritional traits followed the environmental gradient of the collection sites, site-mean nutritional-trait BLUPs were correlated with eight environmental variables (Fig. 1d). Starch was positively correlated with clay, precipitation, and organic carbon, and negatively correlated with sand and pH. Fiber showed the opposite pattern for clay, sand, and precipitation, while β-glucan was positively correlated with sand and electrical conductivity and negatively correlated with clay, precipitation, and organic carbon. Protein was not significantly correlated with any of the measured environmental variables. Plotting the accession BLUPs across the 29 sites ordered by ecological region showed a matching descriptive pattern (Online Resource 13<sup>(comment 327)</sup>): northern-derived genotypes generally had higher starch, whereas desert-derived genotypes tended to have lower starch and higher fiber and β-glucan, and protein showed no clear directional trend. Together, these associations linked the contrasting patterns of starch versus fiber and β-glucan to gradients in soil texture and precipitation, consistent with the broader geographic pattern observed across the collection sites.

---

## Change 75 of 125 · Delete supmat comment 327 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To examine whether variation in grain nutritional traits followed the environmental gradient of the collection sites, site-mean nutritional-trait BLUPs were correlated with eight environmental variables (Fig. 1d). Starch was positively correlated with clay, precipitation, and organic carbon, and negatively correlated with sand and pH. Fiber showed the opposite pattern for clay, sand, and precipitation, while β-glucan was positively correlated with sand and electrical conductivity and negatively correlated with clay, precipitation, and organic carbon. Protein was not significantly correlated with any of the measured environmental variables. Plotting the accession BLUPs across the 29 sites ordered by ecological region showed a matching descriptive pattern (Online Resource 13<del style="color:#c0392b;background:#fdecea"><sup>(comment 327)</sup></del>): northern-derived genotypes generally had higher starch, whereas desert-derived genotypes tended to have lower starch and higher fiber and β-glucan, and protein showed no clear directional trend. Together, these associations linked the contrasting patterns of starch versus fiber and β-glucan to gradients in soil texture and precipitation, consistent with the broader geographic pattern observed across the collection sites.

**Reads after approval:**

> To examine whether variation in grain nutritional traits followed the environmental gradient of the collection sites, site-mean nutritional-trait BLUPs were correlated with eight environmental variables (Fig. 1d). Starch was positively correlated with clay, precipitation, and organic carbon, and negatively correlated with sand and pH. Fiber showed the opposite pattern for clay, sand, and precipitation, while β-glucan was positively correlated with sand and electrical conductivity and negatively correlated with clay, precipitation, and organic carbon. Protein was not significantly correlated with any of the measured environmental variables. Plotting the accession BLUPs across the 29 sites ordered by ecological region showed a matching descriptive pattern (Online Resource 13): northern-derived genotypes generally had higher starch, whereas desert-derived genotypes tended to have lower starch and higher fiber and β-glucan, and protein showed no clear directional trend. Together, these associations linked the contrasting patterns of starch versus fiber and β-glucan to gradients in soil texture and precipitation, consistent with the broader geographic pattern observed across the collection sites.

---

## Change 76 of 125 · Delete supmat comment 328

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 5<sup>(comment 328)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource 3<sup>(comment 329)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

**Reads after approval:**

> Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 5<sup>(comment 328)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource 3<sup>(comment 329)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

---

## Change 77 of 125 · Delete supmat comment 328 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 5<del style="color:#c0392b;background:#fdecea"><sup>(comment 328)</sup></del>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource 3<sup>(comment 329)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

**Reads after approval:**

> Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 5). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource 3<sup>(comment 329)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

---

## Change 78 of 125 · Delete supmat comment 329

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 5<sup>(comment 328)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 3<sup>(comment 329)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

**Reads after approval:**

> Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 5<sup>(comment 328)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource 3<sup>(comment 329)</sup>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

---

## Change 79 of 125 · Delete supmat comment 329 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 5<sup>(comment 328)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource 3<del style="color:#c0392b;background:#fdecea"><sup>(comment 329)</sup></del>), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

**Reads after approval:**

> Overall, 52 SNPs exceeded the genome-wide threshold and resolved into 36 loci, each defined as the lead SNP together with the contiguous block of SNPs in linkage disequilibrium with it (Fig. 2; Online Resource 5<sup>(comment 328)</sup>). Fiber was the most mappable trait, with 24 significant SNPs in 18 loci, followed by β-glucan (20 SNPs, 11 loci), starch (6 SNPs, 5 loci) and protein (2 SNPs, 2 loci). Interestingly, the number of loci did not follow heritability: starch, the most heritable nutritional trait (H² = 0.474; Online Resource 3), yielded only five loci, whereas fiber, with about half its heritability (H² = 0.248), yielded more than three times as many.

---

## Change 80 of 125 · Delete supmat comment 330

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), <del style="color:#c0392b;background:#fdecea"></del>which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><sup>(comment 330)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 81 of 125 · Delete supmat comment 330 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup><del style="color:#c0392b;background:#fdecea"><sup>(comment 330)</sup></del>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

**Reads after approval:**

> The strongest association in the study was a β-glucan peak on chromosome 2H (2H:41,996,626; −log~10~(*P*) = 7.60), and the densest β-glucan signal was a cluster of ten significant SNPs on chromosome 4H (4H:34,825,152–35,062,236), which the 50-kb gap rule split into two adjacent loci<sup>(comment 312)</sup><sup>(comment 311)</sup>. The strongest fiber and starch peaks were both located on the distal arm of chromosome 7H (7H:573,606,306 and 7H:573,606,460; −log~10~(*P*) = 7.25 and 6.44, respectively), whereas both protein loci were located on chromosome 3H, the stronger at 3H:106,623,911 (−log~10~(*P*) = 6.43). Fiber loci were concentrated on chromosome 1H, which carried six of the 18 fiber loci, five of them in two clusters (329.7–330.7 Mb and 440.0–440.6 Mb). Locus size varied widely, from a single SNP (four loci) to 1.73 Mb, but most loci were small, with a median span of 93.5 kb. <sup>(comment 313)</sup>

---

## Change 82 of 125 · Delete supmat comment 331

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci <del style="color:#c0392b;background:#fdecea"></del>(median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

**Reads after approval:**

> Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<sup>(comment 331)</sup>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

---

## Change 83 of 125 · Delete supmat comment 331 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5)<del style="color:#c0392b;background:#fdecea"><sup>(comment 331)</sup></del>. Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

**Reads after approval:**

> Most lead SNPs carried rare alleles, with a minor allele frequency (MAF) below 0.10 at 25 of the 36 loci (median MAF = 0.068; Online Resource 5). Only two lead SNPs, both associated with starch, carried a common minor allele (MAF = 0.42 at 3H:530,660,057 and 6H:525,776,080). <sup>(comment 315)</sup>For β-glucan, fiber and starch, the direction of the allelic effects was consistent across loci. At every β-glucan and fiber locus the minor allele was associated with higher trait values (allelic effects of 0.15–0.36% and 0.12–0.25%, respectively), whereas at every starch locus it was associated with lower starch (−1.31 to −2.69%). The two protein loci showed opposite effects (−0.40% at 3H:106,623,911 and 0.68% at 3H:159,799,839). <sup>(comment 314)</sup>

---

## Change 84 of 125 · Delete supmat comment 332

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. <del style="color:#c0392b;background:#fdecea"></del>Fiber and starch loci co-localized in three genomic regions<sup>(comment 316)</sup><sup>(comment 332)</sup>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

**Reads after approval:**

> To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. Fiber and starch loci co-localized in three genomic regions<sup>(comment 316)</sup><sup>(comment 332)</sup>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

---

## Change 85 of 125 · Delete supmat comment 332 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. Fiber and starch loci co-localized in three genomic regions<sup>(comment 316)</sup><del style="color:#c0392b;background:#fdecea"><sup>(comment 332)</sup></del>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

**Reads after approval:**

> To examine whether the negative relationship between starch and fiber has a genomic counterpart, the positions of the loci were compared across traits. Fiber and starch loci co-localized in three genomic regions<sup>(comment 316)</sup>. On chromosome 7H (573.6 Mb), the fiber locus was contained within the starch locus, and the two lead SNPs were only 154 bp apart. On chromosome 1H, four loci fell within 1.02 Mb (329.7–330.7 Mb): a starch locus flanked by two fiber loci, 121 kb and 58 kb away, and a third fiber locus a further 532 kb beyond the first of them. On chromosome 3H (544.0–546.5 Mb), the widest fiber locus was followed 595 kb downstream by a starch locus. In contrast, none of the β-glucan loci overlapped or neighbored a locus of another trait (the closest, on chromosome 3H, lay 6.0 Mb away), although β-glucan was phenotypically correlated with both starch and fiber (Fig. 1c).

---

## Change 86 of 125 · Delete supmat comment 333

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. Of the 55 candidate genes, <del style="color:#c0392b;background:#fdecea"></del>30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6)<sup>(comment 333)</sup>. The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested. <sup>(comment 317)</sup>

**Reads after approval:**

> Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. Of the 55 candidate genes, 30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6)<sup>(comment 333)</sup>. The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested. <sup>(comment 317)</sup>

---

## Change 87 of 125 · Delete supmat comment 333 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. Of the 55 candidate genes, 30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6)<del style="color:#c0392b;background:#fdecea"><sup>(comment 333)</sup></del>. The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested. <sup>(comment 317)</sup>

**Reads after approval:**

> Candidate genes were further evaluated by haplotype analysis to identify phenotypic divergence between accessions carrying contrasting alleles across haplotype groups. Gene-level haplotype analysis tested whether variation across each gene was associated with the corresponding trait (Kruskal–Wallis, Benjamini–Hochberg correction within each trait, *q* ≤ 0.05) and quantified the direction and magnitude of differences among haplotype groups, information not provided by the single-SNP tests. Of the 55 candidate genes, 30 could be tested and 23 were significant, distributed over 12 of the 36 loci (Online Resource 6). The remaining 25 genes carried no SNPs, too few SNPs, or SNPs too weakly linked to define haplotypes. Within traits, all six testable β-glucan genes and all five testable starch genes were significant, against 12 of 19 for fiber, whereas none of the six protein genes could be tested. <sup>(comment 317)</sup>

---

## Change 88 of 125 · Delete supmat comment 334

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (<del style="color:#c0392b;background:#fdecea"></del>Online Resource 6<sup>(comment 318)</sup><sup>(comment 334)</sup>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

**Reads after approval:**

> A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (Online Resource 6<sup>(comment 318)</sup><sup>(comment 334)</sup>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

---

## Change 89 of 125 · Delete supmat comment 334 (end mark)

*Why:* Its end mark.

**With the changes marked:**

A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (Online Resource 6<sup>(comment 318)</sup><del style="color:#c0392b;background:#fdecea"><sup>(comment 334)</sup></del>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

**Reads after approval:**

> A single locus, the 1.73 Mb fiber locus on chromosome 3H, contained seven of the 12 significant fiber genes. The 23 genes were functionally annotated, and a biological-relevance screen kept only those with a documented or plausible role in the biosynthesis, transport or regulation of the corresponding grain component. The three genes with the clearest such role were two for fiber and one for starch, and none for β-glucan (Online Resource 6<sup>(comment 318)</sup>). Notably, these did not include the two most significant genes: the first (*HORVU.MOREX.r3.6HG0619810*, starch, *q* = 5.5 × 10⁻¹¹) returned no functional call from Swiss-Prot, InterPro or NCBI nr, and the second (*HORVU.MOREX.r3.3HG0301210*, fiber, *q* = 1.1 × 10⁻⁹) carried a transposase-derived domain.

---

## Change 90 of 125 · Delete supmat comment 335

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

The first of the two fiber genes, the glycerol-3-phosphate 2-*O*-acyltransferase *GPAT6* (*HORVU.MOREX.r3.3HG0301300*), resolved into three haplotype groups that differed in grain fiber (*q* = 5.3 × 10⁻⁶, η² = 0.208; <del style="color:#c0392b;background:#fdecea"></del>Online Resources 6, 9<sup>(comment 335)</sup>; Fig. 3a). It belongs to the land-plant GPAT4/6/8 clade, whose members are bifunctional acyltransferase–phosphatases that produce *sn*-2 monoacylglycerol, the committed precursor of cutin (Yang et al. 2010; Petit et al. 2016). Cutin and suberin, such as those of the pericarp and testa, are counted among the constituents of dietary fiber (AACC 2001).

**Reads after approval:**

> The first of the two fiber genes, the glycerol-3-phosphate 2-*O*-acyltransferase *GPAT6* (*HORVU.MOREX.r3.3HG0301300*), resolved into three haplotype groups that differed in grain fiber (*q* = 5.3 × 10⁻⁶, η² = 0.208; Online Resources 6, 9<sup>(comment 335)</sup>; Fig. 3a). It belongs to the land-plant GPAT4/6/8 clade, whose members are bifunctional acyltransferase–phosphatases that produce *sn*-2 monoacylglycerol, the committed precursor of cutin (Yang et al. 2010; Petit et al. 2016). Cutin and suberin, such as those of the pericarp and testa, are counted among the constituents of dietary fiber (AACC 2001).

---

## Change 91 of 125 · Delete supmat comment 335 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The first of the two fiber genes, the glycerol-3-phosphate 2-*O*-acyltransferase *GPAT6* (*HORVU.MOREX.r3.3HG0301300*), resolved into three haplotype groups that differed in grain fiber (*q* = 5.3 × 10⁻⁶, η² = 0.208; Online Resources 6, 9<del style="color:#c0392b;background:#fdecea"><sup>(comment 335)</sup></del>; Fig. 3a). It belongs to the land-plant GPAT4/6/8 clade, whose members are bifunctional acyltransferase–phosphatases that produce *sn*-2 monoacylglycerol, the committed precursor of cutin (Yang et al. 2010; Petit et al. 2016). Cutin and suberin, such as those of the pericarp and testa, are counted among the constituents of dietary fiber (AACC 2001).

**Reads after approval:**

> The first of the two fiber genes, the glycerol-3-phosphate 2-*O*-acyltransferase *GPAT6* (*HORVU.MOREX.r3.3HG0301300*), resolved into three haplotype groups that differed in grain fiber (*q* = 5.3 × 10⁻⁶, η² = 0.208; Online Resources 6, 9; Fig. 3a). It belongs to the land-plant GPAT4/6/8 clade, whose members are bifunctional acyltransferase–phosphatases that produce *sn*-2 monoacylglycerol, the committed precursor of cutin (Yang et al. 2010; Petit et al. 2016). Cutin and suberin, such as those of the pericarp and testa, are counted among the constituents of dietary fiber (AACC 2001).

---

## Change 92 of 125 · Delete supmat comment 336

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

The second, the glucan endo-1,3-β-glucosidase *GH17* (*HORVU.MOREX.r3.5HG0487060*), resolved into five haplotype groups that differed in grain fiber (*q* = 7.5 × 10⁻⁵, η² = 0.222; <del style="color:#c0392b;background:#fdecea"></del>Online Resources 6, 9<sup>(comment 336)</sup>; Fig. 3b). It carries the glycoside hydrolase family 17 catalytic domain together with an X8 carbohydrate-binding module. Barley’s own (1,3;1,4)-β-D-glucan endohydrolases EI and EII belong to this family (Houston et al. 2014), and mixed-linkage β-glucan is an important component of grain dietary fiber (Hrmova et al. 2024). This gene, on 5H, is an uncharacterized member of the family rather than either endohydrolase, which map to 1H and 7H.

**Reads after approval:**

> The second, the glucan endo-1,3-β-glucosidase *GH17* (*HORVU.MOREX.r3.5HG0487060*), resolved into five haplotype groups that differed in grain fiber (*q* = 7.5 × 10⁻⁵, η² = 0.222; Online Resources 6, 9<sup>(comment 336)</sup>; Fig. 3b). It carries the glycoside hydrolase family 17 catalytic domain together with an X8 carbohydrate-binding module. Barley’s own (1,3;1,4)-β-D-glucan endohydrolases EI and EII belong to this family (Houston et al. 2014), and mixed-linkage β-glucan is an important component of grain dietary fiber (Hrmova et al. 2024). This gene, on 5H, is an uncharacterized member of the family rather than either endohydrolase, which map to 1H and 7H.

---

## Change 93 of 125 · Delete supmat comment 336 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The second, the glucan endo-1,3-β-glucosidase *GH17* (*HORVU.MOREX.r3.5HG0487060*), resolved into five haplotype groups that differed in grain fiber (*q* = 7.5 × 10⁻⁵, η² = 0.222; Online Resources 6, 9<del style="color:#c0392b;background:#fdecea"><sup>(comment 336)</sup></del>; Fig. 3b). It carries the glycoside hydrolase family 17 catalytic domain together with an X8 carbohydrate-binding module. Barley’s own (1,3;1,4)-β-D-glucan endohydrolases EI and EII belong to this family (Houston et al. 2014), and mixed-linkage β-glucan is an important component of grain dietary fiber (Hrmova et al. 2024). This gene, on 5H, is an uncharacterized member of the family rather than either endohydrolase, which map to 1H and 7H.

**Reads after approval:**

> The second, the glucan endo-1,3-β-glucosidase *GH17* (*HORVU.MOREX.r3.5HG0487060*), resolved into five haplotype groups that differed in grain fiber (*q* = 7.5 × 10⁻⁵, η² = 0.222; Online Resources 6, 9; Fig. 3b). It carries the glycoside hydrolase family 17 catalytic domain together with an X8 carbohydrate-binding module. Barley’s own (1,3;1,4)-β-D-glucan endohydrolases EI and EII belong to this family (Houston et al. 2014), and mixed-linkage β-glucan is an important component of grain dietary fiber (Hrmova et al. 2024). This gene, on 5H, is an uncharacterized member of the family rather than either endohydrolase, which map to 1H and 7H.

---

## Change 94 of 125 · Delete supmat comment 337

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

For starch, the plastid-localized anion/phosphate transporter *PHT4;3* (*HORVU.MOREX.r3.3HG0301710*), an ortholog of rice PHT4;3, resolved into four haplotype groups that differed in grain starch (*q* = 3.8 × 10⁻⁴, η² = 0.078; <del style="color:#c0392b;background:#fdecea"></del>Online Resources 6, 9<sup>(comment 337)</sup>; Fig. 3c). Phosphate is the allosteric inhibitor of ADP-glucose pyrophosphorylase, which catalyzes the committed step of starch synthesis (Irigoyen et al. 2011), so the rate at which phosphate crosses the plastid envelope can constrain starch accumulation. Consistent with this role, loss of the related plastidic transporter PHT4;2 alters starch accumulation in *Arabidopsis* (Irigoyen et al. 2011).

**Reads after approval:**

> For starch, the plastid-localized anion/phosphate transporter *PHT4;3* (*HORVU.MOREX.r3.3HG0301710*), an ortholog of rice PHT4;3, resolved into four haplotype groups that differed in grain starch (*q* = 3.8 × 10⁻⁴, η² = 0.078; Online Resources 6, 9<sup>(comment 337)</sup>; Fig. 3c). Phosphate is the allosteric inhibitor of ADP-glucose pyrophosphorylase, which catalyzes the committed step of starch synthesis (Irigoyen et al. 2011), so the rate at which phosphate crosses the plastid envelope can constrain starch accumulation. Consistent with this role, loss of the related plastidic transporter PHT4;2 alters starch accumulation in *Arabidopsis* (Irigoyen et al. 2011).

---

## Change 95 of 125 · Delete supmat comment 337 (end mark)

*Why:* Its end mark.

**With the changes marked:**

For starch, the plastid-localized anion/phosphate transporter *PHT4;3* (*HORVU.MOREX.r3.3HG0301710*), an ortholog of rice PHT4;3, resolved into four haplotype groups that differed in grain starch (*q* = 3.8 × 10⁻⁴, η² = 0.078; Online Resources 6, 9<del style="color:#c0392b;background:#fdecea"><sup>(comment 337)</sup></del>; Fig. 3c). Phosphate is the allosteric inhibitor of ADP-glucose pyrophosphorylase, which catalyzes the committed step of starch synthesis (Irigoyen et al. 2011), so the rate at which phosphate crosses the plastid envelope can constrain starch accumulation. Consistent with this role, loss of the related plastidic transporter PHT4;2 alters starch accumulation in *Arabidopsis* (Irigoyen et al. 2011).

**Reads after approval:**

> For starch, the plastid-localized anion/phosphate transporter *PHT4;3* (*HORVU.MOREX.r3.3HG0301710*), an ortholog of rice PHT4;3, resolved into four haplotype groups that differed in grain starch (*q* = 3.8 × 10⁻⁴, η² = 0.078; Online Resources 6, 9; Fig. 3c). Phosphate is the allosteric inhibitor of ADP-glucose pyrophosphorylase, which catalyzes the committed step of starch synthesis (Irigoyen et al. 2011), so the rate at which phosphate crosses the plastid envelope can constrain starch accumulation. Consistent with this role, loss of the related plastidic transporter PHT4;2 alters starch accumulation in *Arabidopsis* (Irigoyen et al. 2011).

---

## Change 96 of 125 · Delete supmat comment 338

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 10<sup>(comment 319)</sup><sup>(comment 338)</sup>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

**Reads after approval:**

> To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; Online Resource 10<sup>(comment 319)</sup><sup>(comment 338)</sup>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

---

## Change 97 of 125 · Delete supmat comment 338 (end mark)

*Why:* Its end mark.

**With the changes marked:**

To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; Online Resource 10<sup>(comment 319)</sup><del style="color:#c0392b;background:#fdecea"><sup>(comment 338)</sup></del>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

**Reads after approval:**

> To ask which of these haplotypes are present in cultivated barley, we compared five European spring malting cultivars released between 2012 and 2018 with the wild haplotype groups at the same three genes (Fig. 3d–f; Online Resource 10<sup>(comment 319)</sup>). At *GPAT6*, the five cultivars matched the two low-fiber haplotype groups at 95–98% and 88–91% of called sites but the high-fiber group at only 39–45%, so the high-fiber wild haplotype is essentially absent from these lines. The cultivars, like the two low-fiber groups, carried the reference (Morex) allele at only 12–14% of these sites, against 65% in the high-fiber group. At *GH17* and *PHT4;3*, by contrast, all five cultivars carried the reference allele at every called site.

---

## Change 98 of 125 · Delete supmat comment 339

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 4; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 9<sup>(comment 339)</sup>). The haplotypes were defined by a single marker group of four SNPs: the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 4a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 4b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 1c). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 4c).

**Reads after approval:**

> In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 4; Online Resource 9<sup>(comment 339)</sup>). The haplotypes were defined by a single marker group of four SNPs: the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 4a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 4b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 1c). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 4c).

---

## Change 99 of 125 · Delete supmat comment 339 (end mark)

*Why:* Its end mark.

**With the changes marked:**

In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 4; Online Resource 9<del style="color:#c0392b;background:#fdecea"><sup>(comment 339)</sup></del>). The haplotypes were defined by a single marker group of four SNPs: the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 4a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 4b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 1c). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 4c).

**Reads after approval:**

> In the haplotype analysis, the gene resolved into two haplotype groups over 246 of the 290 accessions, with the same grouping for both traits (Fig. 4; Online Resource 9). The haplotypes were defined by a single marker group of four SNPs: the three significant SNPs and a fourth in linkage disequilibrium with them (7H:573,606,282; r² = 0.44–0.54). The minority haplotype B comprised exactly the 34 accessions carrying the minor allele at both lead SNPs. It had higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD; Fig. 4a) and lower starch (*P* = 1.1 × 10⁻³, η² = 0.040, 0.56 SD; Fig. 4b) than the major haplotype A (*n* = 212). Thus, a single haplotype at this gene was associated with both traits in opposite directions, matching the negative phenotypic correlation between fiber and starch (Fig. 1c). At the four SNPs that defined the haplotypes, the five elite malting cultivars carried the alleles of the low-fiber, high-starch haplotype A at every called site, and none carried those of the high-fiber, low-starch haplotype B (Fig. 4c).

---

## Change 100 of 125 · Delete supmat comment 423

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Previous studies of Israeli wild barley showed that population structure, flowering time and growth traits are organized along gradients of temperature and precipitation (Hübner et al. 2009, 2013). The present study extends this ecological differentiation to grain nutritional traits. Because the comparisons were based on genotype BLUPs from common-garden experiments, the regional patterns reflect heritable differences among source populations rather than immediate plastic responses to the environments in which the accessions were collected. The contrast between the starch-rich grain of northern genotypes and the fiber- and β-glucan-rich grain of desert genotypes followed soil texture and precipitation at the sites of origin (Fig. 1d; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 13<sup>(comment 423)</sup>), placing the carbohydrate axis on the environmental gradient of the collection.

**Reads after approval:**

> Previous studies of Israeli wild barley showed that population structure, flowering time and growth traits are organized along gradients of temperature and precipitation (Hübner et al. 2009, 2013). The present study extends this ecological differentiation to grain nutritional traits. Because the comparisons were based on genotype BLUPs from common-garden experiments, the regional patterns reflect heritable differences among source populations rather than immediate plastic responses to the environments in which the accessions were collected. The contrast between the starch-rich grain of northern genotypes and the fiber- and β-glucan-rich grain of desert genotypes followed soil texture and precipitation at the sites of origin (Fig. 1d; Online Resource 13<sup>(comment 423)</sup>), placing the carbohydrate axis on the environmental gradient of the collection.

---

## Change 101 of 125 · Delete supmat comment 423 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Previous studies of Israeli wild barley showed that population structure, flowering time and growth traits are organized along gradients of temperature and precipitation (Hübner et al. 2009, 2013). The present study extends this ecological differentiation to grain nutritional traits. Because the comparisons were based on genotype BLUPs from common-garden experiments, the regional patterns reflect heritable differences among source populations rather than immediate plastic responses to the environments in which the accessions were collected. The contrast between the starch-rich grain of northern genotypes and the fiber- and β-glucan-rich grain of desert genotypes followed soil texture and precipitation at the sites of origin (Fig. 1d; Online Resource 13<del style="color:#c0392b;background:#fdecea"><sup>(comment 423)</sup></del>), placing the carbohydrate axis on the environmental gradient of the collection.

**Reads after approval:**

> Previous studies of Israeli wild barley showed that population structure, flowering time and growth traits are organized along gradients of temperature and precipitation (Hübner et al. 2009, 2013). The present study extends this ecological differentiation to grain nutritional traits. Because the comparisons were based on genotype BLUPs from common-garden experiments, the regional patterns reflect heritable differences among source populations rather than immediate plastic responses to the environments in which the accessions were collected. The contrast between the starch-rich grain of northern genotypes and the fiber- and β-glucan-rich grain of desert genotypes followed soil texture and precipitation at the sites of origin (Fig. 1d; Online Resource 13), placing the carbohydrate axis on the environmental gradient of the collection.

---

## Change 102 of 125 · Delete supmat comment 424

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; <del style="color:#c0392b;background:#fdecea"></del>Online Resource 13<sup>(comment 424)</sup>), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

**Reads after approval:**

> Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; Online Resource 13<sup>(comment 424)</sup>), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

---

## Change 103 of 125 · Delete supmat comment 424 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; Online Resource 13<del style="color:#c0392b;background:#fdecea"><sup>(comment 424)</sup></del>), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

**Reads after approval:**

> Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; Online Resource 13), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

---

## Change 104 of 125 · Delete supmat comment 425

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

Most associations rested on rare alleles, and at every fiber, β-glucan and starch locus the minor allele shifted grain composition in the direction of the phenotypic correlations, toward higher fiber and β-glucan and lower starch. This consistency may partly reflect where the carriers came from: at most fiber and β-glucan loci, they originated mainly from desert and coast–desert sites, the populations that already produce low-starch, fiber- and β-glucan-rich grain (<del style="color:#c0392b;background:#fdecea"></del>Online Resources 11, 12<sup>(comment 425)</sup>). Such regional enrichment is expected if these alleles contributed to local adaptation, and does not by itself indicate an artifact. The carriers also came from several desert sites rather than a single one, and at most loci they differed from non-carriers of their own site in the direction of the allelic effect, supporting a genetic effect beyond geographic origin. However, correction for population structure is less effective for rare alleles shared by related accessions, so a contribution of common ancestry cannot be fully excluded. Plausible gene functions further support individual loci. For the candidate genes, however, the allelic effect is harder to separate from geographic origin: the haplotypes associated with high fiber at *GPAT6* and *GH17* came mainly from several desert populations, and the high-starch haplotype of *PHT4;3* was nearly absent from them. The 7H locus and its GDSL esterase/lipase, whose carriers and haplotypes were not concentrated in the desert, therefore provide the clearest case.

**Reads after approval:**

> Most associations rested on rare alleles, and at every fiber, β-glucan and starch locus the minor allele shifted grain composition in the direction of the phenotypic correlations, toward higher fiber and β-glucan and lower starch. This consistency may partly reflect where the carriers came from: at most fiber and β-glucan loci, they originated mainly from desert and coast–desert sites, the populations that already produce low-starch, fiber- and β-glucan-rich grain (Online Resources 11, 12<sup>(comment 425)</sup>). Such regional enrichment is expected if these alleles contributed to local adaptation, and does not by itself indicate an artifact. The carriers also came from several desert sites rather than a single one, and at most loci they differed from non-carriers of their own site in the direction of the allelic effect, supporting a genetic effect beyond geographic origin. However, correction for population structure is less effective for rare alleles shared by related accessions, so a contribution of common ancestry cannot be fully excluded. Plausible gene functions further support individual loci. For the candidate genes, however, the allelic effect is harder to separate from geographic origin: the haplotypes associated with high fiber at *GPAT6* and *GH17* came mainly from several desert populations, and the high-starch haplotype of *PHT4;3* was nearly absent from them. The 7H locus and its GDSL esterase/lipase, whose carriers and haplotypes were not concentrated in the desert, therefore provide the clearest case.

---

## Change 105 of 125 · Delete supmat comment 425 (end mark)

*Why:* Its end mark.

**With the changes marked:**

Most associations rested on rare alleles, and at every fiber, β-glucan and starch locus the minor allele shifted grain composition in the direction of the phenotypic correlations, toward higher fiber and β-glucan and lower starch. This consistency may partly reflect where the carriers came from: at most fiber and β-glucan loci, they originated mainly from desert and coast–desert sites, the populations that already produce low-starch, fiber- and β-glucan-rich grain (Online Resources 11, 12<del style="color:#c0392b;background:#fdecea"><sup>(comment 425)</sup></del>). Such regional enrichment is expected if these alleles contributed to local adaptation, and does not by itself indicate an artifact. The carriers also came from several desert sites rather than a single one, and at most loci they differed from non-carriers of their own site in the direction of the allelic effect, supporting a genetic effect beyond geographic origin. However, correction for population structure is less effective for rare alleles shared by related accessions, so a contribution of common ancestry cannot be fully excluded. Plausible gene functions further support individual loci. For the candidate genes, however, the allelic effect is harder to separate from geographic origin: the haplotypes associated with high fiber at *GPAT6* and *GH17* came mainly from several desert populations, and the high-starch haplotype of *PHT4;3* was nearly absent from them. The 7H locus and its GDSL esterase/lipase, whose carriers and haplotypes were not concentrated in the desert, therefore provide the clearest case.

**Reads after approval:**

> Most associations rested on rare alleles, and at every fiber, β-glucan and starch locus the minor allele shifted grain composition in the direction of the phenotypic correlations, toward higher fiber and β-glucan and lower starch. This consistency may partly reflect where the carriers came from: at most fiber and β-glucan loci, they originated mainly from desert and coast–desert sites, the populations that already produce low-starch, fiber- and β-glucan-rich grain (Online Resources 11, 12). Such regional enrichment is expected if these alleles contributed to local adaptation, and does not by itself indicate an artifact. The carriers also came from several desert sites rather than a single one, and at most loci they differed from non-carriers of their own site in the direction of the allelic effect, supporting a genetic effect beyond geographic origin. However, correction for population structure is less effective for rare alleles shared by related accessions, so a contribution of common ancestry cannot be fully excluded. Plausible gene functions further support individual loci. For the candidate genes, however, the allelic effect is harder to separate from geographic origin: the haplotypes associated with high fiber at *GPAT6* and *GH17* came mainly from several desert populations, and the high-starch haplotype of *PHT4;3* was nearly absent from them. The 7H locus and its GDSL esterase/lipase, whose carriers and haplotypes were not concentrated in the desert, therefore provide the clearest case.

---

## Change 106 of 125 · Delete supmat comment 426

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, <del style="color:#c0392b;background:#fdecea"></del>lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<sup>(comment 414)</sup><sup>(comment 426)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

**Reads after approval:**

> β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<sup>(comment 414)</sup><sup>(comment 426)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

---

## Change 107 of 125 · Delete supmat comment 426 (end mark)

*Why:* Its end mark.

**With the changes marked:**

β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<sup>(comment 414)</sup><del style="color:#c0392b;background:#fdecea"><sup>(comment 426)</sup></del>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

**Reads after approval:**

> β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8)<sup>(comment 414)</sup>, and no SNP on chromosome 7H, which carries *HvCslF6* and *HvGlb2*, was significantly associated with β-glucan (Fig. 2a). The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Houston et al. 2014), is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content. The authors suggested that this conservation reflects the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional: other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021). A panel of European spring barley did detect its region (Shu and Rasmussen 2014), indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

---

## Change 108 of 125 · Delete supmat comment 427

*Why:* Resolved (✅); all supmat comments deleted (user, 2026-10-01).

**With the changes marked:**

The five elite cultivars carried nearly identical genotypes at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every called site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 3f), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 3e). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (<del style="color:#c0392b;background:#fdecea"></del>Online Resources 11, 12<sup>(comment 427)</sup>). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

**Reads after approval:**

> The five elite cultivars carried nearly identical genotypes at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every called site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 3f), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 3e). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resources 11, 12<sup>(comment 427)</sup>). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

---

## Change 109 of 125 · Delete supmat comment 427 (end mark)

*Why:* Its end mark.

**With the changes marked:**

The five elite cultivars carried nearly identical genotypes at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every called site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 3f), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 3e). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resources 11, 12<del style="color:#c0392b;background:#fdecea"><sup>(comment 427)</sup></del>). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

**Reads after approval:**

> The five elite cultivars carried nearly identical genotypes at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every called site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 3f), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 3e). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resources 11, 12). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

---

## Change 110 of 125 · Rewrite TODO 201 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil (Online Resource 1)<sup>(comment 239)</sup>. At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through three rounds of single-seed descent and selfing<sup>(TODO 201)</sup>.

**Reads after approval:**

> We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil (Online Resource 1)<sup>(comment 239)</sup>. At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through three rounds of single-seed descent and selfing<sup>(TODO 201)</sup>.

---

## Change 111 of 125 · Rewrite TODO 202 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026b). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed as described by Potapenko et al. (2026b)<sup>(TODO 202)</sup>. Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow (McKenna et al. 2010; Van der Auwera and O'Connor 2020)<sup>(TODO 203)</sup>. For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, with at most two alternate alleles per site. To reduce the computational load of the large genome, calling was run separately in non-overlapping 60-Mb intervals.

**Reads after approval:**

> Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026b). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed as described by Potapenko et al. (2026b)<sup>(TODO 202)</sup>. Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow (McKenna et al. 2010; Van der Auwera and O'Connor 2020)<sup>(TODO 203)</sup>. For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, with at most two alternate alleles per site. To reduce the computational load of the large genome, calling was run separately in non-overlapping 60-Mb intervals.

---

## Change 112 of 125 · Rewrite TODO 204 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or more than 30% missing data were removed<sup>(TODO 205)</sup>. Applied to all 300 accessions, these filters retained 7,110,996 SNPs. The ten accessions of site 04 (Hachola) were then removed because of an error in the sampling site<sup>(TODO 204)</sup>; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions.

**Reads after approval:**

> Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or more than 30% missing data were removed<sup>(TODO 205)</sup>. Applied to all 300 accessions, these filters retained 7,110,996 SNPs. The ten accessions of site 04 (Hachola) were then removed because of an error in the sampling site<sup>(TODO 204)</sup>; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions.

---

## Change 113 of 125 · Rewrite TODO 205 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or more than 30% missing data were removed<sup>(TODO 205)</sup>. Applied to all 300 accessions, these filters retained 7,110,996 SNPs. The ten accessions of site 04 (Hachola) were then removed because of an error in the sampling site<sup>(TODO 204)</sup>; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions.

**Reads after approval:**

> Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or more than 30% missing data were removed<sup>(TODO 205)</sup>. Applied to all 300 accessions, these filters retained 7,110,996 SNPs. The ten accessions of site 04 (Hachola) were then removed because of an error in the sampling site<sup>(TODO 204)</sup>; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions.

---

## Change 114 of 125 · Rewrite TODO 206 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

The accessions were grown in a common-garden experiment over three growing seasons (2019–20, 2020–21 and 2021–22) in an insect-proof net-house in the Hula Valley, Israel (33°09′08.4″ N 35°37′15.8″ E), in 5-L pots under full irrigation, as described by Potapenko et al. (2026b)<sup>(TODO 206)</sup>. Each season was laid out as a randomized complete block design with one plant per accession in each block: four blocks in 2019–20 and 2020–21, and five in 2021–22.

**Reads after approval:**

> The accessions were grown in a common-garden experiment over three growing seasons (2019–20, 2020–21 and 2021–22) in an insect-proof net-house in the Hula Valley, Israel (33°09′08.4″ N 35°37′15.8″ E), in 5-L pots under full irrigation, as described by Potapenko et al. (2026b)<sup>(TODO 206)</sup>. Each season was laid out as a randomized complete block design with one plant per accession in each block: four blocks in 2019–20 and 2020–21, and five in 2021–22.

---

## Change 115 of 125 · Rewrite TODO 207 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 116 of 125 · Rewrite TODO 208 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 117 of 125 · Rewrite TODO 209 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 118 of 125 · Rewrite TODO 210 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 119 of 125 · Rewrite TODO 211 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

At the end of each season, five spikes per plant were bagged and oven-dried at 35 °C for 30 days. The dried spikes were threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22<sup>(TODO 211)</sup>, and the grain was cleaned with a South Dakota seed blower (Seedburo Equipment Company) and stored for nutritional analysis.

**Reads after approval:**

> At the end of each season, five spikes per plant were bagged and oven-dried at 35 °C for 30 days. The dried spikes were threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22<sup>(TODO 211)</sup>, and the grain was cleaned with a South Dakota seed blower (Seedburo Equipment Company) and stored for nutritional analysis.

---

## Change 120 of 125 · Rewrite TODO 213 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites. Their grain was scanned five times and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), with the protocols adapted to small samples; each accession was assayed in four subsamples of about 200 mg, and their mean was used. Total protein was determined by the Kjeldahl method (nitrogen × 6.25) in three subsamples of about 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College. The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.<sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (Online Resource 2<sup>(comment 240)</sup>).

**Reads after approval:**

> The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites. Their grain was scanned five times and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), with the protocols adapted to small samples; each accession was assayed in four subsamples of about 200 mg, and their mean was used. Total protein was determined by the Kjeldahl method (nitrogen × 6.25) in three subsamples of about 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College. The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.<sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (Online Resource 2<sup>(comment 240)</sup>).

---

## Change 121 of 125 · Rewrite TODO 234 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup><sup>(comment 256)</sup> material-sharing statement<sup>(TODO 235)</sup>

**Reads after approval:**

> The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup><sup>(comment 256)</sup> material-sharing statement<sup>(TODO 235)</sup>

---

## Change 122 of 125 · Rewrite TODO 235 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup><sup>(comment 256)</sup> material-sharing statement<sup>(TODO 235)</sup>

**Reads after approval:**

> The raw sequencing reads of the 300 accessions are available from the European Nucleotide Archive under project PRJEB79623 (Potapenko et al. 2026b); the sample alias of each accession is listed in Online Resource 1.<sup>(TODO 234)</sup><sup>(comment 256)</sup> material-sharing statement<sup>(TODO 235)</sup>

---

## Change 123 of 125 · Rewrite TODO 236 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Generative artificial intelligence tools, Claude Opus and Claude Sonnet (Anthropic; used through Claude Code) and ChatGPT with GPT-5 and GPT-5.5 (OpenAI), were used between February and October 2026 to (i) discuss analysis strategies and the interpretation of results; (ii) write, debug and document the R, Python and shell scripts of the analysis; (iii) check the consistency of the manuscript with the code and its outputs; and (iv) draft and revise the manuscript text. All analytical decisions and interpretations were made by the authors, all code and its outputs were checked by the authors, and all text was reviewed and edited by the authors, who take full responsibility for the content of this publication.<sup>(TODO 236)</sup>

**Reads after approval:**

> Generative artificial intelligence tools, Claude Opus and Claude Sonnet (Anthropic; used through Claude Code) and ChatGPT with GPT-5 and GPT-5.5 (OpenAI), were used between February and October 2026 to (i) discuss analysis strategies and the interpretation of results; (ii) write, debug and document the R, Python and shell scripts of the analysis; (iii) check the consistency of the manuscript with the code and its outputs; and (iv) draft and revise the manuscript text. All analytical decisions and interpretations were made by the authors, all code and its outputs were checked by the authors, and all text was reviewed and edited by the authors, who take full responsibility for the content of this publication.<sup>(TODO 236)</sup>

---

## Change 124 of 125 · Rewrite TODO 401 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; Online Resource 13<sup>(comment 424)</sup>), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

**Reads after approval:**

> Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; Online Resource 13<sup>(comment 424)</sup>), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

---

## Change 125 of 125 · Rewrite TODO 402 (comment text only)

*Why:* Make the comment self-explanatory for other readers (user, 2026-10-01).

**With the changes marked:**

In contrast, protein was not correlated with any of the measured environmental variables, including total soil nitrogen and organic carbon at the source sites, and only weakly and negatively with grain weight (Fig. 1c, d). This pattern suggests that grain protein is governed by processes other than those organizing carbohydrate composition. A negative relationship between grain yield or grain weight and grain protein concentration is well documented (Simmonds 1995) and can be explained by a dilution effect: carbohydrates accumulate faster than nitrogen during grain filling, so protein becomes a smaller proportion of total grain mass (Bogard et al. 2010). Differences between wild and cultivated barley in dry-matter accumulation, nitrogen allocation and nitrogen harvest index further suggested that grain protein depends on the balance between carbon deposition and nitrogen acquisition rather than on direct competition with starch alone (Corke et al. 1989). This separation has a practical implication: the strong starch–fiber trade-off may not constrain protein, so increasing fiber or β-glucan need not cause a reduction in protein. At the genomic level, none of the protein candidate genes could be tested, so the genetic control of protein remains unresolved in this panel<sup>(TODO 402)</sup>.

**Reads after approval:**

> In contrast, protein was not correlated with any of the measured environmental variables, including total soil nitrogen and organic carbon at the source sites, and only weakly and negatively with grain weight (Fig. 1c, d). This pattern suggests that grain protein is governed by processes other than those organizing carbohydrate composition. A negative relationship between grain yield or grain weight and grain protein concentration is well documented (Simmonds 1995) and can be explained by a dilution effect: carbohydrates accumulate faster than nitrogen during grain filling, so protein becomes a smaller proportion of total grain mass (Bogard et al. 2010). Differences between wild and cultivated barley in dry-matter accumulation, nitrogen allocation and nitrogen harvest index further suggested that grain protein depends on the balance between carbon deposition and nitrogen acquisition rather than on direct competition with starch alone (Corke et al. 1989). This separation has a practical implication: the strong starch–fiber trade-off may not constrain protein, so increasing fiber or β-glucan need not cause a reduction in protein. At the genomic level, none of the protein candidate genes could be tested, so the genetic control of protein remains unresolved in this panel<sup>(TODO 402)</sup>.

---

**Your decision:** approve · discard · or tell me what to change.
