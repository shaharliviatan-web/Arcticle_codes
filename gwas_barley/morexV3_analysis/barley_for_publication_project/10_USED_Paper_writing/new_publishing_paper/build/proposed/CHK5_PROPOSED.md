# Review item CHK5: proposed changes

Proposal only. `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item CHK5** — checklist item 5 points to the supplement map and the upload folder (user request: easy to know where each thing is).
>

## Change 1 of 1 · Checklist at the top, item 5

*Why:* Where the files are.

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References:** build the list in TAG style (author–year in the text with no comma; list alphabetical; journal names abbreviated (ISSN LTWA); full DOI links); verify every citation marked "verify" in the TODO comments, and add the β-amylase citation (TODO 411).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). <del style="color:#c0392b;background:#fdecea">List: `supplementary/supplementary.md`.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`.</ins> Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References:** build the list in TAG style (author–year in the text with no comma; list alphabetical; journal names abbreviated (ISSN LTWA); full DOI links); verify every citation marked "verify" in the TODO comments, and add the β-amylase citation (TODO 411).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

**Your decision:** approve · discard · or tell me what to change.
