# Supplementary material (Online Resources): hand-over, 2026-10-01

Checklist item 5 of the manuscript. **Where everything is: `supplementary/README.md`** (folder map, download folder `supplementary/for_submission/`, rebuild command). Questions and answers are in `build/proposed/ESM0_QUESTIONS.md`; the decision log is
at the end of `build/REVIEW_2026-09-27_Results_Discussion.md`. Nothing was committed.

## 1. Final list

The files are numbered by first citation in the manuscript (checked: 1 → 13 in reading order). They are all in
`new_publishing_paper/supplementary/` and are rebuilt with `supplementary/scripts/build_all.sh` (13 s; 52 checks
against the manuscript's numbers pass). The full source paths are in `supplementary/scripts/README.md`.

| ESM | file | content | main sources | cited in |
|---|---|---|---|---|
| 1 | `ESM_1.xlsx` | 30 sites (region, coordinates, 10 environmental variables); 290 accessions (accession ID = ENA alias, VCF name, ENA sample and run accessions, site, region, genotype BLUPs of 10 traits); environmental correlations over 29 sites | `00_THIN_…` (site maps, `all years barley.csv`, 30-site environment file, `BLUP_10_traits.csv`, `C2_*` tables); ENA report (`supplementary/sources/`) | M&M §1 ¶1, §4 ¶2 (BLUPs, new), §4 ¶5, Data and material availability |
| 2 | `ESM_2.xlsx` | NIR calibration (30 calibration accessions: laboratory values, factory-calibration predictions, 2019–20 values; correlations, thesis Table 1 reproduced, β-glucan r = 1 and fiber left out) + inter-block correlations 2019–20 | `supplementary/sources/NIR_calibration/` (N. Pintel's files); `12_…/01_block4_…/Table_OnlineResource_2019_block_correlations.tsv` | M&M §3, §4 ¶1 |
| 3 | `ESM_3.xlsx` | variance components and H², 8 traits | `00_THIN_…/Variance_components_GxE.csv`, `H2_GxE.csv` | M&M §4 ¶3; Results Ch. 1, Ch. 2 |
| 4 | `ESM_4.pdf` | LD decay: a genome-wide, b per chromosome (same run, legitimate) | `02_USED_LD_decay_V2_wholegenome/results/` | M&M §5 ¶3 |
| 5 | `ESM_5.xlsx` | 36 loci (with minor-allele carrier counts); 52 significant SNPs | `01_…/02_loci_FINAL/tables/`, `00_snp_level/significant_snps.tsv` | M&M §6 ¶1; Results Ch. 2 ¶1, ¶3 |
| 6 | `ESM_6.xlsx` | 55 candidate genes | `03_…`, `04_…/Stats`, `05_…/08_…`, `06_…` | M&M §6–§8; Results Ch. 3 (4×) |
| 7 | `ESM_7.xlsx` | choice of ε (summary; 55 genes × 5 ε) | `04_…/Diagnostics/epsilon_choice_*.tsv` | M&M §7 ¶2 |
| 8 | `ESM_8.xlsx` | ten canonical β-glucan genes on MorexV3, with distances | `05_…/09_USED_canonical_betaglucan_gene_check/results/tables/` | M&M §8 ¶3; Discussion §3 (new) |
| 9 | `ESM_9.pdf` | crosshap plots (page 5) of GPAT6, GH17, PHT4;3 (MGmin 2, ε 0.6) and the GDSL gene (fiber and starch; MGmin 3, ε 0.9); **replaces the 7H grid** (user) | `04_…/CombinedPDF/`, `03_01_…/results/figures/crosshap/` | M&M §8 ¶4 (7H); Results Ch. 3 (3×), Ch. 4 ¶3 |
| 10 | `ESM_10.xlsx` | elite cultivars: 5 lines + call rates, 136-line screen, site accounting (4 genes), 757 excluded sites | `07_…/results/tables/`, `03_01_…/results/tables/` | M&M §10 ¶1; Results Ch. 3 ¶6 |
| 11 | `ESM_11.xlsx` | carrier origin: T02 (+T01), T04 (+T03), T08 (`desert_OR` and the extreme-accession check dropped) | `11_…/results/tables/` | M&M §11; Discussion §3, §5 |
| 12 | `ESM_12.pdf` | carrier origin, Figs. S1–S2 | `11_…/results/figures/` | with 11 |
| 13 | `ESM_13.pdf` | site figure (former Fig. 2a) | `08_USED_creating_figures/Figure_ESM_site_BLUPs/` | Results Ch. 1 ¶3; Discussion §1, §2 |

Every file starts with the TAG title block. In a workbook this is the *About* sheet, which also holds the caption and
the column definitions; in a PDF it is page 1, followed by one page per figure.

**For download and upload:** `supplementary/for_submission/`: `xlsx/` (9 files), `pdf/` (4 files),
`ESM_all_Online_Resources.zip` (all 13), and `CAPTIONS.md` (each file with its caption, for the submission system).
The files stay separate, one per Online Resource, as TAG asks.

The supplement that is **not cited in the text** is `supplementary/unreferenced/ESM_extra_phenotype_checks.{pdf,xlsx}`
(Figs. E1–E7 with their tables). It is unnumbered.

## 2. Files and scripts created or changed

- **Supplement:**
  - New scripts in `supplementary/scripts/`: `esm_common.py`, `build_ESM_xlsx.py`, `build_ESM_pdf.py` (page sources in
    `scripts/tex/`), `build_ESM_extra.py`, `fetch_ENA_report.sh`, `build_all.sh`. `README.md` is rewritten.
  - New outputs: `supplementary/ESM_1`–`13`, `supplementary/unreferenced/`, and
    `supplementary/sources/ENA_PRJEB79623_read_run_2026-10-01.tsv`.
  - The two table scripts from 2026-09-30 are unchanged; their TSVs were re-run and are byte-identical.
- **Check folder:**
  - `12_USED_adittional_checks` was renamed to `12_USED_adittional_checks_for_supplementary`. Check 01's path was
    updated and a re-run reproduced every output (only the run date differs).
  - New check `02_raw_phenotype_checks/` (scripts 00–04, `run_all.sh`, README, figures E1–E4, E6 and E7, tables T1–T5,
    `results_numbers.txt`).
  - Both READMEs were updated.
- **Documentation:**
  - `10_USED_Paper_writing/CLAUDE.md`: the 12_ module, the "supmat" comment rule and the next free IDs, the numbered
    Online Resources, the layout, and a session entry.
  - `supplementary/supplementary.md`: the final numbering, with the planning list kept below it.
  - `build/DECLARATIONS_HANDOVER.md`: Online Resource 1 is now named.

## 3. Manuscript changes (all through `review_tool.py`)

| item | what |
|---|---|
| SUPMAT1 | One "supmat" Word comment (❌/✅) on each of the 36 supplementary places (IDs 239–256, 326–338, 423–427; next free 257, 339, 428); also the renamed path in TODO 215 and its hidden note |
| CHK1 | Checklist item 1: the ESM title blocks need the authors, affiliation and e-mail |
| E1 | 34 edits: every placeholder now has its number; 31 resolved TODOs were removed; TODOs 214 and 234 kept (reworded or updated); your comments 311, 312, 316, 318, 319 and 414 kept; "(Online Resource 1)" added for the BLUPs (Q2) and "(Online Resource 8)" in the Discussion (Q3). The Word build now has 114 → 83 comments |
| CHK2, CHK3 | Checklist item 5 is marked (✅ done / ⬜ open) |
| CHK6 | The open-items list rewritten as plain status lines with no paths (item 4 unchanged); the 37 supmat comments (all ✅), the user's answered comments (311–319, 414–416) and TODO 203 deleted; the 16 other open TODOs rewritten to explain each problem. 23 Word comments remain, all open TODOs |
| E3 | NIR facts from N. Pintel's files (Q17): M&M "scanned five times", "four subsamples" (kits), "Kjeldahl method (nitrogen × 6.25) in three subsamples"; TODO 212 and 213 comments get the facts (both stay open: the citation; S. Hübner to inspect the calibration across seasons); TODO 214 resolved; caption of Online Resource 2 extended; checklist ✅. All 37 supmat comments ✅; 83 Word comments |
| E2 | Citations only (Q10, Q13): Results Ch. 2 ¶3 "(median MAF = 0.068; Online Resource 5)"; Ch. 3 ¶1 "(Online Resource 6)"; GPAT6, GH17, PHT4;3 "Online Resources 6, 9"; Ch. 4 ¶3 "(Fig. 4; Online Resource 9)"; new supmat 339. Now 37 supmat comments, 36 ✅; the only ❌ is 240 (NIR calibration) |
| CAP1, CHK4 | **Applied (Q15):** a "Supplementary Information" section after the Conclusions with 13 concise captions (about 270 words); matches published TAG papers (checked in PMC); enter the same captions in the submission system |

## 4. Noticed outside the scope (not changed)

- **Negative raw starch values in 2021–22.** 25 of 1,163 plants are below 0 % (minimum −21.6 %). Four remain negative
  even after the +8.85 that step 00 adds, and the ±3 SD rule removes only values below about −6 %. Worth raising with
  N. Pintel together with TODO 213 (the calibration across seasons).
- **Two hidden notes still carry old Online Resource numbers.** Results Ch. 2 ¶3 says "Online Resource 1" (now 5), and
  the note after Ch. 3 ¶5 says "Online Resource 2" (now 6). Both are dated records and are invisible in Word.
- **Results Ch. 3 ¶1** says the 25 untestable genes "carried no SNPs, too few SNPs, or SNPs too weakly linked". One of
  them (3HG0301250) failed with a crosshap internal error, which ESM 6 shows (already noted 2026-09-30).
- **Two features of the source figures were kept as is.** The per-chromosome LD figure has a title inside the image.
  The crosshap PDFs reference non-embedded fonts; ESM 9 embeds substitutes.
- **LuaLaTeX cannot load system fonts** on fidel (luaotfload is missing), so the PDFs use pdfLaTeX with Helvetica
  (Arial metrics).

## 5. Still open

- **Title blocks:** author names, and the corresponding author's affiliation and e-mail. They are placeholders in
  `scripts/esm_common.py`; rerun `build_all.sh` once they are set (checklist item 1).
- **Q5 / Q17:** done (N. Pintel's NIR.zip; ESM 2 option A; E3).
- **CAP1:** approve the captions and their position.
- Answered and done since: Q10 (both GDSL pages; citations), Q11 (ENA kept; no units for plant height and grain weight; BLUPs centered), Q12 (starch with the analysis correction, not reported), Q13, Q14, Q15 (captions applied), Q16 (carrier counts in ESM 5).
- **NIR files:** upload to `10_USED_Paper_writing/section_materials/02_materials_methods/NIR_calibration_files/` (README there).
- **TODO 234:** the Zenodo DOI and the GitHub repository (for S. Hübner).
- **TODO 213 (S. Hübner):** whether one NIR calibration served all three seasons (the 2020 and 2021 exports name no calibration).
- **Fiber, for the record:** thesis Table 1's fiber r = 0.526 does not reproduce from the source data (r = 0.25, P = 0.19; 5 of 30 rows of the CSV it came from differ); it is not used anywhere.
