# supplementary/scripts — builders of the supplementary material (ESM) of the TAG manuscript

**Since 2026-10-17: two files** (S. Hübner's comment; user decision; manuscript item SUPP2). `../ESM_1.xlsx` = Supplementary Tables S1–S22 and `../ESM_2.pdf` = Supplementary Figures S1–S5, one number per table or figure, by first citation. Pipeline (`build_all.sh`): `build_ESM_xlsx.py` builds the content as before, grouped by former Online Resource, into intermediate workbooks `../parts/part_N.xlsx` (never published) → `build_ESM_tables.py` assembles ESM_1.xlsx (Contents sheet + one sheet per table: caption, description, column definitions, table) and checks the manuscript (Tables S1–S22 and Figs. S1–S5 cited in first-citation order, no "Online Resource" left) → `build_ESM_pdf.py` builds ESM_2.pdf (title page with the list, then each figure with its legend on its first page) → `build_ESM_extra.py` (unnumbered supplement) → `package_for_submission.py` (`../for_submission/`) → `package_for_review.py` (`../for_review_Hubner/`). **The one map of the numbers** is `SUPP_TABLES` / `SUPP_FIGURES` in `esm_common.py` (former part and sheet → Table S; figure key → Fig. S, with the captions). The sections below describe the content and sources **by former Online Resource (2026-10-01)**; they are still valid for the content.

**Numbered and built 2026-10-01** (supplementary session; questions and answers `../../build/proposed/ESM0_QUESTIONS.md`;
hand-over `../../build/SUPPLEMENTARY_DONE_2026-10.md`). Numbering = order of first citation in the manuscript
(checked). **Nothing is re-computed:** every value is copied from a read-only project output; each builder asserts the
numbers the manuscript cites from it and stops if one differs (52 checks in `build_ESM_xlsx.py`, plus those of the two
table scripts). The relevance-screen categories (strong / plausible / unlikely / no annotation) are never written
(user rule, 2026-09-30; `esm_common.check_no_screen_categories` checks every cell).

```bash
bash build_all.sh            # everything, about a minute
python3 build_ESM_xlsx.py 1 5 # only some spreadsheets;  python3 build_ESM_pdf.py 9  # only some PDFs
```

## Files

| script | writes | what |
|---|---|---|
| `esm_common.py` | — | TAG title block (**authors, affiliation, e-mail are placeholders**: set them here and run `build_all.sh`), the *About* sheet, table writer, title and screen-category checks |
| `make_ESM_variance_components.py` | `../tables/ESM_variance_components.{tsv,xlsx}` | checked variance-components table (2026-09-30), input of ESM 3 |
| `make_ESM_candidate_genes.py` | `../tables/ESM_candidate_genes.{tsv,xlsx}` | checked 55-gene table (2026-09-30), input of ESM 5, 6, 7 |
| `build_ESM_xlsx.py` | `../parts/part_{1,2,3,5,6,7,8,10,11}.xlsx` | the content of the tables, grouped by former Online Resource (intermediate since 2026-10-17) |
| `build_ESM_tables.py` | `../ESM_1.xlsx` | **Supplementary Tables S1–S22** (2026-10-17): assembled from the parts with `esm_common.SUPP_TABLES`; also checks the citation order in the manuscript |
| `esm_assemble.py` | — | xlsx helpers (markup, About-sheet parser, table copy, Contents sheet), moved from `package_for_review.py` on 2026-10-17 |
| `build_ESM_pdf.py` | `../ESM_2.pdf` (+ page sources in `tex/`) | **Supplementary Figures S1–S5** (one file since 2026-10-17; `figs()` holds the legends and panels, also read by `package_for_review.py`); its page helpers are reused by `build_ESM_extra.py` |
| `build_ESM_extra.py` | `../unreferenced/ESM_extra_phenotype_checks.{pdf,xlsx}` | the **unnumbered supplement, not cited in the text** (user, 2026-10-01): Figs. E1–E7 and their tables from `12_USED_adittional_checks_for_supplementary/02_raw_phenotype_checks/` (E1–E4, E6, E7) and `01_block4_…/` (E5); same title block |
| `../sources/NIR_calibration/` | — | N. Pintel's calibration files, copied unchanged from the user's `NIR.zip` (2026-10-01); read by ESM 2 |
| `fetch_ENA_report.sh` | `../sources/ENA_PRJEB79623_read_run_<date>.tsv` | public ENA run report of PRJEB79623 (read-only query), used by ESM 1 |
| `package_for_submission.py` | `../for_submission/` (`ESM_1.xlsx`, `ESM_2.pdf`, `ESM_all_supplementary.zip`, `CAPTIONS.md`) | the two files for upload, with the two captions from the manuscript |
| `package_for_review.py` | `../for_review_Hubner/` | S. Hübner's review package: `Supplementary_Tables.xlsx` (= ESM_1.xlsx) and `Supplementary_Figures.docx` (Figs. S1–S5 with legends), the unnumbered supplement as xlsx + docx, and a zip. Run by `build_all.sh` since 2026-10-17 |
| `build_all.sh` | all of the above | one command (about a minute) |

The files in `../tables/` and `../parts/` are intermediate; the files to submit are `../ESM_1.xlsx` and `../ESM_2.pdf`, copied with their captions to `../for_submission/`.

## Content and sources, by former Online Resource (2026-10-01; new numbers: `esm_common.SUPP_TABLES` / `SUPP_FIGURES`, list in `../supplementary.md`)

| # | file | content | sources (read-only, project root) | cited in the manuscript |
|---|---|---|---|---|
| 1 | `ESM_1.xlsx` | sites (30; region, coordinates, 10 environmental variables), accessions (290; ENA alias = accession ID, VCF name, ENA sample and run accessions, site, region, genotype BLUPs of 10 traits), environmental correlation matrix (29 sites) | `00_THIN_…/01_…R` (site names, regions), `all years barley.csv` (coordinates), `merged_Evgeny_env_fitness_environmental_30_sites.csv`, `outputs/subsection_1/tables/BLUP_10_traits.csv`, `C2_*/tables/`, `../sources/ENA_…tsv` | M&M §1 ¶1, §4 ¶2 (BLUPs), §4 ¶5, Data and material availability |
| 2 | `ESM_2.xlsx` | NIR calibration: the 30 calibration accessions (laboratory values, NIR predictions with the factory calibration, 2019–20 values) and their correlations (Q17 option A: no β-glucan lab-vs-calibration r, which is 1 by construction; no fiber); inter-block correlations 2019–20 (4 traits × 6 pairs) | `../sources/NIR_calibration/` (N. Pintel's files, README there); `12_USED_adittional_checks_for_supplementary/01_block4_…/results/tables/Table_OnlineResource_2019_block_correlations.tsv` | M&M §3, §4 ¶1 |
| 3 | `ESM_3.xlsx` | variance components and H², 8 traits | `../tables/ESM_variance_components.tsv` ← `00_THIN_…/tables/Variance_components_GxE.csv`, `H2_GxE.csv` | M&M §4 ¶3; Results Ch. 1 ¶1, Ch. 2 ¶1 |
| 4 | `ESM_4.pdf` | LD decay: a genome-wide, b per chromosome (same run: 7,110,996 SNPs, 290 accessions) | `02_USED_LD_decay_V2_wholegenome/results/figure_ld_decay_*.pdf`, `ld_decay_summary_for_results_chapter.txt` | M&M §5 ¶3 |
| 5 | `ESM_5.xlsx` | 36 loci (+ minor-allele carriers / accessions called, Q16); 52 significant SNPs | `01_…/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv`, `Table_loci_members_full.tsv`, `00_snp_level/significant_snps.tsv`; carriers from `11_…/results/tables/T02_lead_enrichment_tests.tsv` | M&M §6 ¶1; Results Ch. 2 ¶1 |
| 6 | `ESM_6.xlsx` | the 55 candidate genes (replaces Table 2) | `../tables/ESM_candidate_genes.tsv` (sources in `make_ESM_candidate_genes.py`) | M&M §6 ¶2, §7 ¶3, §8 ¶1; Results Ch. 3 ¶2–¶5 |
| 7 | `ESM_7.xlsx` | choice of ε: summary + 55 genes × 5 ε | `04_…/loci_LDspan_eps06_V4/Diagnostics/epsilon_choice_{supplementary,per_gene}.tsv` | M&M §7 ¶2 |
| 8 | `ESM_8.xlsx` | ten canonical β-glucan genes on MorexV3 + distances | `05_…/09_USED_canonical_betaglucan_gene_check/results/tables/canonical_betaglucan_gene_distances.tsv` (mapping only), `canonical_gene_distances_current_loci.tsv` | M&M §8 ¶3; Discussion §3 |
| 9 | `ESM_9.pdf` | crosshap plots (page 5 of each combined PDF) of GPAT6, GH17, PHT4;3 (MGmin 2, ε 0.6) and the GDSL gene (fiber and starch pages; MGmin 3, ε 0.9). Replaces the 7H MGmin × ε grid (dropped, user) | `04_…/loci_LDspan_eps06_V4/CombinedPDF/`, `03_01_7H_…/results/figures/crosshap/` | M&M §8 ¶4 |
| 10 | `ESM_10.xlsx` | elite cultivars: 5 lines + call rates (3 genes, pooled, GDSL), 136-line screen, site accounting (4 genes), 757 excluded sites | `07_…/results/tables/Table_elite_lines.tsv`, `Table_elite_line_screen.tsv`, `Table_site_overlap.tsv`, `Table_removed_sites.tsv`; `03_01_7H_…/results/tables/` (same) | M&M §10 ¶1; Results Ch. 3 ¶6 |
| 11 | `ESM_11.xlsx` | carrier origin: T02 (+T01), T04 (+T03), T08; `desert_OR` and the extreme-accession check dropped (user, Q7) | `11_USED_…/results/tables/T01–T04, T08` | M&M §11; Discussion §3, §5 (with 12) |
| 12 | `ESM_12.pdf` | Figs. S1–S2 of step 11 (600 dpi PNG, captions from `results/CAPTIONS.md`) | `11_USED_…/results/figures/Fig_carrier_origin.png`, `Fig_desert_enrichment.png` | with 11 |
| 13 | `ESM_13.pdf` | site figure (former Fig. 2a) | `08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf` | Results Ch. 1 ¶3; Discussion §1, §2 |

**Layout until 2026-10-17 (superseded; current layout at the top of this file).** Every `.xlsx` opened with the sheet *About*: title block, "Online Resource N", the full caption, the list of
sheets and the column definitions; the data sheets start with their header row. Every PDF opens with the title block
and the full legend; each figure page carries "Online Resource N", its panel label and an identification line. Project
paths are not written into the published files (they are here and in the builders' docstrings).

**Decisions applied** (ESM0, 2026-10-01): EC in mS/cm as in the source and Fig. 1d (Q11d); ENA sample/run accessions
kept (Q11a); no unit for plant height and grain weight, as in Potapenko et al. 2026b, Mol Ecol (Q11b); `desert_OR` and the
extreme-accession check left out of ESM 11 (Q7); the 7H grid dropped and both GDSL crosshap pages shown (Q1, Q10); the
starch values of 2021–22 in the unnumbered supplement carry the analysis correction, never mentioned in the files (Q12a).

**Supplement not cited in the text** (user, 2026-10-01): `../unreferenced/ESM_extra_phenotype_checks.{pdf,xlsx}` (raw-phenotype
checks, `build_ESM_extra.py`); unnumbered, so it is outside TAG's citation rule until the user moves a part of it into the
numbered set (it then needs a number and a citation in the text).

## History

- 2026-09-30: `make_ESM_variance_components.py`, `make_ESM_candidate_genes.py` (S. Hübner's structural comments; Tables 1, 2 → Online Resources).
- 2026-10-01: numbering, title blocks, all 13 files (`esm_common.py`, `build_ESM_xlsx.py`, `build_ESM_pdf.py`, `fetch_ENA_report.sh`, `build_all.sh`); the unnumbered supplement (`build_ESM_extra.py`).
- 2026-10-17 (item MRG4; user approval): the manuscript's Figs. 3 and 4 were combined into one Fig. 3 (row 4 = the GDSL gene). Two texts that cited the old Fig. 4 were updated: the legend of Online Resource 9 (`build_ESM_pdf.py`: "those of Fig. 3; their tests are in Fig. 3 and Online Resource 6") and the "Haplotype groups" definition of Online Resource 11 (`build_ESM_xlsx.py`, sheet About). `build_all.sh` re-run: only ESM_9 and ESM_11 changed in content (every other file, the unnumbered supplement and `CAPTIONS.md` identical in content, checked against a copy taken before the run). The review package for S. Hübner (`../for_review_Hubner/`, `package_for_review.py`, first built 2026-10-04) was rebuilt the same day at the user's request; only the same two texts changed.
- 2026-10-17 (item SUPP2; S. Hübner's comment, user decision): renumbered into two running series, Supplementary Tables S1–S22 (`ESM_1.xlsx`) and Supplementary Figures S1–S5 (`ESM_2.pdf`). New: `build_ESM_tables.py`, `esm_assemble.py`, the map in `esm_common.py`; `build_ESM_xlsx.py` writes `../parts/`; `build_ESM_pdf.py` writes one PDF (`figs()`); `package_for_submission.py` and `package_for_review.py` rewritten for the two files; cross-references inside the files name the new numbers (step 11's draft "Fig. S1, S2" mapped to Figs. S3, S4). Checked: the 22 tables equal the former sheets cell by cell; the manuscript cites Tables S1–S22 and Figs. S1–S5 in first-citation order. The former files `ESM_2.xlsx` … `ESM_13.pdf` were deleted.
