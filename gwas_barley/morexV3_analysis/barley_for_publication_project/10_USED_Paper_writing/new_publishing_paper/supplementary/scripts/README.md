# supplementary/scripts — builders of the Online Resources (ESM) of the TAG manuscript

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
| `build_ESM_xlsx.py` | `../ESM_{1,2,3,5,6,7,8,10,11}.xlsx` | the spreadsheet Online Resources |
| `build_ESM_pdf.py` | `../ESM_{4,9,12,13}.pdf` (+ page sources in `tex/`) | the figure Online Resources: title page + figure pages (pdflatex, pdfunite, Ghostscript for font embedding) |
| `build_ESM_extra.py` | `../unreferenced/ESM_extra_phenotype_checks.{pdf,xlsx}` | the **unnumbered supplement, not cited in the text** (user, 2026-10-01): Figs. E1–E7 and their tables from `12_USED_adittional_checks_for_supplementary/02_raw_phenotype_checks/` (E1–E4, E6, E7) and `01_block4_…/` (E5); same title block |
| `../sources/NIR_calibration/` | — | N. Pintel's calibration files, copied unchanged from the user's `NIR.zip` (2026-10-01); read by ESM 2 |
| `fetch_ENA_report.sh` | `../sources/ENA_PRJEB79623_read_run_<date>.tsv` | public ENA run report of PRJEB79623 (read-only query), used by ESM 1 |
| `package_for_submission.py` | `../for_submission/` (`xlsx/`, `pdf/`, `ESM_all_Online_Resources.zip`, `CAPTIONS.md`) | copies of the 13 files for download / upload, one file per Online Resource (TAG); captions taken from the manuscript; refreshed on every build; the unnumbered supplement is not included |
| `build_all.sh` | all of the above | one command |

The files in `../tables/` are intermediate (checked tables); the files to submit are `../ESM_*.xlsx|pdf`, copied with their captions to `../for_submission/`.

## The Online Resources

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

**Layout.** Every `.xlsx` opens with the sheet *About*: title block, "Online Resource N", the full caption, the list of
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
