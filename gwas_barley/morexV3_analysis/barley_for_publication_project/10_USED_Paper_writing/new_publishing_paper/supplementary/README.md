# supplementary/ — the supplementary material (ESM) of the TAG manuscript: where everything is

**Since 2026-10-17 (S. Hübner; user decision): two files**, `ESM_1.xlsx` = **Supplementary Tables S1–S22** (one sheet per table) and `ESM_2.pdf` = **Supplementary Figures S1–S5** (each with its legend), numbered by first citation and cited "Table S5", "Fig. S2d, e". Until then: 13 "Online Resources" `ESM_1` … `ESM_13` (map former → new: `supplementary.md`).

Built 2026-10-01 (supplementary session). TAG rules for the supplement: `../../TAG_requirements.md`, section
"Supplementary Information (ESM)".

## Quick answers

| I want… | go to |
|---|---|
| **to download or upload the files** | **`for_submission/`**: `ESM_1.xlsx`, `ESM_2.pdf`, `ESM_all_supplementary.zip`, `CAPTIONS.md` (each file with its caption, for the submission system) |
| the final list: what each table and figure holds, and the former Online Resource it came from | `supplementary.md` (top table); the map in code: `scripts/esm_common.py` (`SUPP_TABLES`, `SUPP_FIGURES`) |
| how the files are built, their sources and checks | `scripts/README.md` |
| to rebuild everything after a change (e.g. the author names) | `bash scripts/build_all.sh` (about a minute; refreshes `for_submission/` and `for_review_Hubner/`, and checks the citation order in the manuscript) |
| to set the author names, affiliation and e-mail of the title blocks | `scripts/esm_common.py` (`AUTHORS`, `CORR_NAME`, `CORR_AFFILIATION`, `CORR_EMAIL`), then `build_all.sh` |
| the captions in the manuscript | last section of the working file, "Supplementary Information" (one caption per file); each table's and figure's own caption is inside the files |
| the status of each citation in the manuscript | all resolved: the 37 "supmat" Word comments that tracked them were deleted on 2026-10-01 (rule in `../../CLAUDE.md`) |
| the decisions and their history | `../build/proposed/ESM0_QUESTIONS.md` (questions and answers), `../build/SUPPLEMENTARY_DONE_2026-10.md` (hand-over), end of `../build/REVIEW_2026-09-27_Results_Discussion.md` (log) |
| the supplement that is not cited in the text | `unreferenced/` |
| the review package for S. Hübner | `for_review_Hubner/`: `Supplementary_Tables.xlsx` (= ESM_1.xlsx), `Supplementary_Figures.docx` (Figs. S1–S5, editable), the uncited raw-phenotype files, a zip (`scripts/package_for_review.py`, run by `build_all.sh`) |

## Folder map

```
supplementary/
├── README.md                this map
├── supplementary.md         final list (number, content, status) + the planning list kept for the record
├── ESM_1.xlsx              SUPPLEMENTARY TABLES S1–S22 (one sheet per table; scripts/build_ESM_tables.py)
├── ESM_2.pdf               SUPPLEMENTARY FIGURES S1–S5 (each with its legend; scripts/build_ESM_pdf.py)
├── parts/                  intermediate workbooks part_N.xlsx (content grouped by former Online Resource; never published)
├── for_submission/         copies for upload: ESM_1.xlsx, ESM_2.pdf, ESM_all_supplementary.zip, CAPTIONS.md
│                            (rewritten on every build: do not edit by hand)
├── for_review_Hubner/      S. Hübner's review package (his request, 2026-10-04): Supplementary_Tables.xlsx + Supplementary_Figures.docx, the
│                            unnumbered supplement likewise, and a zip; `scripts/package_for_review.py` (run after build_all.sh)
├── unreferenced/            ESM_extra_phenotype_checks.{pdf,xlsx}: raw-phenotype checks, NOT cited, unnumbered
│                            (the user decides later whether a part joins the numbered set)
├── scripts/                 the builders (README.md there): esm_common.py, build_ESM_xlsx.py, build_ESM_pdf.py,
│                            build_ESM_extra.py, package_for_submission.py, fetch_ENA_report.sh, build_all.sh,
│                            make_ESM_variance_components.py, make_ESM_candidate_genes.py; tex/ = PDF page sources
├── sources/                 inputs that are not project outputs:
│   ├── ENA_PRJEB79623_read_run_2026-10-01.tsv   public ENA run report (Table S5)
│   └── NIR_calibration/     N. Pintel's calibration files, copied from the user's NIR.zip (Tables S2, S3; README there)
└── tables/                  intermediate checked tables (variance components, candidate genes) read by the builders
```

## Where the content comes from

Every value is copied from a read-only project output; nothing is re-estimated (except the Pearson P values of the
NIR-calibration sheet, recomputed from N. Pintel's data and identical to R's `cor.test`). The analysis folders behind
each file are listed in `scripts/README.md`. Supporting analyses written for the supplement:
`../../../12_USED_adittional_checks_for_supplementary/` (01 block-4 correlations → Table S4; 02 raw-phenotype checks →
`unreferenced/`). The original NIR files: `../../section_materials/02_materials_methods/NIR_calibration_files/NIR.zip`.

## Rules that apply here (from `../../CLAUDE.md` and the user's decisions)

- Two files since 2026-10-17: `ESM_1.xlsx` (Tables S1–S22) and `ESM_2.pdf` (Figs. S1–S5), each table and figure numbered by first
  citation; each file carries the TAG title block; each table and figure carries its own caption (and column definitions).
- The relevance-screen categories (strong / plausible / unlikely / no annotation) never appear (checked by the builders).
- The +8.85 adjustment of the 2021–22 starch values is never mentioned in any supplementary file (user, 2026-10-01).
- While a supplementary citation is unresolved, it carries a "supmat" Word comment (❌); it is deleted once resolved. None is open.
- Delete before submission: the open-items list at the top of the manuscript.
