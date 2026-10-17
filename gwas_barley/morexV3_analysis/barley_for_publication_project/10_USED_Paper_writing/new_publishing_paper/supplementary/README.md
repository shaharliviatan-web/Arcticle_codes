# supplementary/ — the Online Resources (ESM) of the TAG manuscript: where everything is

Built 2026-10-01 (supplementary session). TAG rules for the supplement: `../../TAG_requirements.md`, section
"Supplementary Information (ESM)".

## Quick answers

| I want… | go to |
|---|---|
| **to download or upload the files** | **`for_submission/`**: `xlsx/` (9), `pdf/` (4), `ESM_all_Online_Resources.zip` (all 13), `CAPTIONS.md` (each file with its caption, for the submission system) |
| the final list: what each Online Resource holds and where it is cited | `supplementary.md` (top table) |
| how the files are built, their sources and checks | `scripts/README.md` |
| to rebuild everything after a change (e.g. the author names) | `bash scripts/build_all.sh` (about 15 s; refreshes `for_submission/` too) |
| to set the author names, affiliation and e-mail of the title blocks | `scripts/esm_common.py` (`AUTHORS`, `CORR_NAME`, `CORR_AFFILIATION`, `CORR_EMAIL`), then `build_all.sh` |
| the captions in the manuscript | last section of the working file, "Supplementary Information" |
| the status of each citation in the manuscript | all resolved: the 37 "supmat" Word comments that tracked them were deleted on 2026-10-01 (rule in `../../CLAUDE.md`) |
| the decisions and their history | `../build/proposed/ESM0_QUESTIONS.md` (questions and answers), `../build/SUPPLEMENTARY_DONE_2026-10.md` (hand-over), end of `../build/REVIEW_2026-09-27_Results_Discussion.md` (log) |
| the supplement that is not cited in the text | `unreferenced/` |

## Folder map

```
supplementary/
├── README.md                this map
├── supplementary.md         final list (number, content, status) + the planning list kept for the record
├── ESM_1.xlsx … ESM_13.pdf  THE 13 ONLINE RESOURCES (numbered by first citation; built by scripts/)
├── for_submission/          copies for download / upload: xlsx/, pdf/, ESM_all_Online_Resources.zip, CAPTIONS.md
│                            (rewritten on every build: do not edit by hand)
├── unreferenced/            ESM_extra_phenotype_checks.{pdf,xlsx}: raw-phenotype checks, NOT cited, unnumbered
│                            (the user decides later whether a part joins the numbered set)
├── scripts/                 the builders (README.md there): esm_common.py, build_ESM_xlsx.py, build_ESM_pdf.py,
│                            build_ESM_extra.py, package_for_submission.py, fetch_ENA_report.sh, build_all.sh,
│                            make_ESM_variance_components.py, make_ESM_candidate_genes.py; tex/ = PDF page sources
├── sources/                 inputs that are not project outputs:
│   ├── ENA_PRJEB79623_read_run_2026-10-01.tsv   public ENA run report (ESM 1)
│   └── NIR_calibration/     N. Pintel's calibration files, copied from the user's NIR.zip (ESM 2; README there)
└── tables/                  intermediate checked tables (variance components, candidate genes) read by the builders
```

## Where the content comes from

Every value is copied from a read-only project output; nothing is re-estimated (except the Pearson P values of the
NIR-calibration sheet, recomputed from N. Pintel's data and identical to R's `cor.test`). The analysis folders behind
each file are listed in `scripts/README.md`. Supporting analyses written for the supplement:
`../../../12_USED_adittional_checks_for_supplementary/` (01 block-4 correlations → ESM 2; 02 raw-phenotype checks →
`unreferenced/`). The original NIR files: `../../section_materials/02_materials_methods/NIR_calibration_files/NIR.zip`.

## Rules that apply here (from `../../CLAUDE.md` and the user's decisions)

- One file per Online Resource, named `ESM_N`, numbered by first citation; each carries the TAG title block and, on its
  first sheet or page, the caption and definitions. No merged PDF.
- The relevance-screen categories (strong / plausible / unlikely / no annotation) never appear (checked by the builders).
- The +8.85 adjustment of the 2021–22 starch values is never mentioned in any supplementary file (user, 2026-10-01).
- While a supplementary citation is unresolved, it carries a "supmat" Word comment (❌); it is deleted once resolved. None is open
  since 2026-10-01.
- Delete before submission: the open-items list at the top of the manuscript.
