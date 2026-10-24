#!/usr/bin/env bash
# build_all.sh -- rebuild the supplementary material of the TAG manuscript: ESM_1.xlsx (Supplementary Tables S1-S22),
# ESM_2.pdf (Supplementary Figures S1-S5), the unnumbered supplement, ../for_submission/ and ../for_review_Hubner/.
# Renumbered 2026-10-17 (S. Hübner; user decision): until then 13 Online Resources ESM_1 ... ESM_13.
# Order: the two checked tables first (their TSVs feed ESM 3, 5, 6, 7), then the spreadsheets, then the PDFs.
# Re-run after any change to an input or to the title block (esm_common.py: authors, affiliation, e-mail).
# Created 2026-10-01 (supplementary session).
set -euo pipefail
cd "$(dirname "$0")"
export TMPDIR=/mnt/data/shahar/.tmp
PY=python3
$PY make_ESM_variance_components.py
$PY make_ESM_candidate_genes.py
$PY build_ESM_xlsx.py          # content of the tables -> ../parts/part_N.xlsx (intermediate)
$PY build_ESM_tables.py        # ESM_1.xlsx, Tables S1-S22 (also checks the manuscript's citation order)
$PY build_ESM_pdf.py           # ESM_2.pdf, Figs. S1-S5
# the unnumbered supplement (not cited in the text; inputs from 12_USED_adittional_checks_for_supplementary/0{1,2}_*)
$PY build_ESM_extra.py
# copies for download / upload (../for_submission/: ESM_1.xlsx, ESM_2.pdf, zip, CAPTIONS.md)
$PY package_for_submission.py
# the review package for S. Hübner (../for_review_Hubner/; python-docx from /mnt/data/shahar/gwas_barley/tools/pylib)
$PY package_for_review.py
