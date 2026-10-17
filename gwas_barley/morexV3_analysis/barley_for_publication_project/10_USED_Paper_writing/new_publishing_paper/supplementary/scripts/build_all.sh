#!/usr/bin/env bash
# build_all.sh -- rebuild every Online Resource of the TAG manuscript (supplementary/ESM_1.xlsx ... ESM_13.pdf).
# Order: the two checked tables first (their TSVs feed ESM 3, 5, 6, 7), then the spreadsheets, then the PDFs.
# Re-run after any change to an input or to the title block (esm_common.py: authors, affiliation, e-mail).
# Created 2026-10-01 (supplementary session).
set -euo pipefail
cd "$(dirname "$0")"
export TMPDIR=/mnt/data/shahar/.tmp
PY=python3
$PY make_ESM_variance_components.py
$PY make_ESM_candidate_genes.py
$PY build_ESM_xlsx.py
$PY build_ESM_pdf.py
# the unnumbered supplement (not cited in the text; inputs from 12_USED_adittional_checks_for_supplementary/0{1,2}_*)
$PY build_ESM_extra.py
# copies for download / upload (../for_submission/: xlsx/, pdf/, zip, CAPTIONS.md)
$PY package_for_submission.py
