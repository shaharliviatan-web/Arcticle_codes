#!/usr/bin/env bash
# fetch_ENA_report.sh -- download the public ENA run report of project PRJEB79623 (the raw reads of the 300
# accessions, Potapenko et al. 2026b, Mol Ecol) into ../sources/, dated. Read-only public query (ENA Portal API).
# Used by build_ESM_xlsx.py (ESM 1): verifies that every accession ID of the study (= ENA sample alias, e.g. 01_03)
# exists in PRJEB79623 and gives its ENA sample and run accessions.
# Run 2026-10-01: 300 runs, 300 aliases; all 290 analysed accessions present; the other 10 are site 04 (04_01-04_11).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/../sources/ENA_PRJEB79623_read_run_$(date +%F).tsv"
curl -s -m 120 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB79623&result=read_run&fields=run_accession,sample_accession,secondary_sample_accession,sample_alias,sample_title&format=tsv" -o "$OUT"
echo "written $OUT ($(($(wc -l < "$OUT") - 1)) runs)"
