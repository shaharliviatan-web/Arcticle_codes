#!/usr/bin/env bash
# run_all.sh - reproduce 09_USED_canonical_betaglucan_gene_check end to end.
# Idempotent: refs are re-fetched (small), BLAST dbs only built if missing, the
# rest re-parses. Steps 02 (genome makeblastdb + tblastn) is the slow part;
# for a long run launch it under GNU screen (see 02 header).
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT=/usr/bin/Rscript

echo "### 01 fetch canonical references";        bash "$HERE/01_fetch_refs.sh"
echo "### 02 build dbs + tblastn/blastp";          bash "$HERE/02_makedb_and_blast.sh"
echo "### 03 reciprocal best hit";                 bash "$HERE/03_reciprocal_blast.sh"
echo "### 04 CslF/CslH alignment + NJ tree";       bash "$HERE/04_tree.sh"
echo "### 05 GFF3 description cross-check";         bash "$HERE/05_gff3_crosscheck.sh"
echo "### 06 build distance/assignment table";     "$RSCRIPT" "$HERE/06_build_table.R"
echo "### done"
