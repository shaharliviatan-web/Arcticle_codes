#!/usr/bin/env bash
# run_all.sh - reproduce the whole 05_USED_gene_annotation step end to end.
# Each script logs to ../logs/. Safe to re-run: downloads and the proteome copy
# are idempotent.
set -euo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(dirname "$HERE")"
RSCRIPT=/usr/bin/Rscript
LOGS="$BASE/logs"
mkdir -p "$LOGS"

echo "### 00 lock FDR genes"
"$RSCRIPT" "$HERE/00_lock_fdr_genes.R"        2>&1 | tee "$LOGS/00_lock_fdr_genes.log"
echo "### 01 build protein FASTA"
"$RSCRIPT" "$HERE/01_make_protein_fasta.R"    2>&1 | tee "$LOGS/01_make_protein_fasta.log"
echo "### 02 Swiss-Prot DIAMOND blastp"
bash      "$HERE/02_swissprot_diamond.sh"     2>&1 | tee "$LOGS/02_swissprot_diamond.log"
echo "### 03 build annotation table"
"$RSCRIPT" "$HERE/03_build_annotation_table.R" 2>&1 | tee "$LOGS/03_build_annotation_table.log"

echo "### done"
