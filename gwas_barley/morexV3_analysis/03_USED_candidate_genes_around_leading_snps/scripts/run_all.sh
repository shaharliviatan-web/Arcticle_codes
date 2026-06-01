#!/bin/bash
# run_all.sh -- run the candidate-gene pipeline end to end.
#   01_build_loci.R   curated SNPs -> single-linkage loci (200 kb) -> interval BED
#   02_extract_genes.sh  bedtools intersect intervals x GFF genes (1H-7H)
#   03_build_tables.R    parse genes, distances, assemble output tables
# Logs go to logs/; results to results/tables/.

set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/03_USED_candidate_genes_around_leading_snps
S="$BASE/scripts"
L="$BASE/logs"
ts() { date +"%Y-%m-%d %H:%M:%S"; }

echo "[$(ts)] 01_build_loci.R"
Rscript "$S/01_build_loci.R"        2>&1 | tee "$L/01_build_loci.log"
echo "[$(ts)] 02_extract_genes.sh"
bash    "$S/02_extract_genes.sh"    2>&1 | tee "$L/02_extract_genes.log"
echo "[$(ts)] 03_build_tables.R"
Rscript "$S/03_build_tables.R"      2>&1 | tee "$L/03_build_tables.log"
echo "[$(ts)] Done. Tables in $BASE/results/tables/"
