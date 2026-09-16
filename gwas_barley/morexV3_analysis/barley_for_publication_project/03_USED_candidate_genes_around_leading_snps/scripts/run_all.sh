#!/bin/bash
# ============================================================================
# run_all.sh -- the whole of step 03, end to end. ~10 seconds.
#
#   01_build_intervals.R    step-01 locus table -> search intervals (BED)
#   02_extract_genes.sh     bedtools intersect intervals x Morex V3 genes
#   03_build_tables.R       parse, annotate, write every output table
#   04_flank_sensitivity.R  what other flank sizes would have returned
#   05_check_publication_genes.R  regression check against step 04's curated list
#
# Logs -> logs/  (one per step, plus run_all.log for the whole run)
# Tables -> results/tables/
#
# To change the search rule, edit config/params.sh and re-run this. Nothing in
# this pipeline depends on a value written anywhere else.
# ============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/../config/params.sh"

S="$STEP03_BASE/scripts"
L="$STEP03_BASE/logs"
mkdir -p "$L" "$STEP03_BASE/results/tables" "$STEP03_BASE/intermediates" "$STEP03_BASE/inputs"

ts() { date +"%Y-%m-%d %H:%M:%S"; }

{
echo "############################################################"
echo "# step 03 -- candidate genes at the GWAS loci"
echo "# started  : $(ts)"
echo "# host     : $(hostname)"
echo "# FLANK_BP : $FLANK_BP"
echo "# R        : $(Rscript -e 'cat(R.version.string)')"
echo "############################################################"

echo; echo "[$(ts)] 01_build_intervals.R"
Rscript "$S/01_build_intervals.R"   2>&1 | tee "$L/01_build_intervals.log"

echo; echo "[$(ts)] 02_extract_genes.sh"
bash    "$S/02_extract_genes.sh"    2>&1 | tee "$L/02_extract_genes.log"

echo; echo "[$(ts)] 03_build_tables.R"
Rscript "$S/03_build_tables.R"      2>&1 | tee "$L/03_build_tables.log"

echo; echo "[$(ts)] 04_flank_sensitivity.R"
Rscript "$S/04_flank_sensitivity.R" 2>&1 | tee "$L/04_flank_sensitivity.log"

echo; echo "[$(ts)] 05_check_publication_genes.R"
Rscript "$S/05_check_publication_genes.R" 2>&1 | tee "$L/05_check_publication_genes.log"

echo; echo "[$(ts)] DONE. Tables in $STEP03_BASE/results/tables/"
ls -1 "$STEP03_BASE/results/tables/"
} 2>&1 | tee "$L/run_all.log"
