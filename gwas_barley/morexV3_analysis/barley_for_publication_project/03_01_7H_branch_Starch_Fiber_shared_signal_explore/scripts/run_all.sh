#!/usr/bin/env bash
# run_all.sh - reproduce the whole replacement analysis (~3 min, plus the elite fetch).
set -euo pipefail
S="$(dirname "$(readlink -f "$0")")"; source "$S/config.sh"; mkdir -p "$TEMP_ROOT/logs"
Rscript "$S/01_run_mgmin3.R"            2>&1 | tee "$TEMP_ROOT/logs/01.log"
bash    "$S/02_fetch_elite_vcf.sh"      2>&1 | tee "$TEMP_ROOT/logs/02.log"
Rscript "$S/02_elite_comparison.R"      2>&1 | tee "$TEMP_ROOT/logs/02b.log"
Rscript "$S/03_violin_plus_barcode.R"   2>&1 | tee "$TEMP_ROOT/logs/03.log"
Rscript "$S/04_paper_tables.R"          2>&1 | tee "$TEMP_ROOT/logs/04.log"
Rscript "$S/05_results_chapter_numbers.R" 2>&1 | tee "$TEMP_ROOT/logs/05.log"
echo "done -> $TEMP_ROOT/results"
