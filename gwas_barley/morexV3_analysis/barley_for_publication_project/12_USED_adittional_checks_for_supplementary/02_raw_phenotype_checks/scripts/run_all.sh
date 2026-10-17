#!/usr/bin/env bash
# run_all.sh -- check 02 (raw phenotype checks for the supplement). Run from anywhere; logs in ../logs/.
set -euo pipefail
cd "$(dirname "$0")/.."
export TMPDIR=/mnt/data/shahar/.tmp
for s in 01_boxplots 02_correlations 03_distributions; do
  /usr/bin/Rscript "scripts/$s.R" > "logs/$s.log" 2>&1
  tail -n +1 "logs/$s.log" | grep -E '^\[|figure|table' || true
done
python3 scripts/04_results_numbers.py > logs/04_results_numbers.log 2>&1 && echo "  results_numbers.txt"
