#!/usr/bin/env bash
# run_all.sh — rebuild every output of this step, in order (~6 min).
# Read-only on the rest of the project; writes only inside this folder.
# Order matters: 02 -> 03 (independent) -> 04 (needs 02) -> 05, 06 (need 02; 06 needs 04)
# -> 07 -> 08, 09 (need T02, T04) -> 10 (needs every table).
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp TEMP=/mnt/data/shahar/.tmp TMP=/mnt/data/shahar/.tmp
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

bash scripts/01_extract_lead_genotypes.sh
for s in 02_lead_carrier_origin 03_haplotype_group_origin 04_matched_background 05_lead_within_site \
         06_lead_carrier_sharing 07_kinship_by_region 08_figure_carrier_origin 09_figure_desert_enrichment \
         10_results_numbers; do
  Rscript "scripts/$s.R" > "logs/$s.log" 2>&1 || { echo "[run_all] FAILED at $s (logs/$s.log)"; exit 1; }
  tail -n 1 "logs/$s.log"
done
echo "[run_all] done"
