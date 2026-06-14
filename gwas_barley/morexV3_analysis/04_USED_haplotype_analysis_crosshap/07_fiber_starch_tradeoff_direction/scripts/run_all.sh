#!/usr/bin/env bash
# Run the fiber/starch tradeoff-direction test and regenerate all tables.
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
HERE="/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/07_fiber_starch_tradeoff_direction"
/usr/bin/Rscript "$HERE/scripts/01_tradeoff_direction.R" 2>&1 | tee "$HERE/logs/01_tradeoff_direction.log"
