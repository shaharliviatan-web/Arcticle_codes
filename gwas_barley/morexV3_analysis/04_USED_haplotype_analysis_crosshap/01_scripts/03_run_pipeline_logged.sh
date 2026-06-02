#!/usr/bin/env bash
# 03_run_pipeline_logged.sh -- run the crosshap pipeline with logging.
# Usage: bash 03_run_pipeline_logged.sh [config.yaml] [gene_windows.tsv]
# Intended to run inside GNU screen for the full 87-gene job.

set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="${1:-$script_dir/../00_config/config.yaml}"
GW="${2:-$script_dir/../00_config/gene_windows.tsv}"

cfgval() { grep -E "^[[:space:]]*$1:" "$CFG" | head -1 | sed -E "s/^[[:space:]]*$1:[[:space:]]*\"?([^\"]*)\"?[[:space:]]*$/\1/"; }
run_id="$(cfgval run_id)"
output_root="$(cfgval output_root)"

log_dir="${output_root}/04_runs/${run_id}/Logs"
mkdir -p "$log_dir"
ts="$(date +'%Y%m%d_%H%M%S')"
log_file="${log_dir}/run_${run_id}_${ts}.log"

echo "Run ID:       ${run_id}"
echo "Config:       ${CFG}"
echo "Gene windows: ${GW}"
echo "Logging to:   ${log_file}"

/usr/bin/Rscript "$script_dir/03_run_crosshap_pipeline.R" "$CFG" "$GW" 2>&1 | tee "$log_file"
