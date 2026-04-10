#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

run_id="$(grep -E '^[[:space:]]*run_id:' config.yaml | head -n1 | sed -E 's/^[[:space:]]*run_id:[[:space:]]*"?([^"]+)"?/\1/')"
output_root="$(grep -E '^[[:space:]]*output_root:' config.yaml | head -n1 | sed -E 's/^[[:space:]]*output_root:[[:space:]]*"?([^"]+)"?/\1/')"

run_dir="${output_root}/Runs/${run_id}"
log_dir="${run_dir}/Logs"
mkdir -p "$log_dir"

ts="$(date +'%Y%m%d_%H%M%S')"
log_file="${log_dir}/run_${run_id}_${ts}.log"

echo "Logging to: ${log_file}"
echo "Run ID: ${run_id}"
echo "Output root: ${output_root}"

Rscript run_pipeline.R config.yaml targets.yaml 2>&1 | tee "$log_file"