#!/usr/bin/env bash
# run_all.sh - reproduce 06_USED_interpro_domains end to end.
# Idempotent: step 00 only downloads/installs if missing; step 01 skips sequences
# whose JSON result already exists; step 02 re-parses from scratch.
#
# Note: step 01 contacts the EBI InterProScan REST service (45 jobs, ~30 min).
# For a long run prefer launching it under GNU screen, e.g.:
#   screen -dmS iprscan bash -lc 'TMPDIR=/mnt/data/shahar/.tmp bash scripts/01_run_interproscan.sh'
set -euo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT=/usr/bin/Rscript

echo "### 00 setup + versions"
bash "$HERE/00_setup_and_versions.sh"
echo "### 01 InterProScan run (EBI REST)"
bash "$HERE/01_run_interproscan.sh"
echo "### 02 parse -> per-gene table"
"$RSCRIPT" "$HERE/02_parse_interpro.R"
echo "### done"
