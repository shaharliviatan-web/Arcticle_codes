#!/usr/bin/env bash
# run_all.sh - build the consolidated annotation master table (step 08).
# Pure local merge of the step 05/06/07 deliverables + the step-04 lead-SNP class;
# no network. Idempotent (re-parses from scratch).
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "### 00 build master table"
/usr/bin/Rscript "$HERE/00_build_master.R"
echo "### done"
