#!/usr/bin/env bash
# run_all.sh - reproduce 07_USED_BLASTP_genes_with_no_annotation_left end to end.
# Idempotent: step 00 re-extracts the 5 proteins; step 01 skips if the .done
# marker exists; step 02 re-parses from scratch.
#
# Step 01 contacts NCBI (remote blastp vs nr); for a long search prefer launching
# it under GNU screen:
#   screen -dmS nrblast bash -lc 'TMPDIR=/mnt/data/shahar/.tmp bash scripts/01_blastp_nr_remote.sh'
set -euo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT=/usr/bin/Rscript

echo "### 00 extract residual proteins"
bash "$HERE/00_extract_residual_proteins.sh"
echo "### 01 remote BLASTP vs nr (Viridiplantae)"
bash "$HERE/01_blastp_nr_remote.sh"
echo "### 02 build deliverable table"
"$RSCRIPT" "$HERE/02_build_table.R"
echo "### done"
