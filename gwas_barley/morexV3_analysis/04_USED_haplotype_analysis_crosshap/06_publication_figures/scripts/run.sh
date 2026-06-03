#!/usr/bin/env bash
# Build the 4 publication haplotype figures + composite + results-chapter CSVs.
# Fast (seconds); no screen needed. TMPDIR pinned per workspace policy.
set -euo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$HERE/logs/make_figures_$(date +%Y%m%d_%H%M%S).log"

/usr/bin/Rscript "$HERE/scripts/00_make_figures.R" 2>&1 | tee "$LOG"
echo "Log: $LOG"
