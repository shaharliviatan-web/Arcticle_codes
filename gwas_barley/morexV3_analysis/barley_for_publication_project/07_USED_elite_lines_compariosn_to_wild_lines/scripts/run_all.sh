#!/bin/bash
# =============================================================================
# run_all.sh -- the whole step, end to end
# =============================================================================
# Steps 00 and 01 hit public APIs and are both RESUMABLE: 00 caches one JSON per
# BioSamples accession, 01 reuses a padded pool download that is complete and
# still brackets its target window. A warm re-run is seconds; a cold run is
# dominated by the three DivBrowse pool exports.
#
#   FORCE_REFETCH=1 bash scripts/run_all.sh   # ignore cached downloads
#
# Author : Shahar Liviatan
# Created: 2026-09-10   Screen step added and scripts renumbered: 2026-09-14
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/params.sh"
mkdir -p "$DIR_LOGS"

RUN_LOG="${DIR_LOGS}/run_all.log"
: > "$RUN_LOG"

step () {
  local name="$1"; shift
  echo "=== ${name} ===" | tee -a "$RUN_LOG"
  if "$@" >> "$RUN_LOG" 2>&1; then
    echo "    OK" | tee -a "$RUN_LOG"
  else
    echo "    FAILED -- see ${RUN_LOG}" | tee -a "$RUN_LOG"; exit 1
  fi
}

step "00 panel metadata (DivBrowse + EBI BioSamples; verifies configured lines)" \
     bash "${DIR_SCRIPTS}/00_panel_metadata.sh"
step "01 fetch elite candidate pool (DivBrowse API, padded + asserted + trimmed)" \
     bash "${DIR_SCRIPTS}/01_fetch_elite_vcfs.sh"
step "02 screen pool call rate (asserts configured lines pass)" \
     "$RSCRIPT_BIN" "${DIR_SCRIPTS}/02_screen_elite_lines.R"
step "03 build matrices (consensus, allele check, triallelic, 3 versions)" \
     "$RSCRIPT_BIN" "${DIR_SCRIPTS}/03_build_matrices.R"
step "04 figures (violin + aligned barcodes)" \
     "$RSCRIPT_BIN" "${DIR_SCRIPTS}/04_figures.R"
step "05 paper tables" \
     "$RSCRIPT_BIN" "${DIR_SCRIPTS}/05_paper_tables.R"

echo
echo "DONE. Figures  -> ${DIR_FIGURES}/<version>/"
echo "      Tables   -> ${DIR_TABLES}/"
echo "      Prose    -> ${DIR_TABLES}/results_chapter_numbers.txt"
echo "      Full log -> ${RUN_LOG}"
