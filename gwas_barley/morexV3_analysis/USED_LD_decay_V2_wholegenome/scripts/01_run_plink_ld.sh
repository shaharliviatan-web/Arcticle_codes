#!/bin/bash
# =============================================================================
# 01_run_plink_ld.sh
#
# Genome-wide pairwise LD calculation for wild barley (morexV3, 290 accessions).
#
# - Uses the SAME PLINK fileset that the GWAS pipeline uses, so the SNP set is
#   identical to the one on the Manhattan plots.
# - Removes the 10 site-04 accessions (samples_to_remove_V3.txt) -> 290 samples.
# - Runs per chromosome (1H..7H) sequentially. Resumable: re-running skips
#   chromosomes whose .ld.gz already exists.
# - Output: intermediates/ld_pairs_chr{1H..7H}.ld.gz (+ .log) and a run-settings
#   file in logs/.
#
# Tune LDWIN_R2 below depending on free disk space (see comment).
# =============================================================================

set -euo pipefail

# --- 1. Paths -----------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PLINK_PREFIX="/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_data/intermediates_maf0.000001_geno1_thin0_pc5_LDfixed_withIDs/morexV3_barley_full_filtered"
SAMPLES_TO_REMOVE="/mnt/data/shahar/gwas_barley/data/inputs/samples_to_remove_V3.txt"

INTER_DIR="${PROJECT_DIR}/intermediates"
LOG_DIR="${PROJECT_DIR}/logs"
mkdir -p "${INTER_DIR}" "${LOG_DIR}"

# --- 2. Tuning ----------------------------------------------------------------
# r^2 truncation: keep all pairs (LDWIN_R2=0) for unbiased bin averages, or set
# to 0.001 / 0.01 if disk is tight. With --ld-window-kb 2500 and 7.1M SNPs,
# LDWIN_R2=0 can produce ~50-150 GB *per chromosome* (compressed).
#
# NOTE: the publication run used LDWIN_R2=0.001 (passed as an env var override),
# NOT the default of 0. This excludes pairs with r^2 < 0.001 from PLINK's
# output for storage reasons. The effect on the bin-averaged decay curve in
# the region of interest (<1 Mb, where the r^2 = 0.2 crossing lives) is
# negligible because most pairs there have r^2 >> 0.001. The actual value
# used is also recorded in logs/run_settings.txt for each run.
LDWIN_R2="${LDWIN_R2:-0}"

LD_WIN_KB=2500     # 2.5 Mb physical window
LD_WIN_SNPS=999999 # effectively unlimited SNPs per window
THREADS="${THREADS:-32}"

# --- 3. Sanity checks ---------------------------------------------------------
for ext in bed bim fam; do
    [ -s "${PLINK_PREFIX}.${ext}" ] || { echo "ERROR: missing ${PLINK_PREFIX}.${ext}"; exit 1; }
done
[ -s "${SAMPLES_TO_REMOVE}" ] || { echo "ERROR: missing ${SAMPLES_TO_REMOVE}"; exit 1; }
command -v plink >/dev/null || { echo "ERROR: plink not in PATH"; exit 1; }

# --- 4. Record run settings ---------------------------------------------------
SETTINGS_FILE="${LOG_DIR}/run_settings.txt"
{
    echo "Run started:   $(date -Iseconds)"
    echo "Hostname:      $(hostname)"
    echo "PLINK version: $(plink --version 2>&1 | head -1)"
    echo "PLINK_PREFIX:  ${PLINK_PREFIX}"
    echo "SAMPLES_REM:   ${SAMPLES_TO_REMOVE} ($(wc -l < "${SAMPLES_TO_REMOVE}") samples)"
    echo "LDWIN_R2:      ${LDWIN_R2}"
    echo "LD_WIN_KB:     ${LD_WIN_KB}"
    echo "LD_WIN_SNPS:   ${LD_WIN_SNPS}"
    echo "THREADS:       ${THREADS}"
    echo "Disk free:     $(df -h /mnt/data | tail -1)"
} | tee "${SETTINGS_FILE}"

# --- 5. Run PLINK per chromosome ----------------------------------------------
for CHR in 1H 2H 3H 4H 5H 6H 7H; do
    OUT_PREFIX="${INTER_DIR}/ld_pairs_chr${CHR}"
    LOG_FILE="${LOG_DIR}/plink_chr${CHR}.log"

    if [ -s "${OUT_PREFIX}.ld.gz" ]; then
        echo "[${CHR}] .ld.gz already exists — skipping. ($(date -Iseconds))"
        continue
    fi

    echo "============================================================"
    echo "[${CHR}] Starting PLINK r^2  ($(date -Iseconds))"
    echo "============================================================"

    plink --bfile "${PLINK_PREFIX}" \
          --remove "${SAMPLES_TO_REMOVE}" \
          --allow-extra-chr \
          --chr "${CHR}" \
          --maf 0.000001 \
          --r2 gz \
          --ld-window-r2 "${LDWIN_R2}" \
          --ld-window-kb "${LD_WIN_KB}" \
          --ld-window "${LD_WIN_SNPS}" \
          --threads "${THREADS}" \
          --out "${OUT_PREFIX}" \
        2>&1 | tee "${LOG_FILE}"

    if [ ! -s "${OUT_PREFIX}.ld.gz" ]; then
        echo "ERROR: ${OUT_PREFIX}.ld.gz was not produced for ${CHR}" >&2
        exit 2
    fi

    SIZE=$(du -h "${OUT_PREFIX}.ld.gz" | cut -f1)
    echo "[${CHR}] done — ${OUT_PREFIX}.ld.gz (${SIZE})  ($(date -Iseconds))"
done

echo ""
echo "============================================================"
echo "All chromosomes finished."
echo "Output files:"
ls -lh "${INTER_DIR}"/ld_pairs_chr*.ld.gz 2>/dev/null || true
echo "============================================================"
