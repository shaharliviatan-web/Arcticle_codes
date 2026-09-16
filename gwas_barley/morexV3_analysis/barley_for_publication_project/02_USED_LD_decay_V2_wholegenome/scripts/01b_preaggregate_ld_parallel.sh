#!/bin/bash
# =============================================================================
# 01b_preaggregate_ld_parallel.sh
#
# Same job as 01b_preaggregate_ld.sh, but parallelised: pigz decompresses with
# multiple threads, the line stream is split across N mawk workers via
# `parallel --pipe`, and a final mawk pass sums the per-bin partials.
#
# Math: mean r^2 per bin is a simple SUM / COUNT, and SUM/COUNT decompose
# trivially across stream chunks. Each worker produces partial SUM, SUM-OF-
# SQUARES, COUNT per bin; the aggregator adds them up. Result is identical
# to the serial version (modulo float-ordering noise <1e-15).
#
# Output format (matches the serial version exactly, so 02_compute_decay_profile.R
# does not need to change):
#   pre_aggregated_chr<CHR>.tsv with columns: Bin, Sum_R2, SumSq_R2, N_Pairs
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

INTER_DIR="${PROJECT_DIR}/intermediates"
LOG_DIR="${PROJECT_DIR}/logs"
mkdir -p "${LOG_DIR}"

BIN_PAIRS_AWK="${SCRIPT_DIR}/bin_pairs.awk"
AGGREGATE_AWK="${SCRIPT_DIR}/aggregate_partials.awk"

MAX_DIST_BP=2500000
BIN_SIZE_BP=1000
PIGZ_THREADS="${PIGZ_THREADS:-4}"     # cores for parallel decompression
AWK_WORKERS="${AWK_WORKERS:-8}"        # mawk workers fed by `parallel --pipe`
PARALLEL_BLOCK="${PARALLEL_BLOCK:-100M}"  # chunk size for `parallel --pipe`

# --- sanity ---
for f in "${BIN_PAIRS_AWK}" "${AGGREGATE_AWK}"; do
    [ -s "$f" ] || { echo "ERROR: missing awk script $f"; exit 1; }
done
for cmd in pigz parallel mawk; do
    command -v "$cmd" >/dev/null || { echo "ERROR: $cmd not in PATH"; exit 1; }
done

echo "============================================================"
echo "Pre-aggregating LD pairs into ${BIN_SIZE_BP}-bp bins, up to ${MAX_DIST_BP} bp."
echo "Mode:           parallel"
echo "pigz threads:   ${PIGZ_THREADS}"
echo "mawk workers:   ${AWK_WORKERS}"
echo "parallel block: ${PARALLEL_BLOCK}"
echo "Started:        $(date -Iseconds)"
echo "============================================================"

for CHR in 1H 2H 3H 4H 5H 6H 7H; do
    IN_FILE="${INTER_DIR}/ld_pairs_chr${CHR}.ld.gz"
    OUT_FILE="${INTER_DIR}/pre_aggregated_chr${CHR}.tsv"
    LOG_FILE="${LOG_DIR}/preagg_parallel_chr${CHR}.log"

    if [ ! -s "${IN_FILE}" ]; then
        echo "ERROR: missing ${IN_FILE}" >&2
        exit 2
    fi
    if [ -s "${OUT_FILE}" ]; then
        echo "[${CHR}] pre-aggregated file already exists - skipping"
        continue
    fi

    echo "[${CHR}] $(date -Iseconds) - aggregating (parallel) ${IN_FILE}"

    # Pipeline:
    #   pigz decompresses with ${PIGZ_THREADS} threads
    #   tail strips the PLINK header line
    #   parallel --pipe splits the stream into ${PARALLEL_BLOCK}-sized chunks,
    #     each fed to one of ${AWK_WORKERS} mawk processes running bin_pairs.awk
    #   Each worker emits its per-bin partials (no header)
    #   Final mawk pass adds the partials and emits the proper header + sorted body
    /usr/bin/time -v -o "${LOG_FILE}" \
    bash -c "
        set -o pipefail
        pigz -dc -p ${PIGZ_THREADS} '${IN_FILE}' \
          | tail -n +2 \
          | parallel --pipe --block ${PARALLEL_BLOCK} -j ${AWK_WORKERS} \
                mawk -v MAXD=${MAX_DIST_BP} -v BSZ=${BIN_SIZE_BP} -f '${BIN_PAIRS_AWK}' \
          | mawk -f '${AGGREGATE_AWK}' \
          | (head -1; sort -k1,1n) \
          > '${OUT_FILE}'
    "

    SIZE=$(du -h "${OUT_FILE}" | cut -f1)
    ROWS=$(wc -l < "${OUT_FILE}")
    echo "[${CHR}] $(date -Iseconds) - done. Output: ${OUT_FILE} (${SIZE}, ${ROWS} rows)"
done

echo ""
echo "============================================================"
echo "All chromosomes pre-aggregated (parallel)."
ls -lh "${INTER_DIR}"/pre_aggregated_chr*.tsv
echo "============================================================"
