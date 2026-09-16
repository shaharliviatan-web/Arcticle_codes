#!/bin/bash
# =============================================================================
# 01b_preaggregate_ld.sh
#
# True single-pass aggregation of pairwise LD into 1-kb physical-distance bins.
# For each chromosome:
#   - zcat streams the .ld.gz through stdout (no temp file on disk)
#   - awk computes |BP_A - BP_B|, bins into 1-kb buckets up to 2.5 Mb
#     and accumulates SUM(r^2), SUM(r^2^2), and COUNT per bin
#   - sort sorts by bin
#   - output: a small ~50 KB TSV per chromosome
#
# These pre-aggregated tables are read by 02_compute_decay_profile.R, which
# then computes mean / SD / SE genome-wide, applies a pair-count-weighted
# LOESS smoother, and finds the r^2 = {0.1, 0.2, 0.3} crossings.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

INTER_DIR="${PROJECT_DIR}/intermediates"
LOG_DIR="${PROJECT_DIR}/logs"
mkdir -p "${LOG_DIR}"

MAX_DIST_BP=2500000   # 2.5 Mb (matches --ld-window-kb 2500 from step 1)
BIN_SIZE_BP=1000      # 1 kb bins

# --- choice of decompressor + awk ---------------------------------------------
# pigz is multi-threaded gunzip (much faster than single-thread zcat once awk
# is fast). mawk is typically 5-10x faster than gawk on simple stream tasks.
# Fall back to zcat / awk if the faster tools aren't installed.
if command -v pigz >/dev/null; then
    DECOMPRESS="pigz -dc -p 4"
else
    DECOMPRESS="zcat"
fi
if command -v mawk >/dev/null; then
    AWK_BIN="mawk"
else
    AWK_BIN="awk"
fi
echo "Using decompressor: ${DECOMPRESS}"
echo "Using awk binary:   $(command -v ${AWK_BIN})  ($(${AWK_BIN} -W version 2>&1 | head -1 || true))"

echo "============================================================"
echo "Pre-aggregating LD pairs into ${BIN_SIZE_BP}-bp bins, up to ${MAX_DIST_BP} bp."
echo "Started: $(date -Iseconds)"
echo "============================================================"

for CHR in 1H 2H 3H 4H 5H 6H 7H; do
    IN_FILE="${INTER_DIR}/ld_pairs_chr${CHR}.ld.gz"
    OUT_FILE="${INTER_DIR}/pre_aggregated_chr${CHR}.tsv"
    LOG_FILE="${LOG_DIR}/preagg_chr${CHR}.log"

    if [ ! -s "${IN_FILE}" ]; then
        echo "ERROR: missing ${IN_FILE}" >&2
        exit 2
    fi
    if [ -s "${OUT_FILE}" ]; then
        echo "[${CHR}] pre-aggregated file already exists - skipping"
        continue
    fi

    echo "[${CHR}] $(date -Iseconds) - aggregating ${IN_FILE}"

    # NR==1 skips PLINK header.
    # PLINK .ld columns: 1 CHR_A  2 BP_A  3 SNP_A  4 CHR_B  5 BP_B  6 SNP_B  7 R2
    # The PLINK output is space-separated with leading whitespace.
    /usr/bin/time -v -o "${LOG_FILE}" \
    bash -c "${DECOMPRESS} '${IN_FILE}' | ${AWK_BIN} -v MAXD=${MAX_DIST_BP} -v BSZ=${BIN_SIZE_BP} '
        NR == 1 { next }
        {
            d = (\$5 > \$2) ? \$5 - \$2 : \$2 - \$5;
            if (d == 0 || d > MAXD) next;
            bin = int(d / BSZ) * BSZ;
            sum[bin]   += \$7;
            sumsq[bin] += \$7 * \$7;
            n[bin]++;
        }
        END {
            OFS = \"\\t\";
            print \"Bin\", \"Sum_R2\", \"SumSq_R2\", \"N_Pairs\";
            for (b in n) print b, sum[b], sumsq[b], n[b];
        }
    ' | (head -1; sort -k1,1n) > '${OUT_FILE}'"

    SIZE=$(du -h "${OUT_FILE}" | cut -f1)
    ROWS=$(wc -l < "${OUT_FILE}")
    echo "[${CHR}] $(date -Iseconds) - done. Output: ${OUT_FILE} (${SIZE}, ${ROWS} rows)"
done

echo ""
echo "============================================================"
echo "All chromosomes pre-aggregated."
ls -lh "${INTER_DIR}"/pre_aggregated_chr*.tsv
echo "============================================================"
