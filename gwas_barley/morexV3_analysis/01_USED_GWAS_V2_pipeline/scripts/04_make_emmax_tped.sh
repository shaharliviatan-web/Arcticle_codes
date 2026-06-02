#!/usr/bin/env bash
# 04_make_emmax_tped.sh
# Recode the 290-sample PLINK trio into the FULL 7.1M-SNP tped/tfam that EMMAX scans.
# This is the genotype matrix the 24 EMMAX runs all share.

set -euo pipefail

# ---- Environment ----
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

# ---- Paths ----
PIPE_ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/01_USED_GWAS_V2_pipeline
INTER_DIR="$PIPE_ROOT/intermediates"
LOG_DIR="$PIPE_ROOT/logs"

IN_PREFIX="$INTER_DIR/morexV3_290"
OUT_PREFIX="$INTER_DIR/morexV3_gwas_genotypes"
LOG_FILE="$LOG_DIR/step04_make_emmax_tped.log"

# ---- Tool path (absolute) ----
PLINK=/usr/local/bin/plink

# ---- Parameters ----
THREADS=32
EXPECTED_SNPS=7110996
EXPECTED_SAMPLES=290

# ---- Idempotency gate ----
if [[ -f "${OUT_PREFIX}.tped" && -f "${OUT_PREFIX}.tfam" ]]; then
    echo "[04] Outputs already exist at ${OUT_PREFIX}.{tped,tfam} -- skipping."
    echo "[04] Delete those files to force re-run."
else
    [[ -f "${IN_PREFIX}.bed" ]] || { echo "[04] FAIL: missing ${IN_PREFIX}.bed (run step 01)." >&2; exit 1; }

    echo "[04] Recoding to transposed 12 format (full 7.1M-SNP tped for EMMAX scan)..."
    echo "[04] IN_PREFIX  = $IN_PREFIX"
    echo "[04] OUT_PREFIX = $OUT_PREFIX"
    echo "[04] THREADS    = $THREADS"

    "$PLINK" \
        --bfile "$IN_PREFIX" \
        --recode12 transpose \
        --allow-extra-chr \
        --threads "$THREADS" \
        --out "$OUT_PREFIX" \
        2>&1 | tee "$LOG_FILE"
fi

# ---- Checkpoints ----
TPED_ROWS=$(wc -l < "${OUT_PREFIX}.tped")
TFAM_ROWS=$(wc -l < "${OUT_PREFIX}.tfam")

echo "[04] CHECKPOINT: tped rows = $TPED_ROWS  (expected $EXPECTED_SNPS, tolerance ±1%)"
echo "[04] CHECKPOINT: tfam rows = $TFAM_ROWS  (expected $EXPECTED_SAMPLES)"

if [[ "$TFAM_ROWS" -ne "$EXPECTED_SAMPLES" ]]; then
    echo "[04] FAIL: tfam row count ($TFAM_ROWS) != $EXPECTED_SAMPLES. Halting." >&2
    exit 1
fi

# tped tolerance: 7,039,886 .. 7,182,106 (±1%)
TPED_LO=7039886
TPED_HI=7182106
if [[ "$TPED_ROWS" -lt "$TPED_LO" || "$TPED_ROWS" -gt "$TPED_HI" ]]; then
    echo "[04] FAIL: tped row count ($TPED_ROWS) outside tolerance [$TPED_LO, $TPED_HI]. Halting." >&2
    exit 1
fi

# Sample-order consistency: this tfam must match the pruned tfam (and thus the kinship + covariates)
PRUNED_TFAM="$INTER_DIR/morexV3_pruned_covariates_data.tfam"
if [[ -f "$PRUNED_TFAM" ]]; then
    if ! diff -q <(awk '{print $1, $2}' "${OUT_PREFIX}.tfam") <(awk '{print $1, $2}' "$PRUNED_TFAM") >/dev/null; then
        echo "[04] FAIL: tfam sample order does NOT match pruned tfam. EMMAX requires identical order." >&2
        exit 1
    fi
    echo "[04] CHECKPOINT: tfam sample order matches pruned tfam (kinship + covariates aligned)"
fi

echo "[04] OK: full EMMAX tped ready at ${OUT_PREFIX}.{tped,tfam}"
ls -lh "${OUT_PREFIX}.tped" "${OUT_PREFIX}.tfam"
