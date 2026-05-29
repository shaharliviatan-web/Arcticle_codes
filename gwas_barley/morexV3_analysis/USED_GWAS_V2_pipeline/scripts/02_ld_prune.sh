#!/usr/bin/env bash
# 02_ld_prune.sh
# LD-prune the 290-sample bed/bim/fam with `--indep-pairwise 50 5 0.2`.
# Produces the ~590K independent SNP list used for PCA + kinship downstream.

set -euo pipefail

# ---- Environment ----
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

# ---- Paths ----
PIPE_ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline
INTER_DIR="$PIPE_ROOT/intermediates"
LOG_DIR="$PIPE_ROOT/logs"

IN_PREFIX="$INTER_DIR/morexV3_290"
OUT_PREFIX="$INTER_DIR/morexV3_pruned_for_covs"
LOG_FILE="$LOG_DIR/step02_ld_prune.log"

# ---- Tool paths (absolute; non-interactive bash on fidel resolves bare `plink` to v0.76) ----
PLINK=/usr/local/bin/plink

# ---- Locked LD-pruning parameters ----
WINDOW=50      # SNPs
STEP=5         # SNPs
R2=0.2
THREADS=32

# ---- Idempotency gate ----
if [[ -f "${OUT_PREFIX}.prune.in" && -f "${OUT_PREFIX}.prune.out" ]]; then
    echo "[02] Outputs already exist at ${OUT_PREFIX}.prune.{in,out} -- skipping LD pruning."
    echo "[02] Delete those files to force re-run."
else
    # Sanity-check input exists
    if [[ ! -f "${IN_PREFIX}.bed" ]]; then
        echo "[02] FAIL: missing input ${IN_PREFIX}.bed -- step 01 must run first." >&2
        exit 1
    fi

    echo "[02] LD pruning with --indep-pairwise $WINDOW $STEP $R2"
    echo "[02] IN_PREFIX   = $IN_PREFIX"
    echo "[02] OUT_PREFIX  = $OUT_PREFIX"
    echo "[02] THREADS     = $THREADS"

    "$PLINK" \
        --bfile "$IN_PREFIX" \
        --indep-pairwise "$WINDOW" "$STEP" "$R2" \
        --allow-extra-chr \
        --threads "$THREADS" \
        --out "$OUT_PREFIX" \
        2>&1 | tee "$LOG_FILE"
fi

# ---- Checkpoints ----
PRUNE_IN_ROWS=$(wc -l < "${OUT_PREFIX}.prune.in")
PRUNE_OUT_ROWS=$(wc -l < "${OUT_PREFIX}.prune.out")
TOTAL=$((PRUNE_IN_ROWS + PRUNE_OUT_ROWS))

echo "[02] CHECKPOINT: prune.in  rows = $PRUNE_IN_ROWS  (expected ~590,462; tolerance [580K, 600K])"
echo "[02] CHECKPOINT: prune.out rows = $PRUNE_OUT_ROWS"
echo "[02] CHECKPOINT: total kept+dropped = $TOTAL (sanity vs. 7,110,996 input)"

if [[ "$PRUNE_IN_ROWS" -lt 580000 || "$PRUNE_IN_ROWS" -gt 600000 ]]; then
    echo "[02] FAIL: prune.in row count ($PRUNE_IN_ROWS) outside tolerance [580000, 600000]. Halting." >&2
    exit 1
fi

# Recompute exact Bonferroni-pruned threshold for downstream Manhattan scripts
# v2: report both alpha=0.10 thresholds (the values actually used downstream in step 07/08).
PRUNED_BONF010=$(awk -v n="$PRUNE_IN_ROWS" 'BEGIN { printf "%.4f\n", -log(0.10 / n) / log(10) }')
ALL_BONF010=$(  awk          'BEGIN { printf "%.4f\n", -log(0.10 / 7110996) / log(10) }')
echo "[02] BonfAll010    = -log10(0.10 / 7,110,996)     = $ALL_BONF010 (expect ~7.85)"
echo "[02] BonfPruned010 = -log10(0.10 / $PRUNE_IN_ROWS) = $PRUNED_BONF010 (expect ~6.77)"
echo "[02] (expected ~7.07)"

echo "[02] OK: ${OUT_PREFIX}.prune.in ready ($PRUNE_IN_ROWS independent SNPs)"
