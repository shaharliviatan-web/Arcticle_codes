#!/usr/bin/env bash
# 03_pca_kinship.sh
# Compute population-structure covariates ONCE on the LD-pruned ~590K SNP set:
#   (a) build pruned tped/tfam via --extract prune.in
#   (b) PLINK PCA -> 10 PCs (eigenvec + eigenval)
#   (c) EMMAX aIBS kinship matrix (290x290)
#   (d) subset eigenvec into pc3 / pc5 / pc10 covariate files for EMMAX -c
#   (e) scree-plot data TSV from eigenvalues
# These outputs are reused unchanged across all 24 EMMAX runs.

set -euo pipefail

# ---- Environment ----
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

# ---- Paths ----
PIPE_ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline
INTER_DIR="$PIPE_ROOT/intermediates"
COV_DIR="$INTER_DIR/covariates"
LOG_DIR="$PIPE_ROOT/logs"
mkdir -p "$COV_DIR"

IN_PREFIX="$INTER_DIR/morexV3_290"
PRUNE_IN="$INTER_DIR/morexV3_pruned_for_covs.prune.in"
PRUNED_PREFIX="$INTER_DIR/morexV3_pruned_covariates_data"
PCA_PREFIX="$INTER_DIR/morexV3_pca"
KINSHIP="$INTER_DIR/morexV3_kinship.aIBS.kinf"
SCREE_TSV="$INTER_DIR/morexV3_pca_scree_data.tsv"
LOG_FILE="$LOG_DIR/step03_pca_kinship.log"

# ---- Tool paths (absolute; non-interactive bash on fidel resolves bare `plink` to v0.76) ----
PLINK=/usr/local/bin/plink
EMMAX_KIN=/mnt/data/shahar/gwas_barley/tools/emmax-kin-intel64

# ---- Parameters ----
N_PCS_FULL=10
THREADS=32

# Fresh log for this step
: > "$LOG_FILE"

# ---- (a) Build pruned tped/tfam ----
if [[ -f "${PRUNED_PREFIX}.tped" && -f "${PRUNED_PREFIX}.tfam" ]]; then
    echo "[03a] Pruned tped/tfam already exist -- skipping." | tee -a "$LOG_FILE"
else
    [[ -f "${IN_PREFIX}.bed" ]] || { echo "[03a] FAIL: missing ${IN_PREFIX}.bed (run step 01)." >&2; exit 1; }
    [[ -f "$PRUNE_IN" ]]        || { echo "[03a] FAIL: missing $PRUNE_IN (run step 02)." >&2; exit 1; }
    echo "[03a] Building pruned tped via --extract prune.in ..." | tee -a "$LOG_FILE"
    "$PLINK" \
        --bfile "$IN_PREFIX" \
        --extract "$PRUNE_IN" \
        --recode12 transpose \
        --allow-extra-chr \
        --threads "$THREADS" \
        --out "$PRUNED_PREFIX" \
        2>&1 | tee -a "$LOG_FILE"
fi

# ---- (b) PCA: 10 PCs ----
if [[ -f "${PCA_PREFIX}.eigenvec" && -f "${PCA_PREFIX}.eigenval" ]]; then
    echo "[03b] PCA eigenvec/eigenval already exist -- skipping." | tee -a "$LOG_FILE"
else
    echo "[03b] Running PLINK --pca $N_PCS_FULL on pruned set ..." | tee -a "$LOG_FILE"
    "$PLINK" \
        --tfile "$PRUNED_PREFIX" \
        --pca "$N_PCS_FULL" \
        --allow-extra-chr \
        --threads "$THREADS" \
        --out "$PCA_PREFIX" \
        2>&1 | tee -a "$LOG_FILE"
fi

# ---- (c) aIBS kinship via EMMAX ----
# -s : IBS kinship (-> .aIBS.kinf);  -x : include non-autosomal chr (barley 1H-7H are "extra");  -d 10 : precision
if [[ -f "$KINSHIP" ]]; then
    echo "[03c] Kinship matrix already exists -- skipping." | tee -a "$LOG_FILE"
else
    echo "[03c] Building aIBS kinship (emmax-kin -v -s -d 10 -x) ..." | tee -a "$LOG_FILE"
    "$EMMAX_KIN" -v -s -d 10 -x "$PRUNED_PREFIX" 2>&1 | tee -a "$LOG_FILE"
    mv "${PRUNED_PREFIX}.aIBS.kinf" "$KINSHIP"
    echo "[03c] Renamed ${PRUNED_PREFIX}.aIBS.kinf -> $KINSHIP" | tee -a "$LOG_FILE"
fi

# ---- (d) Subset eigenvec -> pc3 / pc5 / pc10 covariate files ----
# eigenvec format: FID IID PC1 PC2 ... PC10 (no header). Keep FID IID + first N PCs.
for N in 3 5 10; do
    OUT="$COV_DIR/pca_pc${N}.cov"
    if [[ -f "$OUT" ]]; then
        echo "[03d] $OUT already exists -- skipping." | tee -a "$LOG_FILE"
    else
        awk -v n="$N" '{printf "%s %s", $1, $2; for (i=3; i<=2+n; i++) printf " %s", $i; printf "\n"}' \
            "${PCA_PREFIX}.eigenvec" > "$OUT"
        echo "[03d] Wrote $OUT" | tee -a "$LOG_FILE"
    fi
done

# ---- (e) Scree data TSV from eigenvalues ----
if [[ -f "$SCREE_TSV" ]]; then
    echo "[03e] Scree TSV already exists -- skipping." | tee -a "$LOG_FILE"
else
    awk 'BEGIN{print "PC\teigenval\tpct_variance\tcum_pct"}
         {v[NR]=$1; s+=$1}
         END{c=0; for(i=1;i<=NR;i++){p=100*v[i]/s; c+=p; printf "%d\t%.6f\t%.4f\t%.4f\n", i, v[i], p, c}}' \
        "${PCA_PREFIX}.eigenval" > "$SCREE_TSV"
    echo "[03e] Wrote $SCREE_TSV" | tee -a "$LOG_FILE"
fi

# ====================================================================
# Checkpoints
# ====================================================================
echo "---------------------------------" | tee -a "$LOG_FILE"

EIGVEC_ROWS=$(wc -l < "${PCA_PREFIX}.eigenvec")
EIGVEC_COLS=$(awk 'NR==1{print NF; exit}' "${PCA_PREFIX}.eigenvec")
EIGVAL_ROWS=$(wc -l < "${PCA_PREFIX}.eigenval")
KIN_ROWS=$(wc -l < "$KINSHIP")
KIN_COLS=$(awk 'NR==1{print NF; exit}' "$KINSHIP")

echo "[03] CHECKPOINT: eigenvec = ${EIGVEC_ROWS} x ${EIGVEC_COLS}  (expect 290 x 12)" | tee -a "$LOG_FILE"
echo "[03] CHECKPOINT: eigenval rows = ${EIGVAL_ROWS}  (expect 10)" | tee -a "$LOG_FILE"
echo "[03] CHECKPOINT: kinship  = ${KIN_ROWS} x ${KIN_COLS}  (expect 290 x 290)" | tee -a "$LOG_FILE"

[[ "$EIGVEC_ROWS" -eq 290 && "$EIGVEC_COLS" -eq 12 ]] || { echo "[03] FAIL: eigenvec shape != 290x12. Halting." >&2; exit 1; }
[[ "$EIGVAL_ROWS" -eq 10 ]]                            || { echo "[03] FAIL: eigenval rows != 10. Halting." >&2; exit 1; }
[[ "$KIN_ROWS" -eq 290 && "$KIN_COLS" -eq 290 ]]       || { echo "[03] FAIL: kinship shape != 290x290. Halting." >&2; exit 1; }

# Covariate-file shape checks
for N in 3 5 10; do
    OUT="$COV_DIR/pca_pc${N}.cov"
    R=$(wc -l < "$OUT")
    C=$(awk 'NR==1{print NF; exit}' "$OUT")
    EXP=$((2 + N))
    echo "[03] CHECKPOINT: $(basename "$OUT") = ${R} x ${C}  (expect 290 x ${EXP})" | tee -a "$LOG_FILE"
    [[ "$R" -eq 290 && "$C" -eq "$EXP" ]] || { echo "[03] FAIL: $OUT shape != 290x${EXP}. Halting." >&2; exit 1; }
done

# Sample-order consistency: covariate IDs must match pruned .tfam order (EMMAX requirement)
if ! diff -q <(awk '{print $1, $2}' "$COV_DIR/pca_pc10.cov") <(awk '{print $1, $2}' "${PRUNED_PREFIX}.tfam") >/dev/null; then
    echo "[03] FAIL: covariate sample order does not match ${PRUNED_PREFIX}.tfam. Halting." >&2
    exit 1
fi
echo "[03] CHECKPOINT: covariate sample order matches pruned .tfam (EMMAX-safe)" | tee -a "$LOG_FILE"

echo "[03] OK: PCA (10 PCs) + aIBS kinship + pc3/pc5/pc10 covariates + scree data ready" | tee -a "$LOG_FILE"
