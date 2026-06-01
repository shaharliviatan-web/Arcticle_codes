#!/usr/bin/env bash
# 01_vcf_to_plink.sh
# VCF -> canonical 290-sample PLINK trio for the GWAS V2 rebuild.
# Applies sample removal (10 site-04 IDs) and the locked filters (MAF=1e-6, GENO=1)
# on the same PLINK invocation, mirroring the working pattern in the current pipeline.

set -euo pipefail

# ---- Environment ----
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

# ---- Paths ----
PROJECT_ROOT=/mnt/data/shahar/gwas_barley
PIPE_ROOT="$PROJECT_ROOT/morexV3_analysis/USED_GWAS_V2_pipeline"
INTER_DIR="$PIPE_ROOT/intermediates"
LOG_DIR="$PIPE_ROOT/logs"

VCF_IN="$PROJECT_ROOT/data/inputs/morexV3_with_ids.vcf.gz"
SAMPLES_REMOVE="$PROJECT_ROOT/data/inputs/samples_to_remove_V3.txt"

OUT_PREFIX="$INTER_DIR/morexV3_290"
SNP_MAP_COPY="$INTER_DIR/snp_map.bim"
LOG_FILE="$LOG_DIR/step01_vcf_to_plink.log"

# ---- Tool paths (absolute; bare `plink` resolves to v0.76 in non-interactive bash on fidel) ----
PLINK=/usr/local/bin/plink

# ---- Locked filter parameters ----
MAF=0.000001
GENO=1
THREADS=32

# ---- Idempotency gate ----
if [[ -f "${OUT_PREFIX}.bed" && -f "${OUT_PREFIX}.bim" && -f "${OUT_PREFIX}.fam" ]]; then
    echo "[01] Outputs already exist at ${OUT_PREFIX}.{bed,bim,fam} -- skipping PLINK conversion."
    echo "[01] Delete those files to force re-run."
else
    echo "[01] Building canonical 290-sample PLINK trio from VCF..."
    echo "[01] VCF_IN          = $VCF_IN"
    echo "[01] SAMPLES_REMOVE  = $SAMPLES_REMOVE"
    echo "[01] OUT_PREFIX      = $OUT_PREFIX"
    echo "[01] MAF / GENO      = $MAF / $GENO"
    echo "[01] THREADS         = $THREADS"

    "$PLINK" \
        --vcf "$VCF_IN" \
        --remove "$SAMPLES_REMOVE" \
        --maf "$MAF" \
        --geno "$GENO" \
        --allow-extra-chr \
        --make-bed \
        --threads "$THREADS" \
        --out "$OUT_PREFIX" \
        2>&1 | tee "$LOG_FILE"
fi

# ---- Snapshot snp_map for downstream Manhattan / summary scripts ----
if [[ ! -f "$SNP_MAP_COPY" ]]; then
    cp "${OUT_PREFIX}.bim" "$SNP_MAP_COPY"
    echo "[01] Copied ${OUT_PREFIX}.bim -> $SNP_MAP_COPY"
fi

# ---- Checkpoints ----
FAM_ROWS=$(wc -l < "${OUT_PREFIX}.fam")
BIM_ROWS=$(wc -l < "${OUT_PREFIX}.bim")

echo "[01] CHECKPOINT: fam rows = $FAM_ROWS  (expected 290)"
echo "[01] CHECKPOINT: bim rows = $BIM_ROWS  (expected ~7,110,996)"

if [[ "$FAM_ROWS" -ne 290 ]]; then
    echo "[01] FAIL: fam row count ($FAM_ROWS) != 290. Halting." >&2
    exit 1
fi

# Allow 1% tolerance around 7,110,996 (~7,039,886 .. 7,182,106)
BIM_LO=7039886
BIM_HI=7182106
if [[ "$BIM_ROWS" -lt "$BIM_LO" || "$BIM_ROWS" -gt "$BIM_HI" ]]; then
    echo "[01] FAIL: bim row count ($BIM_ROWS) outside tolerance [$BIM_LO, $BIM_HI]. Halting." >&2
    exit 1
fi

echo "[01] OK: 290-sample PLINK trio ready at $OUT_PREFIX"
