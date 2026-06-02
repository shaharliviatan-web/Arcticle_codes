#!/usr/bin/env bash
# 05_run_emmax_grid.sh
# Run the 24-cell EMMAX sensitivity grid:
#   4 traits  x  {BLUP, BLUE} phenotype correction  x  {3, 5, 10} PCs
# All 24 cells reuse the same kinship matrix and same full 7.1M-SNP tped.
# Per-cell .ps output is written to results/emmax_ps/.

set -euo pipefail

# ---- Environment ----
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

# ---- BLAS thread cap (EMMAX is BLAS-linked: by default uses 7-9 threads/process) ----
# Cap to 4 threads per EMMAX so total core usage = JOBS x 4 stays well under 100 (fidel is shared)
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

# ---- Paths ----
PIPE_ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/01_USED_GWAS_V2_pipeline
INTER_DIR="$PIPE_ROOT/intermediates"
COV_DIR="$INTER_DIR/covariates"
LOG_DIR="$PIPE_ROOT/logs"
PS_DIR="$PIPE_ROOT/results/emmax_ps"
mkdir -p "$PS_DIR"

INPUTS_DIR=/mnt/data/shahar/gwas_barley/data/inputs
TPED_PREFIX="$INTER_DIR/morexV3_gwas_genotypes"
KINSHIP="$INTER_DIR/morexV3_kinship.aIBS.kinf"

JOBLOG="$LOG_DIR/step05_joblog.tsv"
MASTER_LOG="$LOG_DIR/step05_emmax_grid.log"

# ---- Tool path ----
EMMAX=/mnt/data/shahar/gwas_barley/tools/emmax-intel64

# ---- Parallelism (overridable via env) ----
# Default JOBS=20 with OMP=4 -> 80 cores max (safe under 100-core shared-server cap).
# 20 jobs in 1st wave + 4 in 2nd wave; ~25-40 min total estimated.
JOBS=${JOBS:-20}

# ---- Pre-flight: input files exist ----
[[ -f "${TPED_PREFIX}.tped" ]] || { echo "[05] FAIL: missing ${TPED_PREFIX}.tped (run step 04)." >&2; exit 1; }
[[ -f "${TPED_PREFIX}.tfam" ]] || { echo "[05] FAIL: missing ${TPED_PREFIX}.tfam (run step 04)." >&2; exit 1; }
[[ -f "$KINSHIP" ]]             || { echo "[05] FAIL: missing $KINSHIP (run step 03)." >&2; exit 1; }
for N in 3 5 10; do
    [[ -f "$COV_DIR/pca_pc${N}.cov" ]] || { echo "[05] FAIL: missing $COV_DIR/pca_pc${N}.cov (run step 03)." >&2; exit 1; }
done

# ---- Pre-flight: RAM guard (each EMMAX run uses ~4-8 GB; budget 10 GB/job + 20 GB buffer) ----
FREE_GB=$(free -g | awk '/^Mem:/ {print $7}')
RAM_NEEDED=$(( JOBS * 10 + 20 ))
if [[ "$FREE_GB" -lt "$RAM_NEEDED" ]]; then
    echo "[05] FAIL: only ${FREE_GB} GB available RAM; need >= ${RAM_NEEDED} GB for $JOBS concurrent EMMAX jobs." >&2
    echo "[05] Override with: JOBS=N $0  (or wait for other jobs to finish)" >&2
    exit 1
fi
echo "[05] RAM check: ${FREE_GB} GB available (>= ${RAM_NEEDED} GB needed for $JOBS concurrent jobs)"
echo "[05] Core budget: $JOBS jobs x $OMP_NUM_THREADS threads = $(( JOBS * OMP_NUM_THREADS )) cores"

# ---- Pre-flight: sample ID order in all 8 phenotype files matches tfam (EMMAX requires identical order) ----
declare -a PHENO_FILES=(
    "$INPUTS_DIR/protein_corrected_V3.pheno"
    "$INPUTS_DIR/starch_corrected_V3.pheno"
    "$INPUTS_DIR/fiber_corrected_V3.pheno"
    "$INPUTS_DIR/betaglucan_corrected_V3.pheno"
    "$INPUTS_DIR/BLUE_protein_V3.pheno"
    "$INPUTS_DIR/BLUE_starch_V3.pheno"
    "$INPUTS_DIR/BLUE_fiber_V3.pheno"
    "$INPUTS_DIR/BLUE_betaglucan_V3.pheno"
)
for pf in "${PHENO_FILES[@]}"; do
    [[ -f "$pf" ]] || { echo "[05] FAIL: missing phenotype file $pf" >&2; exit 1; }
    if ! diff -q <(awk '{print $1, $2}' "${TPED_PREFIX}.tfam") <(awk '{print $1, $2}' "$pf") >/dev/null; then
        echo "[05] FAIL: $pf sample order != tfam order (EMMAX requires identical order)." >&2
        exit 1
    fi
done
echo "[05] Sample-order check: all 8 phenotype files align with ${TPED_PREFIX}.tfam"

# ====================================================================
# Generate the 24-cell job list, skipping cells whose .ps already exists
# Columns: trait  pheno_type  N_PCs  pheno_file  out_prefix  cov_file
# ====================================================================
gen_jobs() {
    local trait pheno_type pheno_file N out_prefix cov
    for trait in protein starch fiber betaglucan; do
        for pheno_type in BLUP BLUE; do
            case "$pheno_type" in
                BLUP) pheno_file="$INPUTS_DIR/${trait}_corrected_V3.pheno" ;;
                BLUE) pheno_file="$INPUTS_DIR/BLUE_${trait}_V3.pheno" ;;
            esac
            for N in 3 5 10; do
                out_prefix="$PS_DIR/morexV3__${trait}__${pheno_type}__pc${N}"
                cov="$COV_DIR/pca_pc${N}.cov"
                if [[ -f "${out_prefix}.ps" ]]; then
                    echo "[05] Skipping (already exists): ${out_prefix}.ps" >&2
                    continue
                fi
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$trait" "$pheno_type" "$N" "$pheno_file" "$out_prefix" "$cov"
            done
        done
    done
}

JOBS_TSV=$(mktemp -p "$TMPDIR" step05_jobs.XXXXXX.tsv)
gen_jobs > "$JOBS_TSV"
N_TODO=$(wc -l < "$JOBS_TSV")

if [[ "$N_TODO" -eq 0 ]]; then
    echo "[05] All 24 .ps files already exist -- nothing to run."
else
    echo "[05] Launching $N_TODO EMMAX jobs (${JOBS} concurrent)..."
    echo "[05] tped     = ${TPED_PREFIX}"
    echo "[05] kinship  = ${KINSHIP}"
    echo "[05] joblog   = ${JOBLOG}"

    # Run via GNU parallel; each EMMAX writes its own .ps/.reml/.log alongside out_prefix
    # Cols: {1}=trait {2}=pheno_type {3}=N {4}=pheno_file {5}=out_prefix {6}=cov_file
    parallel -j "$JOBS" --colsep '\t' --joblog "$JOBLOG" --line-buffer \
        "$EMMAX" -v -d 10 \
            -t "$TPED_PREFIX" \
            -p "{4}" \
            -k "$KINSHIP" \
            -c "{6}" \
            -o "{5}" \
        < "$JOBS_TSV" \
        2>&1 | tee "$MASTER_LOG"
fi

rm -f "$JOBS_TSV"

# ====================================================================
# Checkpoints
# ====================================================================
echo "---------------------------------"
echo "[05] Verifying 24 .ps outputs..."

EXPECTED=24
N_PS=$(ls "$PS_DIR"/morexV3__*__*__pc*.ps 2>/dev/null | wc -l)
echo "[05] CHECKPOINT: .ps files present = $N_PS  (expected $EXPECTED)"

if [[ "$N_PS" -ne "$EXPECTED" ]]; then
    echo "[05] FAIL: expected $EXPECTED .ps files, found $N_PS. Halting." >&2
    ls "$PS_DIR"/morexV3__*.ps 2>/dev/null | sort >&2
    exit 1
fi

# Each .ps must have 7,110,996 rows and 4 columns
FAIL=0
for f in "$PS_DIR"/morexV3__*__*__pc*.ps; do
    rows=$(wc -l < "$f")
    cols=$(awk 'NR==1{print NF; exit}' "$f")
    if [[ "$rows" -ne 7110996 || "$cols" -ne 4 ]]; then
        echo "[05] FAIL: $f has ${rows}x${cols} (expected 7110996x4)" >&2
        FAIL=1
    fi
done
if [[ "$FAIL" -ne 0 ]]; then
    exit 1
fi
echo "[05] CHECKPOINT: all 24 .ps files = 7,110,996 rows x 4 cols (SNP, beta, SE, p)"

# Joblog: any nonzero exit code?
if [[ -f "$JOBLOG" ]]; then
    BAD=$(awk 'NR>1 && $7 != 0 {print $0}' "$JOBLOG" | wc -l)
    if [[ "$BAD" -ne 0 ]]; then
        echo "[05] FAIL: $BAD job(s) had nonzero exit code in $JOBLOG" >&2
        awk 'NR==1 || $7 != 0' "$JOBLOG" >&2
        exit 1
    fi
    echo "[05] CHECKPOINT: all joblog exit codes = 0"
fi

echo "[05] OK: 24 EMMAX runs complete; .ps files in $PS_DIR"
echo ""
echo "[05] Summary of .ps files:"
ls -1 "$PS_DIR"/morexV3__*__*__pc*.ps | sort | sed 's|.*/||'
