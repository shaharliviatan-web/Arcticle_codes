#!/bin/bash
# ====================================================================
# GCTA Heritability Sensitivity Analysis Script - V3 Inbred
#
# Compares narrow-sense SNP heritability (h²) across three model
# variants to assess what drives differences in estimates:
#
#   Model 1: Inbred GRM, no PCs
#   Model 2: Inbred GRM + 2 PCs
#   Model 3: Inbred GRM + 5 PCs
#
# All models use the same inbred GRM and PCA run (computed once).
# Inputs are read-only (analysis_V3 folder is never modified).
#
# Output base: /mnt/data/shahar/GCTA/analysis_V3_inbred_sensitivity/
# ====================================================================

set -euo pipefail

echo "====================================================="
echo " GCTA Sensitivity Analysis: V3 Inbred, 0 / 2 / 5 PCs"
echo "====================================================="

# ====================================================================
# --- CONFIGURATION ---
# ====================================================================

# -- Multithreading --
THREADS=16

# -- Tool path --
GCTA="/mnt/data/shahar/GCTA/tools/gcta-1.95.0-linux-kernel-3-x86_64/gcta64"

# -- Read-only inputs --
PLINK_PREFIX_IN="/mnt/data/shahar/GCTA/analysis_V3/intermediates/barley_V3_filtered"
GCTA_INPUTS_DIR="/mnt/data/shahar/GCTA/inputs"

# -- Output directories --
ANALYSIS_DIR="/mnt/data/shahar/GCTA/analysis_V3_inbred_sensitivity"
INTERMEDIATES="${ANALYSIS_DIR}/intermediates"
RESULTS="${ANALYSIS_DIR}/results"
LOGS="${ANALYSIS_DIR}/logs"

mkdir -p "$INTERMEDIATES" "$RESULTS" "$LOGS"

# -- Phenotypes --
PHENOTYPES=( "Protein" "Starch" "Betaglucans" "Fiber" )

# -- Model definitions: label and number of PCs (0 = no qcovar) --
declare -A MODEL_PC
MODEL_PC["noPC"]=0
MODEL_PC["2PC"]=2
MODEL_PC["5PC"]=5
MODEL_ORDER=( "noPC" "2PC" "5PC" )

echo ""
echo "Configuration:"
echo "  GCTA binary : ${GCTA}"
echo "  PLINK input : ${PLINK_PREFIX_IN}"
echo "  Phenotypes  : ${PHENOTYPES[*]}"
echo "  Models      : ${MODEL_ORDER[*]}"
echo "  Threads     : ${THREADS}"
echo "  Output dir  : ${ANALYSIS_DIR}"
echo ""
echo "====================================================="

# ====================================================================
# --- STAGE 1: Build Inbred GRM (computed once) ---
#
# Uses --make-grm-inbred to account for the high homozygosity
# of selfing/inbred barley lines.
# --autosome-num 7 tells GCTA barley has 7 autosomes.
# --thread-num enables multithreading.
# ====================================================================
echo ""
echo "--- STAGE 1: Building Inbred GRM ---"
GRM_PREFIX="${INTERMEDIATES}/barley_V3_inbred_grm"

if [ ! -f "${GRM_PREFIX}.grm.bin" ]; then
    echo "  Running: --make-grm-inbred with ${THREADS} threads..."
    "${GCTA}" \
        --bfile "${PLINK_PREFIX_IN}" \
        --autosome-num 7 \
        --make-grm-inbred \
        --thread-num "${THREADS}" \
        --out "${GRM_PREFIX}" \
        2>&1 | tee "${LOGS}/stage1_make_grm_inbred.log"
    echo "  Inbred GRM created."
else
    echo "  Existing inbred GRM found. Skipping."
fi

echo "--- STAGE 1 COMPLETE ---"
echo "====================================================="

# ====================================================================
# --- STAGE 2: PCA on Inbred GRM (computed once, 5 PCs) ---
#
# Computes the top 5 PCs from the inbred GRM.
# The full eigenvec file is used for the 5-PC model.
# A trimmed version (FID, IID, PC1, PC2) is derived for the 2-PC model.
# The 0-PC model does not use any qcovar file.
# ====================================================================
echo ""
echo "--- STAGE 2: Running PCA (5 PCs) on Inbred GRM ---"
PCA_PREFIX="${INTERMEDIATES}/barley_V3_inbred_pca5"
EIGENVEC_5PC="${PCA_PREFIX}.eigenvec"
EIGENVEC_2PC="${INTERMEDIATES}/barley_V3_inbred_pca2.eigenvec"

if [ ! -f "${EIGENVEC_5PC}" ]; then
    echo "  Running: --pca 5..."
    "${GCTA}" \
        --grm "${GRM_PREFIX}" \
        --pca 5 \
        --thread-num "${THREADS}" \
        --out "${PCA_PREFIX}" \
        2>&1 | tee "${LOGS}/stage2_pca5.log"
    echo "  PCA complete: ${EIGENVEC_5PC}"
else
    echo "  Existing 5-PC eigenvec found. Skipping PCA."
fi

# Derive 2-PC qcovar file: keep FID, IID, PC1, PC2 (columns 1-4)
if [ ! -f "${EIGENVEC_2PC}" ]; then
    echo "  Generating 2-PC qcovar file from 5-PC eigenvec..."
    awk '{print $1, $2, $3, $4}' "${EIGENVEC_5PC}" > "${EIGENVEC_2PC}"
    echo "  2-PC qcovar file written: ${EIGENVEC_2PC}"
else
    echo "  Existing 2-PC qcovar file found. Skipping."
fi

echo "--- STAGE 2 COMPLETE ---"
echo "====================================================="

# ====================================================================
# --- STAGE 3: REML Heritability for All Traits × All Models ---
#
# For each phenotype and each model variant:
#   - noPC : --reml only (no --qcovar)
#   - 2PC  : --reml + --qcovar with 2-PC eigenvec
#   - 5PC  : --reml + --qcovar with 5-PC eigenvec
#
# Output naming: result_<Trait>_h2_inbred_<Model>.hsq
# Each run is logged separately in logs/.
# ====================================================================
echo ""
echo "--- STAGE 3: REML Heritability Analysis ---"

for MODEL in "${MODEL_ORDER[@]}"; do
    N_PC="${MODEL_PC[$MODEL]}"
    echo ""
    echo "  == Model: inbred_${MODEL} (PCs = ${N_PC}) =="

    # Determine qcovar argument
    if [ "${N_PC}" -eq 0 ]; then
        QCOVAR_ARG=""
    elif [ "${N_PC}" -eq 2 ]; then
        QCOVAR_ARG="--qcovar ${EIGENVEC_2PC}"
    else
        QCOVAR_ARG="--qcovar ${EIGENVEC_5PC}"
    fi

    for TRAIT in "${PHENOTYPES[@]}"; do
        PHENOTYPE_FILE="${GCTA_INPUTS_DIR}/${TRAIT}_for_gcta.txt"
        RESULT_PREFIX="${RESULTS}/result_${TRAIT}_h2_inbred_${MODEL}"
        LOG_FILE="${LOGS}/stage3_reml_${TRAIT}_${MODEL}.log"

        if [ ! -f "${RESULT_PREFIX}.hsq" ]; then
            if [ -f "${PHENOTYPE_FILE}" ]; then
                echo "    Running REML: ${TRAIT} | ${MODEL}..."
                "${GCTA}" \
                    --grm "${GRM_PREFIX}" \
                    --pheno "${PHENOTYPE_FILE}" \
                    ${QCOVAR_ARG} \
                    --reml \
                    --thread-num "${THREADS}" \
                    --out "${RESULT_PREFIX}" \
                    2>&1 | tee "${LOG_FILE}"
                echo "    Done: ${RESULT_PREFIX}.hsq"
            else
                echo "    WARNING: Phenotype file not found: ${PHENOTYPE_FILE}. Skipping."
            fi
        else
            echo "    Existing result found for ${TRAIT} | ${MODEL}. Skipping."
        fi
    done
done

echo ""
echo "--- STAGE 3 COMPLETE ---"
echo "====================================================="

# ====================================================================
# --- STAGE 4: Build Summary Table ---
#
# Parses every .hsq file in results/ and extracts:
#   Trait, Model, Heritability_h2 (V(G)/Vp), SE, P_value
#
# Writes a long-format tab-delimited summary file.
# ====================================================================
echo ""
echo "--- STAGE 4: Generating Heritability Summary ---"
SUMMARY="${RESULTS}/heritability_summary_V3_inbred_sensitivity.txt"

printf "Trait\tModel\tHeritability_h2\tSE\tP_value\n" > "${SUMMARY}"

for MODEL in "${MODEL_ORDER[@]}"; do
    for TRAIT in "${PHENOTYPES[@]}"; do
        HSQ_FILE="${RESULTS}/result_${TRAIT}_h2_inbred_${MODEL}.hsq"
        if [ -f "${HSQ_FILE}" ]; then
            H2=$(awk   '/^V\(G\)\/Vp\t/ {print $2}' "${HSQ_FILE}")
            SE=$(awk   '/^V\(G\)\/Vp\t/ {print $3}' "${HSQ_FILE}")
            PVAL=$(awk '/^Pval\t/        {print $2}' "${HSQ_FILE}")
            printf "%s\t%s\t%s\t%s\t%s\n" \
                "${TRAIT}" "inbred_${MODEL}" "${H2}" "${SE}" "${PVAL}" \
                >> "${SUMMARY}"
            echo "  ${TRAIT} | inbred_${MODEL}: h2=${H2}  SE=${SE}  P=${PVAL}"
        else
            echo "  WARNING: Missing .hsq for ${TRAIT} | ${MODEL}. Not included in summary."
        fi
    done
done

echo ""
echo "====================================================="
echo "### Sensitivity analysis complete! ###"
echo ""
echo "  Results  : ${RESULTS}"
echo "  Summary  : ${SUMMARY}"
echo "  Logs     : ${LOGS}"
echo "====================================================="
