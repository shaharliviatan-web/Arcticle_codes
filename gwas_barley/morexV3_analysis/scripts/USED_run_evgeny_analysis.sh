#!/bin/bash
# THIS SCRIPT ACTIVATING THE SCRIPT: snps_gff_overlap.sh
# ===================================================================
# Wrapper Script for Running Evgeny's snps_gff_overlap.sh
#
# Author: Shahar & Gemini
# Version: 1.1 - All comments and messages are in English.
#
# This script acts as a user-friendly wrapper. It takes the same
# arguments but automatically creates an output directory and a
# descriptive output filename.
# ===================================================================

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 SNPs.list N genes.gff feature_type" >&2
    echo "This script will save the output to the 'evgeny_script_outputs' directory." >&2
    exit 1
fi

SNP_LIST=$1
WINDOW_SIZE=$2
GFF_FILE=$3
FEATURE=$4

# --- Directory and File Setup ---
SCRIPT_DIR=$(dirname $(readlink -f "$0"))
ANALYSIS_DIR=$(dirname "$SCRIPT_DIR")
OUTPUT_DIR="${ANALYSIS_DIR}/USED_fixed_window_search_results"
mkdir -p "$OUTPUT_DIR"

# Create a descriptive output filename
INPUT_BASENAME=$(basename "$SNP_LIST" .txt)
OUTPUT_FILENAME="${INPUT_BASENAME}_${FEATURE}_${WINDOW_SIZE}bp.txt"
FULL_OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_FILENAME}"

echo "--- Running Evgeny's SNP Overlap Analysis ---"
echo "SNP List:    $(basename $SNP_LIST)"
echo "Window Size: ${WINDOW_SIZE} bp"
echo "Feature:     ${FEATURE}"
echo "-------------------------------------------------"
echo "Saving results to: ${FULL_OUTPUT_PATH}"

# --- Execute the core script and redirect output ---
"$SCRIPT_DIR/USED_snps_gff_overlap.sh" "$SNP_LIST" "$WINDOW_SIZE" "$GFF_FILE" "$FEATURE" > "$FULL_OUTPUT_PATH"
echo "--- Analysis Complete! ---"

