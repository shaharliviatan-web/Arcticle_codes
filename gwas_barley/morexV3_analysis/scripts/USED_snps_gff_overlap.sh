#!/bin/bash
# THIS SCRIPT IS ACTIVATED BY THE SCRIPT: run_evgeny_analysis.sh

# USAGE: snps_gff_overlap.sh SNPs.list N genes.gff feature_type
# EXAMPLE: snps_gff_overlap.sh SNPs.list 1000 annotations.gff gene
#
# VERSION 1.2: All comments and messages translated to English.
# Now handles compressed (.gz) GFF files automatically.
# DEBUGGING: The final cleanup step is disabled.

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 SNPs.list N genes.gff feature_type" >&2
    exit 1
fi

SNP_LIST=$1      # list of SNPs in chr1H:6409022 format
N=$2             # window size (e.g., 250000)
GFF=$3           # GFF file (can be a .gz file)
FEATURE=$4       # feature type (e.g., gene, mRNA, exon)

# Create temporary files and print their names for debugging
TMP_BED=$(mktemp)
TMP_GFF=$(mktemp)
echo "--- DEBUGGING MODE ---" >&2
echo "Temp BED file is at: $TMP_BED" >&2
echo "Temp GFF file is at: $TMP_GFF" >&2
echo "----------------------" >&2

# --- Step 1: Convert SNP list into BED file with a +/- N window ---
awk -v N="$N" '
BEGIN {OFS="\t"}
{
  if ($0 ~ /:/ && $0 !~ /^#/) {
    split($1, a, ":");
    chr=a[1]; pos=a[2];
    start=pos - N; if (start < 0) start = 0;
    end=pos + N;
    snpid=$1;
    print chr, start, end, snpid;
  }
}' "$SNP_LIST" > "$TMP_BED"

# --- Step 2: Filter the GFF file by feature type (supports .gz files) ---
# This awk command correctly checks ONLY the 3rd column for the $FEATURE
if [[ "$GFF" == *.gz ]]; then
    gunzip -c "$GFF" | awk -v F="$FEATURE" 'BEGIN{OFS="\t"} $3 == F' > "$TMP_GFF"
else
    awk -v F="$FEATURE" 'BEGIN{OFS="\t"} $3 == F' "$GFF" > "$TMP_GFF"
fi

# --- Step 3: Find overlaps using bedtools ---
bedtools intersect -a "$TMP_BED" -b "$TMP_GFF" -wa -wb \
| awk 'BEGIN{OFS="\t"} {snpid=$4; range=$1":"$2"-"$3; print snpid, range, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}'

# --- DEBUG: Cleanup is disabled to allow inspection of temp files ---
# rm -f "$TMP_BED" "$TMP_GFF"
echo "--- DEBUGGING: Cleanup step was SKIPPED. ---" >&2