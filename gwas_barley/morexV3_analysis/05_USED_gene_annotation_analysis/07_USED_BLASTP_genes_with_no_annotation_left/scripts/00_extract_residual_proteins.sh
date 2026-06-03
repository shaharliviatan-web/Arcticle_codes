#!/usr/bin/env bash
# 00_extract_residual_proteins.sh
# Step 00 of 07_USED_BLASTP_genes_with_no_annotation_left.
# Pull the 5 still-unannotated residual proteins from the step-05 FASTA into
# inputs/residual_5.faa. These are the genes left after Swiss-Prot (05) and
# InterProScan (06):
#   1HG0079220, 1HG0079290, 2HG0112840, 4HG0339140 (member-db signatures only)
#   2HG0112890 (no match, 78 aa)
set -euo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/07_USED_BLASTP_genes_with_no_annotation_left
SRC=/mnt/data/shahar/gwas_barley/morexV3_analysis/05_USED_gene_annotation/inputs/fdr_proteins.faa
OUT="$BASE/inputs/residual_5.faa"
IDS="$BASE/inputs/residual_5_ids.txt"

GENES="HORVU.MOREX.r3.1HG0079220 HORVU.MOREX.r3.1HG0079290 HORVU.MOREX.r3.2HG0112840 HORVU.MOREX.r3.4HG0339140 HORVU.MOREX.r3.2HG0112890"

[ -s "$SRC" ] || { echo "ERROR: source FASTA missing: $SRC" >&2; exit 1; }
printf "%s\n" $GENES > "$IDS"

# Extract each gene's record (headers are <gene>.<isoform>) into a single FASTA.
: > "$OUT"
for g in $GENES; do
  awk -v g="$g" '
    /^>/ { keep = ($0 ~ "^>"g"\\.") }
    keep { print }
  ' "$SRC" >> "$OUT"
done

n=$(grep -c '^>' "$OUT")
echo "=== Step 00: extract residual proteins ==="
echo "source: $SRC"
echo "requested genes: 5"
echo "sequences written: $n"
grep '^>' "$OUT" | sed 's/^/  /'
echo "any stray '*' (stop codon) in seqs: $(grep -c '\*' "$OUT")"
[ "$n" -eq 5 ] || { echo "ERROR: expected 5 sequences, got $n" >&2; exit 1; }
echo "Wrote: $OUT"
