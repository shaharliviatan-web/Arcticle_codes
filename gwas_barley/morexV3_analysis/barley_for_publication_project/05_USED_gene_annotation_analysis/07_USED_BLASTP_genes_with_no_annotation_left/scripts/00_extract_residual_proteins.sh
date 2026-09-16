#!/bin/bash
# =============================================================================
# 00_extract_residual_proteins.sh
#
# WHAT   Extract the proteins that still have NO functional call after Swiss-Prot
#        (step 05) and InterProScan (step 06), for a last-resort NCBI nr BLASTP.
#
# REWRITTEN 2026-09-09. v1 hard-coded a list of 5 gene IDs, so it aborted as soon
# as the gene set changed. The residual set is now DERIVED:
#
#   residual = {Swiss-Prot WEAK_OR_NO_HIT}  INTERSECT  {no InterPro match}
#
#   A gene rescued by either line of evidence is not residual. Taking the
#   intersection (not the union) is the point of running two sources.
#
# READS  05_.../results/tables/residual_weak_or_no_hit.txt
#        06_.../results/tables/fdr_gene_interpro.tsv
#        05_.../inputs/fdr_proteins.faa
# WRITES inputs/residual_ids.txt, inputs/residual.faa
#        (v1 named these residual_5*; the count is no longer fixed)
# =============================================================================
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/05_USED_gene_annotation_analysis
BASE="$ROOT/07_USED_BLASTP_genes_with_no_annotation_left"
SRC="$ROOT/05_USED_gene_annotation/inputs/fdr_proteins.faa"
SP="$ROOT/05_USED_gene_annotation/results/tables/residual_weak_or_no_hit.txt"
IP="$ROOT/06_USED_interpro_domains/results/tables/fdr_gene_interpro.tsv"
OUT="$BASE/inputs/residual.faa"
IDS="$BASE/inputs/residual_ids.txt"
mkdir -p "$BASE/inputs"

for f in "$SRC" "$SP" "$IP"; do
  [ -s "$f" ] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

echo "=== Step 00: derive and extract residual proteins ==="
# genes with no InterPro match (status column != INTERPRO_MATCH)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="gene_id")g=i; if($i=="status")s=i}; next} $s!="INTERPRO_MATCH"{print $g}' "$IP" \
  | sort -u > "$TMPDIR/no_interpro.$$"
sort -u "$SP" > "$TMPDIR/sp_weak.$$"
comm -12 "$TMPDIR/sp_weak.$$" "$TMPDIR/no_interpro.$$" > "$IDS"
rm -f "$TMPDIR/no_interpro.$$" "$TMPDIR/sp_weak.$$"

n_ids=$(wc -l < "$IDS")
echo "Swiss-Prot weak/no hit : $(wc -l < "$SP")"
echo "no InterPro match      : (see $IP)"
echo "residual (neither)     : $n_ids"

if [ "$n_ids" -eq 0 ]; then
  : > "$OUT"
  echo "Nothing residual - every gene has a call from Swiss-Prot or InterPro."
  echo "Steps 01/02 of this directory have nothing to do; that is a valid outcome."
  exit 0
fi

: > "$OUT"
while read -r g; do
  [ -n "$g" ] || continue
  awk -v g="$g" '/^>/ { keep = ($0 ~ "^>"g"\\.") } keep { print }' "$SRC" >> "$OUT"
done < "$IDS"

n=$(grep -c '^>' "$OUT")
echo "sequences written      : $n"
grep '^>' "$OUT" | sed 's/^/  /'
echo "stray '*' stop codons  : $(grep -c '\*' "$OUT" || true)"
[ "$n" -eq "$n_ids" ] || { echo "ERROR: expected $n_ids sequences, got $n" >&2; exit 1; }
