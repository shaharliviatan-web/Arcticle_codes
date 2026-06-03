#!/usr/bin/env bash
# 03_reciprocal_blast.sh
# Step 03: reciprocal-best-hit support for the blastp assignment.
# Take the top forward HC hits per reference, pull their protein sequences, and
# blastp them BACK against the canonical references. A reference's assignment is
# reciprocal-confirmed if its best HC hit's best reverse hit is that same reference
# (resolves CslF paralogs). Also extracts the candidate proteins for the tree.
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/09_USED_canonical_betaglucan_gene_check
INT="$BASE/intermediates"; DB="$INT/blastdb"; LOGS="$BASE/logs"
FWD="$INT/blastp_refs_vs_hc.tsv"
HCAA=/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.aa.fa
BLASTP=/usr/bin/blastp
THREADS=8

[ -s "$FWD" ] || { echo "ERROR: forward blastp missing: $FWD" >&2; exit 1; }

echo "=== Step 03: reciprocal best hit ==="
# top-3 HC candidates per reference (by bitscore)
sort -t$'\t' -k1,1 -k7,7nr "$FWD" | awk -F'\t' '{c[$1]++; if(c[$1]<=3) print $2}' | sort -u > "$INT/candidate_ids.txt"
echo "candidate HC proteins (top-3/ref, unique): $(wc -l < "$INT/candidate_ids.txt")"

# extract candidate protein sequences from the HC proteome
awk 'NR==FNR{want[">"$1]=1; next} /^>/{p=($1 in want)} p' "$INT/candidate_ids.txt" "$HCAA" > "$INT/candidates.faa"
echo "extracted candidate sequences: $(grep -c '^>' "$INT/candidates.faa")"

# reverse blastp: candidates vs refs db
"$BLASTP" -query "$INT/candidates.faa" -db "$DB/refs" \
  -evalue 1e-10 -max_target_seqs 5 -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
  -out "$INT/blastp_hc_vs_refs.tsv" 2> "$LOGS/03_blastp_rev.log"
echo "reverse blastp rows: $(wc -l < "$INT/blastp_hc_vs_refs.tsv")"
echo "DONE"
