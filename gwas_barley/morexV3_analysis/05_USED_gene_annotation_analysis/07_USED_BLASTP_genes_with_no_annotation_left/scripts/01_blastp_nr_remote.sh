#!/usr/bin/env bash
# 01_blastp_nr_remote.sh
# Step 01 of 07_USED_BLASTP_genes_with_no_annotation_left.
# Last-resort BLASTP of the 5 residual proteins against NCBI nr, restricted to
# Viridiplantae, run REMOTELY on NCBI servers (nr is not installed locally).
#
# Idempotent: if a .done marker exists the search is skipped. The search is
# retried a few times because -remote can hit transient NCBI errors.
set -uo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/07_USED_BLASTP_genes_with_no_annotation_left
BLASTP=/usr/bin/blastp
QUERY="$BASE/inputs/residual_5.faa"
OUT="$BASE/intermediates/residual_5_vs_nr.tsv"
DONE="$BASE/intermediates/residual_5_vs_nr.done"
LOG="$BASE/logs/01_blastp.log"
VERS="$BASE/logs/versions.txt"

[ -s "$QUERY" ] || { echo "ERROR: query missing: $QUERY" >&2; exit 1; }

# outfmt 6 fields: query, subject id, full title, %id, query coverage per subject,
# e-value, bitscore, subject taxids and scientific names (may be blank on -remote).
FMT="6 qseqid sseqid stitle pident qcovs evalue bitscore staxids sscinames"

echo "=== Step 01: remote BLASTP vs nr (Viridiplantae) ==="

{
  echo "tool:        $($BLASTP -version 2>&1 | head -1)"
  echo "database:    NCBI nr (remote)"
  echo "entrez_query: Viridiplantae[ORGN]"
  echo "evalue:      1e-5    max_target_seqs: 5"
  echo "accessed:    $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  echo "queries:     $(grep -c '^>' "$QUERY") proteins"
} > "$VERS"

if [ -f "$DONE" ] && [ -f "$OUT" ]; then
  echo "Already completed (found $DONE); skipping remote search."
  echo "rows: $(wc -l < "$OUT")  queries with hits: $(cut -f1 "$OUT" 2>/dev/null | sort -u | grep -c . || echo 0)/5"
  exit 0
fi

attempt=0; max=3; ok=0
while [ "$attempt" -lt "$max" ]; do
  attempt=$((attempt+1))
  echo "[attempt $attempt/$max] submitting remote blastp @ $(date -u '+%H:%M:%S UTC') ..."
  if "$BLASTP" \
        -query "$QUERY" \
        -db nr \
        -remote \
        -entrez_query "Viridiplantae[ORGN]" \
        -evalue 1e-5 \
        -max_target_seqs 5 \
        -outfmt "$FMT" \
        -out "$OUT" \
        2>> "$LOG"; then
    ok=1; break
  else
    echo "  attempt $attempt failed (see $LOG); waiting 30s before retry ..."
    sleep 30
  fi
done

if [ "$ok" -ne 1 ]; then
  echo "ERROR: remote blastp failed after $max attempts. See $LOG" >&2
  exit 1
fi

touch "$DONE"
echo "Done. raw output -> $OUT"
echo "rows: $(wc -l < "$OUT")"
echo "queries with >=1 hit: $(cut -f1 "$OUT" 2>/dev/null | sort -u | grep -c . || echo 0)/5"
echo "queries with no hit:  $(comm -23 <(grep '^>' "$QUERY" | sed 's/^>//' | sort -u) <(cut -f1 "$OUT" 2>/dev/null | sort -u) | wc -l)"
