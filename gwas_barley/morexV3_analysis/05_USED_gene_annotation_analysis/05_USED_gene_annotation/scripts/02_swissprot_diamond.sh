#!/usr/bin/env bash
# 02_swissprot_diamond.sh
# Step 02 of 05_USED_gene_annotation.
# Download UniProt Swiss-Prot, build a DIAMOND db, and BLASTP the FDR proteins.
#
# Evidence per hit: qseqid sseqid stitle pident qcovhsp length evalue bitscore
# (stitle carries the Swiss-Prot protein name + OS= organism + OX= taxid).
set -euo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/05_USED_gene_annotation
INTER="$BASE/intermediates"
INPUTS="$BASE/inputs"
LOGS="$BASE/logs"
mkdir -p "$INTER" "$LOGS"

DIAMOND=/usr/bin/diamond
THREADS=32                       # well under the 100-core cap; job is tiny anyway

QUERY="$INPUTS/fdr_proteins.faa"
SPROT_GZ="$INTER/uniprot_sprot.fasta.gz"
SPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
RELDATE_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt"
DB="$INTER/swissprot"            # DIAMOND appends .dmnd
OUT="$INTER/fdr_vs_swissprot.diamond.tsv"
VERS="$INTER/versions.txt"

[ -s "$QUERY" ] || { echo "ERROR: query FASTA missing: $QUERY" >&2; exit 1; }

echo "=== Step 02: Swiss-Prot DIAMOND BLASTP ==="

# --- Record versions / release ---
{
  echo "DIAMOND version:"
  "$DIAMOND" version
  echo ""
  echo "Swiss-Prot release (UniProt reldate.txt):"
} > "$VERS"
curl -sS --connect-timeout 30 "$RELDATE_URL" >> "$VERS" 2>&1 || echo "(reldate fetch failed)" >> "$VERS"
echo "Recorded versions ->" "$VERS"

# --- Download Swiss-Prot (idempotent: skip if size matches remote) ---
REMOTE_LEN=$(curl -sSI --connect-timeout 30 "$SPROT_URL" | awk 'BEGIN{IGNORECASE=1}/content-length/{print $2}' | tr -d '\r')
if [ -s "$SPROT_GZ" ] && [ "$(stat -c%s "$SPROT_GZ")" = "$REMOTE_LEN" ]; then
  echo "Swiss-Prot already present and size matches ($REMOTE_LEN bytes); skipping download."
else
  echo "Downloading Swiss-Prot ($REMOTE_LEN bytes) ..."
  curl -sS --connect-timeout 30 -o "$SPROT_GZ" "$SPROT_URL"
  echo "Downloaded $(stat -c%s "$SPROT_GZ") bytes."
fi
gzip -t "$SPROT_GZ"
echo "gzip integrity OK."
echo "Swiss-Prot sequences: $(zcat "$SPROT_GZ" | grep -c '^>')"

# --- Build DIAMOND database ---
echo "Building DIAMOND db ..."
"$DIAMOND" makedb --in "$SPROT_GZ" -d "$DB" --threads "$THREADS" 2> "$LOGS/02_makedb.log"
echo "DB built -> ${DB}.dmnd ($(stat -c%s "${DB}.dmnd") bytes)"

# --- BLASTP ---
echo "Running DIAMOND blastp (--very-sensitive -e 1e-5 --max-target-seqs 5) ..."
"$DIAMOND" blastp \
  --very-sensitive \
  -q "$QUERY" \
  -d "$DB" \
  -e 1e-5 \
  --max-target-seqs 5 \
  -p "$THREADS" \
  --outfmt 6 qseqid sseqid stitle pident qcovhsp length evalue bitscore \
  -o "$OUT" \
  2> "$LOGS/02_blastp.log"

echo ""
echo "Raw blast output -> $OUT"
echo "Total HSP rows:       $(wc -l < "$OUT")"
echo "Queries with >=1 hit: $(cut -f1 "$OUT" | sort -u | wc -l) of 45"
echo "Done."
