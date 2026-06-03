#!/usr/bin/env bash
# 02_makedb_and_blast.sh
# Step 02: locate the canonical references in Morex V3 by sequence, two ways.
#   (A) tblastn refs vs the Morex V3 genome  -> genomic coordinates (annotation-free)
#   (B) blastp  refs vs the Morex V3 HC proteome -> HORVU.MOREX.r3 id by best hit
# Also builds a protein db of the refs themselves (for the reciprocal-best-hit step).
# Idempotent: BLAST dbs are only built if missing.
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/09_USED_canonical_betaglucan_gene_check
DB="$BASE/intermediates/blastdb"; mkdir -p "$DB"
LOGS="$BASE/logs"
REFS="$BASE/inputs/canonical_refs.faa"

GENOME=/mnt/data/Barley_2021/morexV3/Barley_MorexV3_pseudomolecules.fasta
HCAA=/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.aa.fa
MAKEDB=/usr/bin/makeblastdb
TBLASTN=/usr/bin/tblastn
BLASTP=/usr/bin/blastp
THREADS=16

[ -s "$REFS" ] || { echo "ERROR: refs missing: $REFS" >&2; exit 1; }

echo "=== Step 02: build dbs + locate by sequence @ $(date -u '+%H:%M:%S UTC') ==="

# (1) genome nucleotide db
if [ ! -s "$DB/genome.nsq" ] && [ ! -s "$DB/genome.nin" ]; then
  echo "[makeblastdb] genome (nucl) ..."
  "$MAKEDB" -in "$GENOME" -dbtype nucl -out "$DB/genome" -title MorexV3 > "$LOGS/02_makedb_genome.log" 2>&1
fi
echo "genome db ready."

# (2) HC proteome protein db
if [ ! -s "$DB/hc_prot.psq" ] && [ ! -s "$DB/hc_prot.pin" ]; then
  echo "[makeblastdb] HC proteome (prot) ..."
  "$MAKEDB" -in "$HCAA" -dbtype prot -out "$DB/hc_prot" -title MorexV3_HC > "$LOGS/02_makedb_hcprot.log" 2>&1
fi
echo "HC proteome db ready."

# (3) refs protein db (for reciprocal best hit)
"$MAKEDB" -in "$REFS" -dbtype prot -out "$DB/refs" -title canonical_refs > "$LOGS/02_makedb_refs.log" 2>&1
echo "refs db ready."

# (A) tblastn refs vs genome
echo "[tblastn] refs vs genome ..."
"$TBLASTN" -query "$REFS" -db "$DB/genome" \
  -evalue 1e-10 -max_target_seqs 50 -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out "$BASE/intermediates/tblastn_refs_vs_genome.tsv" 2> "$LOGS/02_tblastn.log"
echo "  tblastn rows: $(wc -l < "$BASE/intermediates/tblastn_refs_vs_genome.tsv")"

# (B) blastp refs vs HC proteome (forward)
echo "[blastp] refs vs HC proteome (forward) ..."
"$BLASTP" -query "$REFS" -db "$DB/hc_prot" \
  -evalue 1e-10 -max_target_seqs 10 -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
  -out "$BASE/intermediates/blastp_refs_vs_hc.tsv" 2> "$LOGS/02_blastp_fwd.log"
echo "  blastp rows: $(wc -l < "$BASE/intermediates/blastp_refs_vs_hc.tsv")"

echo "DONE_MARKER $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
