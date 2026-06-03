#!/usr/bin/env bash
# 01_fetch_refs.sh
# Step 01 of 09_USED_canonical_betaglucan_gene_check.
# Fetch the reference protein sequence for each canonical (1,3;1,4)-beta-glucan
# gene from UniProt / NCBI, by accession (accessions were identified via UniProt
# REST + NCBI esearch by gene name + Hordeum vulgare; fetching by accession makes
# this deterministic and reproducible). Provenance recorded in logs/versions.txt.
#
# Output: inputs/canonical_refs.faa  (one representative protein per canonical gene,
#         header ">GENE accession source").
set -uo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/09_USED_canonical_betaglucan_gene_check
MAP="$BASE/inputs/canonical_refs_accessions.tsv"
OUT="$BASE/inputs/canonical_refs.faa"
PROV="$BASE/inputs/canonical_refs_provenance.tsv"
VERS="$BASE/logs/versions.txt"
mkdir -p "$BASE/logs"

[ -s "$MAP" ] || { echo "ERROR: accession map missing: $MAP" >&2; exit 1; }

echo "=== Step 01: fetch canonical reference proteins ===" | tee "$VERS"
{
  echo "date_run:    $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  echo "UniProt REST: https://rest.uniprot.org/uniprotkb/<ACC>.fasta"
  echo "NCBI efetch:  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi (db=protein,rettype=fasta)"
  echo ""
  echo "gene	accession	source	seq_length	header"
} >> "$VERS"

: > "$OUT"
printf "gene\taccession\tsource\tseq_length\torig_header\n" > "$PROV"

tmp="$TMPDIR/ref_one.$$"
n=0
# skip header line of the map
tail -n +2 "$MAP" | while IFS=$'\t' read -r gene acc source note; do
  [ -z "$gene" ] && continue
  if [ "$source" = "uniprot" ]; then
    url="https://rest.uniprot.org/uniprotkb/${acc}.fasta"
  else
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text"
  fi
  curl -sS --connect-timeout 30 --max-time 90 "$url" > "$tmp" 2>/dev/null
  # validate: must be FASTA with sequence
  if ! grep -q '^>' "$tmp" || [ "$(grep -v '^>' "$tmp" | tr -d '\n[:space:]' | wc -c)" -eq 0 ]; then
    echo "ERROR: failed to fetch $gene ($acc) from $source" >&2; exit 1
  fi
  orig_hdr=$(grep -m1 '^>' "$tmp" | sed 's/^>//' | cut -c1-90)
  seq=$(grep -v '^>' "$tmp" | tr -d '\n[:space:]' | tr -d '*')
  len=${#seq}
  # write normalized record: id = gene, then accession + source
  printf ">%s %s %s\n" "$gene" "$acc" "$source" >> "$OUT"
  printf "%s\n" "$seq" | fold -w 60 >> "$OUT"
  printf "%s\t%s\t%s\t%s\t%s\n" "$gene" "$acc" "$source" "$len" "$orig_hdr" >> "$PROV"
  printf "%s\t%s\t%s\t%s\t%s\n" "$gene" "$acc" "$source" "$len" "$orig_hdr" >> "$VERS"
  n=$((n+1))
done
rm -f "$tmp"

cnt=$(grep -c '^>' "$OUT")
echo ""
echo "sequences written: $cnt"
echo "provenance -> $PROV"
column -t -s $'\t' "$PROV" | cut -c1-150
[ "$cnt" -eq 10 ] || { echo "ERROR: expected 10 references, got $cnt" >&2; exit 1; }
echo "OK"
