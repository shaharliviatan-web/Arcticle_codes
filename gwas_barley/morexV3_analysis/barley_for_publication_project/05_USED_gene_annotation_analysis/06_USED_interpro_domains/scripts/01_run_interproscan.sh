#!/usr/bin/env bash
# 01_run_interproscan.sh
# Step 01 of 06_USED_interpro_domains.
# Submit each of the 45 FDR proteins to EBI InterProScan 5 via the official
# iprscan5 REST client. One sequence per job, in a loop, with polite polling.
# Idempotent / resumable: a sequence whose JSON result already exists is skipped.
#
# Per-sequence results land in intermediates/raw/<seqid>.{tsv.tsv,json.json,gff.txt}.
# GO terms requested (--goterms) and pathways (--pathways); default application set
# (includes Pfam) is used (no --appl given).
set -uo pipefail   # NOTE: no -e; one failed job must not abort the whole batch.

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/05_USED_gene_annotation_analysis/06_USED_interpro_domains
SCRIPTS="$BASE/scripts"
PYLIB="$SCRIPTS/pylib"
SEQS="$BASE/intermediates/seqs"
RAW="$BASE/intermediates/raw"
LOGS="$BASE/logs"
CLIENT="$SCRIPTS/iprscan5.py"
mkdir -p "$RAW" "$LOGS"

EMAIL="shaharliviatan@gmail.com"   # required by EBI; identifies the submitter
POLLFREQ=5                          # client status poll interval (s)
SLEEP_BETWEEN=3                     # polite pause between sequence submissions (s)

export PYTHONPATH="$PYLIB"

PROGRESS="$LOGS/01_progress.log"
echo "=== Step 01: InterProScan run @ $(date -u '+%Y-%m-%d %H:%M:%S UTC') ===" | tee -a "$PROGRESS"

# ---------------------------------------------------------------------------
# Split the input FASTA into one file per sequence (2026-09-09).
#
# v1 did this by hand and left the result in intermediates/seqs/, so the script
# silently kept submitting the OLD gene set after the input FASTA was replaced.
# Regenerating here makes the step self-contained and prevents that class of bug.
# The directory is rebuilt from scratch every run; intermediates/raw/ is NOT
# touched, so completed jobs are still skipped and the run stays resumable.
# ---------------------------------------------------------------------------
IN_FAA="$BASE/inputs/fdr_proteins.faa"
[ -s "$IN_FAA" ] || { echo "ERROR: missing input FASTA: $IN_FAA" >&2; exit 1; }
rm -rf "$SEQS"; mkdir -p "$SEQS"
awk -v out="$SEQS" '
  /^>/ { id=substr($1,2); gsub(/[^A-Za-z0-9._-]/,"_",id); f=out "/" id ".fasta" }
  { if (f != "") print > f }
' "$IN_FAA"
n_split=$(ls -1 "$SEQS"/*.fasta 2>/dev/null | wc -l)
n_faa=$(grep -c '^>' "$IN_FAA")
if [ "$n_split" -ne "$n_faa" ]; then
  echo "ERROR: split produced $n_split files from $n_faa sequences." >&2; exit 1
fi
echo "Split $n_faa sequences into $SEQS" | tee -a "$PROGRESS"
N_TOTAL=$n_faa

total=0; done_already=0; ran=0; failed=0
fails=""

for fa in "$SEQS"/*.fasta; do
  total=$((total+1))
  seqid=$(basename "$fa" .fasta)
  json="$RAW/$seqid.json.json"

  if [ -s "$json" ] && python3 -c "import json,sys; json.load(open('$json'))" 2>/dev/null; then
    done_already=$((done_already+1))
    continue
  fi

  echo "[$total/$N_TOTAL] submitting $seqid ..." | tee -a "$PROGRESS"
  python3 "$CLIENT" \
      --email "$EMAIL" \
      --stype p \
      --goterms \
      --pathways \
      --sequence "$fa" \
      --outfile "$RAW/$seqid" \
      --outformat tsv,json,gff \
      --pollFreq "$POLLFREQ" \
      > "$LOGS/${seqid}.iprscan.log" 2>&1

  if [ -s "$json" ] && python3 -c "import json,sys; json.load(open('$json'))" 2>/dev/null; then
    ran=$((ran+1))
    echo "      OK $seqid" | tee -a "$PROGRESS"
  else
    failed=$((failed+1)); fails="$fails $seqid"
    echo "      FAILED $seqid (see logs/${seqid}.iprscan.log)" | tee -a "$PROGRESS"
  fi
  sleep "$SLEEP_BETWEEN"
done

echo "--- summary: total=$total already_done=$done_already ran_now=$ran failed=$failed ---" | tee -a "$PROGRESS"
[ -n "$fails" ] && echo "FAILED:$fails" | tee -a "$PROGRESS"

# Record the InterProScan engine version from a retrieved JSON.
J=$(ls "$RAW"/*.json.json 2>/dev/null | head -1 || true)
if [ -n "${J:-}" ]; then
  iprv=$(python3 -c "import json;print(json.load(open('$J')).get('interproscan-version'))" 2>/dev/null)
  echo "interproscan-version (from results): $iprv" | tee -a "$PROGRESS"
  grep -q '^interproscan-version:' "$LOGS/versions.txt" 2>/dev/null \
    && sed -i "s|^interproscan-version:.*|interproscan-version:    $iprv|" "$LOGS/versions.txt" || true
fi

echo "DONE_MARKER $(date -u '+%Y-%m-%d %H:%M:%S UTC')" | tee -a "$PROGRESS"
