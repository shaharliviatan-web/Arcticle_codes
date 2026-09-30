#!/usr/bin/env bash
# 03b_interpro_and_nr.sh - annotation checks 2 and 3, mirroring the main chain
# (05_USED_gene_annotation_analysis steps 06 and 07).
#
#   check 2 = EBI InterProScan 5 (remote) : domains from ALL member databases + GO + pathways
#   check 3 = NCBI nr BLASTP     (remote) : independent homology evidence
#
# NOTE: the PGSB proteome carries a trailing '*' stop character. InterProScan REJECTS
# sequences containing '*', so the stripped copy protein_7HG0729030_nostop.faa is used.
# Outputs -> results/tables/{gene_annotation_interpro.tsv, gene_annotation_nr.tsv}
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"
IPDIR="$PROJECT_ROOT/05_USED_gene_annotation_analysis/06_USED_interpro_domains/scripts"
export PYTHONPATH="$IPDIR/pylib"
mkdir -p "$BRANCH/intermediates/interpro" "$BRANCH/intermediates/nr"
sed 's/\*//g' "$BRANCH/inputs/protein_7HG0729030.faa" > "$BRANCH/inputs/protein_7HG0729030_nostop.faa"

# ---- check 2: InterProScan (~3-5 min for one sequence) ----
if [ ! -s "$BRANCH/intermediates/interpro/g30.tsv.tsv" ]; then
  ( cd "$BRANCH/intermediates/interpro" && python3 "$IPDIR/iprscan5.py" \
      --email shaharliviatan@gmail.com --stype p --goterms --pathways --pollFreq 10 \
      --outfile g30 --outformat tsv "$BRANCH/inputs/protein_7HG0729030_nostop.faa" ) \
    > "$BRANCH/logs/interpro.log" 2>&1
fi
{ printf "member_db\tsignature_acc\tsignature_desc\tstart\tend\tevalue\tinterpro_acc\tinterpro_desc\tgo_terms\n"
  awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$4,$5,$6,$7,$8,$9,$12,$13,$14}' \
      "$BRANCH/intermediates/interpro/g30.tsv.tsv"; } > "$BRANCH/results/tables/gene_annotation_interpro.tsv"

# ---- check 3: NCBI nr BLASTP (remote; minutes) ----
if [ ! -s "$BRANCH/intermediates/nr/g30_vs_nr.tsv" ]; then
  /usr/bin/blastp -query "$BRANCH/inputs/protein_7HG0729030_nostop.faa" -db nr -remote \
    -evalue 1e-10 -max_target_seqs 25 \
    -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore stitle' \
    -out "$BRANCH/intermediates/nr/g30_vs_nr.tsv" > "$BRANCH/logs/nr_blast.log" 2>&1 || true
fi
if [ -s "$BRANCH/intermediates/nr/g30_vs_nr.tsv" ]; then
  { printf "query\tsubject\tpident\tlength\tqcovs\tevalue\tbitscore\tsubject_title\n"
    cat "$BRANCH/intermediates/nr/g30_vs_nr.tsv"; } > "$BRANCH/results/tables/gene_annotation_nr.tsv"
fi
echo "03b OK"
