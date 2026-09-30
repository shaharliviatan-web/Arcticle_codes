#!/usr/bin/env bash
# 03_annotate_gene.sh - what protein is 7HG0729030?
# Method mirrors step 05 check 1: DIAMOND BLASTP of the gene's protein against the
# SAME UniProt Swiss-Prot DIAMOND database the main annotation chain uses.
# RESULT: GDSL esterase/lipase family (SGNH hydrolase). All top hits are GDSL
# esterase/lipases at 43-45% identity over ~88% coverage -- a confident FAMILY
# assignment, NOT a specific ortholog. GDSL is a large family (100+ in Arabidopsis).
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"
OUT="$BRANCH/results/tables"; INT="$BRANCH/intermediates"
# extract the protein from the PGSB proteome (HC+LC; note this gene is HC)
awk -v g="$GENE_ID" '/^>/{k=($0 ~ g)} k' "$PGSB_PROT" > "$BRANCH/inputs/protein_7HG0729030.faa"
"$DIAMOND" blastp -q "$BRANCH/inputs/protein_7HG0729030.faa" -d "$SWISSPROT_DMND" \
  -o "$INT/swissprot_diamond_hits.tsv" --very-sensitive --max-target-seqs 5 \
  --outfmt 6 qseqid sseqid pident length qcovhsp scovhsp evalue bitscore stitle --quiet
{ printf "query\tsubject\tpident\tlength\tqcovhsp\tscovhsp\tevalue\tbitscore\tsubject_title\n"
  cat "$INT/swissprot_diamond_hits.tsv"; } > "$OUT/gene_annotation_swissprot.tsv"
echo "03 OK: $(wc -l < "$INT/swissprot_diamond_hits.tsv") Swiss-Prot hits -> results/tables/gene_annotation_swissprot.tsv"
head -3 "$INT/swissprot_diamond_hits.tsv" | cut -f3,5,7,9 | cut -c1-100
