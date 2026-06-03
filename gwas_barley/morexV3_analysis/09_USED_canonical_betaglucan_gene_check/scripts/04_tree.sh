#!/usr/bin/env bash
# 04_tree.sh
# Step 04: confirm CslF/CslH family assignments are not paralog swaps, by building
# a quick MSA (mafft) + neighbour-joining tree (R/ape) of the CslF/CslH references
# together with their candidate Morex V3 hits. The two beta-glucanases (Glb1/Glb2)
# are a different protein family and are excluded from this tree.
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/09_USED_canonical_betaglucan_gene_check
INT="$BASE/intermediates"; RES="$BASE/results"; LOGS="$BASE/logs"
REFS="$BASE/inputs/canonical_refs.faa"
FWD="$INT/blastp_refs_vs_hc.tsv"
MAFFT=/usr/bin/mafft
mkdir -p "$RES/tables"

echo "=== Step 04: CslF/CslH alignment + NJ tree ==="

# Csl references (exclude the glucanases)
awk '/^>/{p=($0 ~ /^>HvCsl/)} p' "$REFS" > "$INT/tree_input.faa"

# Csl candidate V3 proteins = forward HC hits of any Csl reference
awk -F'\t' '$1 ~ /^HvCsl/ {print $2}' "$FWD" | sort -u > "$INT/csl_candidate_ids.txt"
HCAA=/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.aa.fa
awk 'NR==FNR{want[">"$1]=1; next} /^>/{p=($1 in want)} p' "$INT/csl_candidate_ids.txt" "$HCAA" >> "$INT/tree_input.faa"
echo "tree taxa (refs + V3 candidates): $(grep -c '^>' "$INT/tree_input.faa")"

# MSA
"$MAFFT" --auto --anysymbol --quiet "$INT/tree_input.faa" > "$INT/tree_align.fasta" 2> "$LOGS/04_mafft.log"
echo "alignment columns: $(awk '/^>/{next}{print length($0); exit}' "$INT/tree_align.fasta")"

# NJ tree in R (seqinr identity distance + ape::nj); report each V3 tip's nearest reference
/usr/bin/Rscript - "$INT/tree_align.fasta" "$RES/tables/csl_family_tree.nwk" "$RES/tables/csl_tree_nearest_ref.tsv" <<'RS' 2> "$LOGS/04_tree_R.log"
args <- commandArgs(trailingOnly=TRUE)
suppressMessages({library(seqinr); library(ape)})
aln <- read.alignment(args[1], format="fasta")
d <- dist.alignment(aln, matrix="identity")      # sqrt(1-%id); smaller = more similar
tr <- nj(as.matrix(d))
write.tree(tr, args[2])
# for each V3 (HORVU) tip, nearest reference tip by patristic distance
M <- as.matrix(d); labs <- rownames(M)
refs <- grep("^HvCsl", labs, value=TRUE); v3 <- grep("^HORVU", labs, value=TRUE)
out <- data.frame(v3_gene=character(), nearest_ref=character(), identity_dist=numeric())
for (g in v3) {
  dd <- M[g, refs]; nr <- names(which.min(dd))
  out <- rbind(out, data.frame(v3_gene=g, nearest_ref=nr, identity_dist=round(min(dd),4)))
}
out <- out[order(out$nearest_ref, out$identity_dist),]
write.table(out, args[3], sep="\t", quote=FALSE, row.names=FALSE)
cat("V3 candidate -> nearest CslF/CslH reference (by alignment identity distance):\n")
print(out, row.names=FALSE)
RS
echo "tree -> $RES/tables/csl_family_tree.nwk"
echo "nearest-ref table -> $RES/tables/csl_tree_nearest_ref.tsv"
echo "DONE"
