#!/usr/bin/env bash
# 05_gff3_crosscheck.sh
# Step 05, cross-check (a): independently of BLAST, grep the Morex V3 HC GFF3
# functional description= fields for cellulose-synthase / glucan-synthase /
# glucanase terms and list the HORVU.MOREX.r3 IDs found.
#
# CAVEAT: only ~6-7% of HC genes carry a functional description= in this GFF3, so
# this grep can CORROBORATE a BLAST assignment but its SILENCE means nothing
# (absence of a description != absence of the gene). BLAST is the primary method.
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/09_USED_canonical_betaglucan_gene_check
GFF=/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.gff3
OUT="$BASE/intermediates/gff3_description_hits.tsv"

echo "=== Step 05: GFF3 description cross-check (caveat: ~6-7% genes annotated) ==="
# annotation coverage for context
tot=$(awk -F'\t' '$3=="mRNA"' "$GFF" | wc -l)
desc=$(awk -F'\t' '$3=="mRNA" && /description=/ && !/description=NA/' "$GFF" | wc -l)
echo "HC mRNAs: $tot ; with a non-NA description: $desc ($(awk "BEGIN{printf \"%.1f\", 100*$desc/$tot}")%)"

printf "gene_id\tchr\tstart\tend\tstrand\tdescription\n" > "$OUT"
awk -F'\t' '$3=="mRNA"' "$GFF" | grep -iE 'description=[^;]*(cellulose synthase|glucan synthase|glucan-synthase|glucanase|1,3;1,4|1,3-1,4|lichenase|beta-glucan)' \
 | awk -F'\t' '{
     attr=$9; id=attr; sub(/.*ID=/,"",id); sub(/;.*/,"",id); sub(/\.[0-9]+$/,"",id);
     d=attr; sub(/.*description=/,"",d); sub(/;.*/,"",d);
     print id"\t"$1"\t"$4"\t"$5"\t"$7"\t"d}' | sort -u >> "$OUT"

echo "matching HORVU genes: $(($(wc -l < "$OUT")-1))"
column -t -s $'\t' "$OUT" | cut -c1-140
echo "-> $OUT"
echo "DONE"
