#!/usr/bin/env Rscript
# 01_gene_position_and_distances.R
# Where 7HG0729030 sits relative to the two 7H signals, and how isolated it is.
# WHY THIS BRANCH EXISTS: fiber_L17 and starch_L05 lead SNPs are 154 bp apart --
# the only place in the study where two traits share one signal. Both loci have
# near-zero LD span (185 bp / 1,563 bp), so step 03 (FLANK_BP = 0) returned NO gene.
# This script asks what is actually there.
# Outputs -> results/tables/{gene_position.tsv, distances_to_signals.tsv, neighbouring_genes.tsv}
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("BRANCH"); stopifnot(nzchar(B))
GS <- as.numeric(Sys.getenv("GENE_START")); GE <- as.numeric(Sys.getenv("GENE_END"))
TSS <- as.numeric(Sys.getenv("GENE_TSS")); OUT <- file.path(B,"results","tables")

gp <- data.table(gene_id=Sys.getenv("GENE_ID"), chr=Sys.getenv("GENE_CHR"),
  gene_start=GS, gene_end=GE, strand=Sys.getenv("GENE_STRAND"), length_bp=GE-GS+1,
  TSS=TSS, tss_note="minus strand: TSS is the HIGH coordinate; promoter lies above it")
fwrite(gp, file.path(OUT,"gene_position.tsv"), sep="\t")

# distance to every 7H locus of the two traits
lo <- fread(Sys.getenv("LOCI_TABLE"))[chr=="7H"]
lo[, gap_to_gene_body := pmax(0, pmax(start-GE, GS-end))]
lo[, lead_to_gene_body := pmax(0, pmax(lead_bp-GE, GS-lead_bp))]
lo[, lead_offset_from_TSS := lead_bp - TSS]     # +ve = promoter side
setorder(lo, gap_to_gene_body)
fwrite(lo[,.(search_id,trait,lead_SNP,lead_bp,lead_neg_log10_p,start,end,span_kb,
             n_member_SNPs,gap_to_gene_body,lead_to_gene_body,lead_offset_from_TSS)],
       file.path(OUT,"distances_to_signals.tsv"), sep="\t")

# How isolated is the gene? Neighbours are taken with the SAME method and the SAME
# annotation as step 03: Ensembl Plants r62 GFF3, col3 == "gene", chromosomes 1H-7H,
# via `bedtools intersect` -- NOT a bespoke awk scan, so the branch stays on the
# project's existing gene-extraction method.
INT <- file.path(B,"intermediates")
gff <- file.path(INT,"genes_7H_ensembl.gff")
if (!file.exists(gff))
  system(sprintf("zcat %s | awk -F'\t' 'BEGIN{OFS=\"\t\"} $3==\"gene\" && $1 ~ /^[1-7]H$/' > %s",
                 Sys.getenv("ENSEMBL_GFF"), gff))
bed <- file.path(INT,"neighbour_window.bed")
writeLines(sprintf("%s\t%d\t%d\tneighbour_window", Sys.getenv("GENE_CHR"),
                   GS-400000-1L, GE+400000), bed)          # 1-based inclusive -> 0-based BED
raw <- file.path(INT,"neighbour_intersect_raw.tsv")
system(sprintf("%s intersect -a %s -b %s -wa -wb > %s", Sys.getenv("BEDTOOLS"), bed, gff, raw))
r <- fread(raw, header=FALSE)
g <- data.table(gene_id=sub("^gene:","",sub(";.*","",sub(".*ID=","",r$V13))),
                start=as.numeric(r$V8), end=as.numeric(r$V9), strand=r$V11)
g[, dist_to_our_gene := pmax(0, pmax(start-GE, GS-end))]
g[, dist_to_fiber_lead := pmax(0, pmax(start-as.numeric(Sys.getenv("FIBER_LEAD")),
                                        as.numeric(Sys.getenv("FIBER_LEAD"))-end))]
setorder(g, dist_to_fiber_lead)
fwrite(g, file.path(OUT,"neighbouring_genes.tsv"), sep="\t")
cat("01 OK: gene_position / distances_to_signals / neighbouring_genes\n")
