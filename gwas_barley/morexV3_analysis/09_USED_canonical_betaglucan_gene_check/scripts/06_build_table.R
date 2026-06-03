#!/usr/bin/env Rscript
# 06_build_table.R
# Step 06: assign each canonical gene to its Morex V3 ortholog (identity-based
# reciprocal best hit), pull genomic coordinates, cross-check tblastn vs blastp,
# and measure distance to every beta-glucan / fiber / starch lead SNP.
# Deliverable: results/tables/canonical_betaglucan_gene_distances.tsv
#              one row per (canonical gene x trait).
#
# Paralog rule: references and Morex V3 are both Hordeum vulgare, so true
# orthologs align at ~99-100% identity and paralogs at <85%. The ortholog is
# therefore the forward hit with the HIGHEST % identity (NOT highest bitscore -
# bitscore can favour a full-length lower-identity paralog over a partial-model
# true ortholog, e.g. CslF4). Reciprocal best hit + the NJ tree confirm this.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
options(stringsAsFactors = FALSE)

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis"
B09  <- file.path(ROOT, "09_USED_canonical_betaglucan_gene_check")
INT  <- file.path(B09, "intermediates")
GFF  <- "/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.gff3"
LEAD <- file.path(ROOT, "03_USED_candidate_genes_around_leading_snps/results/tables/lead_loci.tsv")
WIN  <- file.path(ROOT, "04_USED_haplotype_analysis_crosshap/00_config/gene_windows.tsv")
HAP  <- file.path(ROOT, "04_USED_haplotype_analysis_crosshap/05_results/tables/gene_annotation_review.csv")
ACC  <- file.path(B09, "inputs/canonical_refs_accessions.tsv")
OUT  <- file.path(B09, "results/tables/canonical_betaglucan_gene_distances.tsv")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("betaglucan", "fiber", "starch")
norm_chr <- function(x) sub("^chr", "", x)
gid <- function(x) sub("\\.[0-9]+$", "", x)

# ---- forward / reverse blastp ----
bn <- c("qseqid","sseqid","pident","length","qcovs","evalue","bitscore")
fwd <- read.delim(file.path(INT,"blastp_refs_vs_hc.tsv"), header=FALSE, col.names=bn)
rev <- read.delim(file.path(INT,"blastp_hc_vs_refs.tsv"), header=FALSE, col.names=bn)
fwd$gene <- gid(fwd$sseqid); rev$qgene <- gid(rev$qseqid)

acc <- read.delim(ACC)   # gene, accession, source, note

# reverse best ref per HC gene, by identity then bitscore
rev_ord <- rev[order(rev$qgene, -rev$pident, -rev$bitscore), ]
rev_best <- rev_ord[!duplicated(rev_ord$qgene), c("qgene","sseqid","pident")]
names(rev_best) <- c("hc_gene","rev_best_ref","rev_best_pid")

assign_rows <- list()
for (g in acc$gene) {
  h <- fwd[fwd$qseqid == g, ]
  if (nrow(h) == 0) { next }
  h <- h[order(-h$pident, -h$bitscore), ]
  best <- h[1, ]                       # highest-identity forward hit = ortholog
  hc <- best$gene
  rb <- rev_best[rev_best$hc_gene == hc, ]
  rbh <- if (nrow(rb) && rb$rev_best_ref == g) "Y" else "N"
  assign_rows[[g]] <- data.frame(
    canonical_gene = g, horvu_id = hc,
    pct_identity = round(best$pident,2), pct_coverage = best$qcovs,
    evalue = best$evalue, reciprocal_best_hit = rbh)
}
asn <- do.call(rbind, assign_rows)
asn <- merge(acc[,c("gene","accession")], asn, by.x="gene", by.y="canonical_gene")
names(asn)[1:2] <- c("canonical_gene","source_accession")

# ---- GFF3 gene coordinates for the assigned HC genes ----
gff <- read.delim(GFF, header=FALSE, comment.char="#", quote="",
                  col.names=c("chr","src","feat","start","end","score","strand","frame","attr"))
gg <- gff[gff$feat=="gene", ]
gg$id <- sub(";.*","", sub(".*ID=","", gg$attr))
gg$chrn <- norm_chr(gg$chr)
gmap <- gg[match(asn$horvu_id, gg$id), c("chrn","start","end","strand")]
names(gmap) <- c("chr","gene_start","gene_end","strand")
asn <- cbind(asn, gmap)

# ---- tblastn genomic cross-check (if available) ----
tb_file <- file.path(INT,"tblastn_refs_vs_genome.tsv")
asn$tblastn_horvu <- NA_character_; asn$tblastn_chr <- NA_character_
asn$tblastn_start <- NA_integer_; asn$tblastn_end <- NA_integer_
if (file.exists(tb_file) && file.info(tb_file)$size > 0) {
  tn <- c("qseqid","sseqid","pident","length","qcovs","mismatch","gapopen",
          "qstart","qend","sstart","send","evalue","bitscore")
  tb <- read.delim(tb_file, header=FALSE, col.names=tn)
  tb$chrn <- norm_chr(tb$sseqid)
  tb$lo <- pmin(tb$sstart, tb$send); tb$hi <- pmax(tb$sstart, tb$send)
  for (i in seq_len(nrow(asn))) {
    g <- asn$canonical_gene[i]
    hh <- tb[tb$qseqid==g, ]
    if (!nrow(hh)) next
    # anchor on the highest-IDENTITY HSP (length-filtered), NOT highest bitscore -
    # bitscore can favour a full-length paralog elsewhere over the true ortholog
    # (same rule as the blastp ortholog assignment).
    hh2 <- hh[hh$length >= 60, ]; if (!nrow(hh2)) hh2 <- hh
    anchor <- hh2[order(-hh2$pident, -hh2$bitscore), ][1, ]
    cl <- hh[hh$chrn==anchor$chrn & hh$lo <= anchor$hi+50000 & hh$hi >= anchor$lo-50000, ]
    asn$tblastn_chr[i] <- anchor$chrn
    asn$tblastn_start[i] <- min(cl$lo); asn$tblastn_end[i] <- max(cl$hi)
    # overlap with GFF3 genes on that chr -> tblastn_horvu (max overlap with anchor midpoint)
    mid <- (anchor$lo+anchor$hi)/2
    gc <- gg[gg$chrn==anchor$chrn & gg$start<=mid & gg$end>=mid, ]
    if (nrow(gc)) asn$tblastn_horvu[i] <- gc$id[1]
  }
}
asn$tblastn_vs_blastp_agree <- ifelse(is.na(asn$tblastn_horvu), "NA",
                                ifelse(asn$tblastn_horvu==asn$horvu_id, "Y", "N"))

# ---- lead SNPs, tested list, haplotype FDR ----
lead <- read.delim(LEAD)
win  <- read.delim(WIN)
hap  <- read.csv(HAP)
hap_key <- paste(hap$trait, hap$gene_id)

dist_to_leads <- function(chr, gstart, gend, trait) {
  ld <- lead[lead$trait==trait & lead$chr==chr, ]
  if (!nrow(ld)) return(list(snp=NA, d=NA_real_))
  d <- ifelse(ld$lead_pos>=gstart & ld$lead_pos<=gend, 0,
              pmin(abs(ld$lead_pos-gstart), abs(ld$lead_pos-gend)))
  j <- which.min(d); list(snp=ld$lead_SNP[j], d=d[j])
}
win_flag <- function(d) ifelse(is.na(d), "FAR",
                        ifelse(d<=200000, "IN_200kb_WINDOW",
                        ifelse(d<=1000000, "NEAR_200kb_to_1Mb", "FAR")))

rows <- list()
for (i in seq_len(nrow(asn))) {
  a <- asn[i, ]
  for (tr in TRAITS) {
    nd <- dist_to_leads(a$chr, a$gene_start, a$gene_end, tr)
    tested <- a$horvu_id %in% win$gene_id[win$trait==tr]
    k <- paste(tr, a$horvu_id)
    sig <- NA; pfdr <- NA
    if (k %in% hap_key) {
      hr <- hap[hap_key==k, ][1, ]
      sig <- as.character(hr$significant_fdr); pfdr <- hr$trait_fdr_p
    }
    rows[[length(rows)+1]] <- data.frame(
      canonical_gene=a$canonical_gene, source_accession=a$source_accession,
      horvu_id=a$horvu_id, chr=a$chr, gene_start=a$gene_start, gene_end=a$gene_end,
      strand=a$strand, pct_identity=a$pct_identity, pct_coverage=a$pct_coverage,
      evalue=a$evalue, reciprocal_best_hit=a$reciprocal_best_hit,
      tblastn_horvu=a$tblastn_horvu, tblastn_vs_blastp_agree=a$tblastn_vs_blastp_agree,
      trait=tr, nearest_lead_SNP=nd$snp, distance_bp=nd$d,
      window_flag=win_flag(nd$d),
      tested_in_candidate_list=ifelse(tested,"Y","N"),
      haplotype_significant_fdr=sig, haplotype_trait_fdr_p=pfdr)
  }
}
res <- do.call(rbind, rows)
res <- res[order(res$canonical_gene, res$trait), ]
write.table(res, OUT, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

# ---- report ----
cat("=== Step 06: canonical gene distances ===\n")
cat("rows:", nrow(res), " (", length(unique(res$canonical_gene)), "genes x", length(TRAITS),"traits )\n\n")
u <- res[!duplicated(res$canonical_gene), c("canonical_gene","horvu_id","chr","gene_start","pct_identity","reciprocal_best_hit","tblastn_vs_blastp_agree")]
cat("assignment summary:\n"); print(u, row.names=FALSE)
cat("\nwindow_flag counts:\n"); print(table(res$window_flag))
cat("\nany IN_200kb_WINDOW:\n")
iw <- res[res$window_flag=="IN_200kb_WINDOW", c("canonical_gene","trait","nearest_lead_SNP","distance_bp","tested_in_candidate_list","haplotype_significant_fdr")]
if (nrow(iw)) print(iw, row.names=FALSE) else cat("  none\n")
cat("\nnearest approach per gene (min distance across the 3 traits):\n")
mn <- aggregate(distance_bp~canonical_gene, res, function(x) min(x,na.rm=TRUE))
print(mn[order(mn$distance_bp),], row.names=FALSE)
cat("\nWrote:", OUT, "\n")
