#!/usr/bin/env Rscript
# 02_local_association_profile.R
# Per-SNP association in and around the gene, for BOTH traits, plus MAF.
# KEY RESULT: the gene BODY is null in both traits; all signal sits 554-763 bp
# upstream of the TSS (promoter side). That is a cis-regulatory signature, and it
# is why a coding-variant interpretation is not supported.
# Outputs -> results/tables/{local_association_profile.tsv, signal_snps.tsv}
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("BRANCH"); OUT <- file.path(B,"results","tables")
GS <- as.numeric(Sys.getenv("GENE_START")); GE <- as.numeric(Sys.getenv("GENE_END"))
TSS <- as.numeric(Sys.getenv("GENE_TSS")); THR <- as.numeric(Sys.getenv("BONF_THRESHOLD"))
A <- Sys.getenv("ASSOC_DIR"); lo <- GS-1000; hi <- GE+3000

grab <- function(tr) {
  d <- fread(cmd=sprintf("awk -F'\\t' 'NR>1{split($1,x,\":\"); if(x[1]==\"7H\" && x[2]+0>=%d && x[2]+0<=%d) print x[2]\"\\t\"$2}' %s/%s.assoc", lo, hi, A, tr),
             header=FALSE, col.names=c("pos","p"))
  d[, neglog10p := -log10(p)][, trait := tr][]
}
prof <- dcast(rbindlist(lapply(c("fiber","starch"), grab)), pos ~ trait, value.var="neglog10p")
setnames(prof, c("fiber","starch"), c("fiber_neglog10p","starch_neglog10p"))
prof[, region := fifelse(pos>=GS & pos<=GE, "GENE_BODY", fifelse(pos>TSS, "PROMOTER_upstream", "downstream_of_3prime"))]
prof[, offset_from_TSS := pos - TSS]
frq <- fread(Sys.getenv("FREQ_FRQ"))[, .(SNP, MAF, NCHROBS)]
frq[, pos := as.numeric(sub(".*:","",SNP))]
prof <- merge(prof, frq[,.(pos,MAF,NCHROBS)], by="pos", all.x=TRUE)
setorder(prof, pos)
prof[, sig_fiber := fiber_neglog10p > THR][, sig_starch := starch_neglog10p > THR]
fwrite(prof, file.path(OUT,"local_association_profile.tsv"), sep="\t")

sig <- prof[fiber_neglog10p > 2 | starch_neglog10p > 2]
fwrite(sig, file.path(OUT,"signal_snps.tsv"), sep="\t")

cat(sprintf("02 OK: %d SNPs profiled | GENE BODY max -log10p: fiber %.2f, starch %.2f\n",
   nrow(prof), max(prof[region=="GENE_BODY"]$fiber_neglog10p), max(prof[region=="GENE_BODY"]$starch_neglog10p)))
cat(sprintf("        promoter-side max: fiber %.2f, starch %.2f\n",
   max(prof[region=="PROMOTER_upstream"]$fiber_neglog10p),
   max(prof[region=="PROMOTER_upstream"]$starch_neglog10p)))
