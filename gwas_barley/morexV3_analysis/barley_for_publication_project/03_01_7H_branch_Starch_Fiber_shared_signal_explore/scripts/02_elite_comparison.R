#!/usr/bin/env Rscript
# 02_elite_comparison.R -- EXPLORATORY: which haplotype of 7HG0729030 do modern elite
# malting cultivars carry?  Method follows step 07 (03_build_matrices.R / 04_figures.R):
#   * wild haplotype groups are CONSUMED from the cached HapObject, never recomputed,
#     and the elite lines are never shown to crosshap;
#   * each wild group is reduced to one row by per-SNP MAJORITY CONSENSUS; an exact tie
#     resolves to REF when REF is one of the tied states, otherwise NA -- step 07's rule
#     since 2026-09-24 (aligned here 2026-09-27, user approval; was: all ties -> NA. This
#     gene has no tie, so no output changed);
#   * shared sites are matched on CHROM:POS:REF:ALT and REF/ALT identity is asserted.
# Difference from step 07: haplotypes come from MGmin = 3 (not the pipeline's MGmin = 2,
# under which this gene has no usable grouping at all).
suppressPackageStartupMessages({library(vcfR); library(data.table); library(ggplot2)})
T <- Sys.getenv("TEMP_ROOT"); P <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
E <- file.path(P,"07_USED_elite_lines_compariosn_to_wild_lines")
GENE <- "HORVU.MOREX.r3.7HG0729030"; SIG <- c(573606282,573606306,573606460,573606491)
COL <- c(REF="#FFFACD", ALT="#2F4F4F", MISS="grey70", HET="#C46210")   # step 07's palette
elite_cfg <- fread(cmd=paste("grep -v '^#'", shQuote(file.path(E,"config/elite_lines.tsv"))))
gt_code <- function(gt){o <- rep(NA_character_,length(gt))
  o[gt %in% c("0/0","0|0")] <- "REF"; o[gt %in% c("1/1","1|1")] <- "ALT"
  o[gt %in% c("0/1","1/0","0|1","1|0")] <- "HET"; o}
read_gt <- function(p){ v <- suppressMessages(read.vcfR(p, verbose=FALSE))
  fx <- as.data.frame(getFIX(v), stringsAsFactors=FALSE)
  g  <- extract.gt(v, element="GT", as.numeric=FALSE)
  keep <- nchar(fx$REF)==1 & nchar(fx$ALT)==1     # biallelic SNVs only
  list(key=paste(fx$CHROM[keep],fx$POS[keep],fx$REF[keep],fx$ALT[keep],sep=":"),
       pos=as.integer(fx$POS[keep]), gt=g[keep,,drop=FALSE]) }
wild  <- read_gt(file.path(T,"work/raw/fiber",paste0(GENE,".vcf.gz")))
elite <- read_gt(file.path(T,"work/elite/7HG0729030.elite.vcf.gz"))
elite$key <- sub("^chr","",elite$key)
shared <- intersect(wild$key, elite$key)
cat(sprintf("wild SNVs %d | elite SNVs %d | SHARED %d | of the 4 signal SNPs shared: %d\n",
    length(wild$key), length(elite$key), length(shared),
    sum(as.integer(sub("^[^:]+:([0-9]+):.*","\\1",shared)) %in% SIG)))
# allele concordance: same key => same REF/ALT by construction; report positions differing
wpos <- wild$pos[match(shared, wild$key)]
out <- list()
for (tr in c("fiber","starch")) for (EPS in c(0.9)) {   # epsilon 0.6 dropped 2026-09-24 (user)
  ho <- readRDS(file.path(T,"Cache",tr,GENE,"MGmin_3",paste0("eps_",EPS),"HapObject.rds"))
  ind <- as.data.frame(ho$HapObject[[paste0("Haplotypes_MGmin3_E",EPS)]]$Indfile)
  ind$Ind <- as.character(ind$Ind); ind$hap <- as.character(ind$hap)
  ind <- ind[ind$hap!="0",]
  wm <- wild$gt[match(shared, wild$key), , drop=FALSE]
  cons <- sapply(sort(unique(ind$hap)), function(h){
    sub <- wm[, colnames(wm) %in% ind$Ind[ind$hap==h], drop=FALSE]
    apply(sub, 1, function(r){ tb <- table(gt_code(r)); if(!length(tb)) return(NA_character_)
      top <- names(tb)[tb==max(tb)]
      if (length(top)>1) return(if ("REF" %in% top) "REF" else NA_character_); top }) })
  cons <- as.data.frame(cons, stringsAsFactors=FALSE)
  em <- elite$gt[match(shared, elite$key), , drop=FALSE]
  em <- em[, colnames(em) %in% elite_cfg$sample_id, drop=FALSE]
  ecode <- apply(em, 2, gt_code)
  colnames(ecode) <- elite_cfg$line_name[match(colnames(em), elite_cfg$sample_id)]
  n_h <- ncol(cons)
  D <- rbind(
    data.table(pos=wpos, row=rep(paste0("wild hap ", colnames(cons),
                 " (n=", sapply(colnames(cons), function(h) sum(ind$hap==h)), ")"), each=length(shared)),
               call=unlist(cons), block="wild"),
    data.table(pos=wpos, row=rep(colnames(ecode), each=length(shared)),
               call=as.vector(ecode), block="elite"))
  D[is.na(call), call := "MISS"]
  D[, signal := pos %in% SIG]
  lev <- c(rev(sort(unique(D[block=="wild"]$row))), rev(colnames(ecode)))
  D[, row := factor(row, levels=lev)]
  D[, x := as.integer(factor(pos, levels=sort(unique(pos))))]
  D[, trait := tr][, eps := EPS]
  out[[length(out)+1]] <- D
  g <- ggplot(D, aes(x, row, fill=call)) +
    geom_tile(colour="white", linewidth=0.25, height=0.72) +
    geom_point(data=unique(D[signal==TRUE, .(x, y=0.3)]), aes(x=x, y=y), inherit.aes=FALSE,
               shape=17, size=2.2, colour="red") +
    scale_fill_manual(values=COL, name="genotype") +
    annotate("segment", x=0.5, xend=length(shared)+0.5, y=n_h+0.5, yend=n_h+0.5, linewidth=0.7) +
    labs(title=sprintf("7HG0729030 (GDSL esterase) - wild haplotypes vs elite cultivars | %s", tr),
         subtitle=sprintf("MGmin=3, eps=%s [EXPLORATORY]  |  %d shared SNVs  |  red triangles = the 4 GWAS signal SNPs",
                          EPS, length(shared)),
         x=sprintf("shared SNV (7H:%s-%s)", min(wpos), max(wpos)), y=NULL) +
    theme_minimal(base_size=9) +
    theme(axis.text.x=element_blank(), panel.grid=element_blank(),
          axis.text.y=element_text(family="mono"), plot.title=element_text(face="bold"))
  ggsave(file.path(T,"results/figures/elite",sprintf("elite_vs_wild__7HG0729030__%s__MGmin3_eps%s.pdf",tr,EPS)),
         g, width=11, height=2.2+0.32*(n_h+ncol(ecode)), device=grDevices::cairo_pdf)
  ggsave(file.path(T,"results/figures/elite",sprintf("elite_vs_wild__7HG0729030__%s__MGmin3_eps%s.png",tr,EPS)),
         g, width=11, height=2.2+0.32*(n_h+ncol(ecode)), dpi=300)
}
A <- rbindlist(out); fwrite(A, file.path(T,"results/tables/elite_vs_wild_barcodes.tsv"), sep="\t")
cat("\n=== elite genotypes AT THE 4 GWAS SIGNAL SNPs (eps=0.9) ===\n")
z <- A[eps==0.9 & trait=="fiber" & signal==TRUE, .(row, pos, call)]
print(dcast(z, row ~ pos, value.var="call"))
