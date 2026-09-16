#!/usr/bin/env Rscript
# 37_review_raw_no_gap.R
# REVIEW ONLY -- not part of the analysis set.
#
# Shows what the loci look like with the SAME clumping but NO gap rule: each locus is
# simply the full min->max span of every SNP --clump absorbed, however scattered.
# Provided so the effect of the gap rule can be judged by eye against
# 02_loci_FINAL/figures/.
#
# Clumping is identical to the live definition -- p1=9.008e-07, p2=1 (membership on
# LD alone), r2=0.5, kb=2000 (+/-2 Mb) -- and is READ FROM plink/pass1/*.clumped of
# that run, so no PLINK re-run is needed and the two views are guaranteed to come
# from the same clumping.
#
# No iteration is required here: without severing, every SNP --clump touched stays in
# its clump, so no significant SNP can be orphaned.
#
# Legend matches the live figures: locus id, span coordinates (chr:start-end), span
# in kb, member count, lead SNP position.
#
# Outputs -> 02b_REVIEW_raw_no_gap/{tables/loci_raw_summary.tsv,
#                                   figures/manhattan_<trait>__raw_nogap.{pdf,png}}
# Created 2026-09-09.
Sys.setenv(TMPDIR="/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(png) })
PIPE  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
INTER <- file.path(PIPE,"intermediates"); PS <- file.path(PIPE,"results","emmax_ps")
SET   <- file.path(PIPE,"results","00_FINAL_BLUP_3PC")
CLP   <- file.path(SET,"02_loci_FINAL","plink","pass1")
OUT   <- file.path(SET,"02b_REVIEW_raw_no_gap")
TAB   <- file.path(OUT,"tables"); FIG <- file.path(OUT,"figures")
dir.create(TAB, showWarnings=FALSE, recursive=TRUE); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

BONF <- -log10(0.10/length(readLines(file.path(INTER,"morexV3_pruned_for_covs.prune.in"))))
TRAITS <- c("betaglucan","fiber","protein","starch")
TRAIT_LAB <- c(betaglucan="beta-glucan", fiber="Fiber", protein="Protein", starch="Starch")
BG <- c("grey72","grey55"); H_MAN <- 7; RDPI <- 900; PDPI <- 400
W_IN <- 18; CEX_LEG <- 1.05; CEX_FOOT <- 0.85; ROW_IN <- 0.42
PALETTE <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF",
             "#1B9E77","#666666","#66A61E","#E6AB02","#7570B3","#00CED1","#B15928","#8DD3C7")

sm <- fread(file.path(INTER,"morexV3_290.bim"), header=FALSE,
            col.names=c("CHR_raw","SNP","cm","BP","A1","A2"))
sm[, CHR := as.integer(sub("H$","",CHR_raw))]; setorder(sm, CHR, BP)
cmx <- sm[, .(mx=max(BP)), by=CHR][order(CHR)]; cmx[, off := cumsum(as.numeric(shift(mx, fill=0)))]
sm <- merge(sm, cmx[, .(CHR,off)], by="CHR"); setorder(sm, CHR, BP)
sm[, x := as.numeric(BP)+off][, col := ifelse(CHR%%2L==0L, BG[2], BG[1])]
cc <- sm[, .(c=(min(x)+max(x))/2), by=CHR][order(CHR)]
XL <- c(min(sm$x), max(sm$x)); setkey(sm, SNP)

allS <- list()
for (tr in TRAITS) {
  gw <- fread(file.path(PS, sprintf("morexV3__%s__BLUP__pc3.ps", tr)), header=FALSE,
              col.names=c("SNP","beta","SE","P"), showProgress=FALSE)
  d <- data.table(SNP=sm$SNP, x=sm$x, col=sm$col,
                  nlp=-log10(gw$P[match(sm$SNP, gw$SNP)]))[is.finite(nlp)]
  ymax <- ceiling(max(max(d$nlp), BONF, 6))+1; setkey(d, SNP)

  cl <- fread(file.path(CLP, paste0(tr,".clumped")), fill=TRUE); cl <- cl[!is.na(BP) & SNP!=""]
  setorder(cl, CHR, BP)
  K <- list(); MM <- list()
  for (i in seq_len(nrow(cl))) {
    sp <- cl$SP2[i]
    ms <- if (is.na(sp) || sp %chin% c("NONE","")) character(0) else trimws(gsub("\\(1\\)","",strsplit(sp,",")[[1]]))
    ms <- ms[grepl(":", ms, fixed=TRUE)]
    bp <- suppressWarnings(as.numeric(sub(".*:","", c(cl$SNP[i], ms))))
    bp <- sort(unique(bp[is.finite(bp)]))                 # NO gap rule: keep everything
    id <- sprintf("%s_R%02d", tr, i)
    K[[i]] <- data.table(trait=tr, locus_id=id, lead_SNP=cl$SNP[i], chr=as.character(cl$CHR[i]),
      lead_bp=as.numeric(cl$BP[i]), start=min(bp), end=max(bp),
      span_kb=(max(bp)-min(bp))/1000, n_members=length(bp),
      max_internal_gap_kb=if (length(bp)>1) max(diff(bp))/1000 else 0)
    MM[[i]] <- data.table(locus_id=id, SNP=sprintf("%s:%d", as.character(cl$CHR[i]), bp))
  }
  K <- rbindlist(K); K[, locus_id := sprintf("%s_R%02d", trait, seq_len(.N))]
  K[, colour := PALETTE[((seq_len(.N)-1L) %% length(PALETTE))+1L]]
  K[, lead_x := sm[lead_SNP, x]]
  K[, leg := sprintf("%-14s %s:%s-%s  |  %9s kb, %5d SNPs  |  lead %s",
       locus_id, chr, format(start, big.mark=","), format(end, big.mark=","),
       formatC(span_kb, format="f", digits=1, big.mark=","), n_members,
       format(lead_bp, big.mark=","))]
  M <- rbindlist(MM); M[, locus_id := rep(K$locus_id, K$n_members)]
  M <- merge(M, K[, .(locus_id, colour)], by="locus_id")
  M[, `:=`(x=sm[SNP,x], nlp=d[SNP,nlp])]; M <- M[is.finite(nlp) & is.finite(x)]
  allS[[tr]] <- K[, .(trait, locus_id, lead_SNP, chr, lead_bp, start, end, span_kb,
                      n_members, max_internal_gap_kb)]

  n <- nrow(K); NC <- if (n > 8) 2L else 1L; nr <- ceiling(n/NC)
  H_LEG <- ROW_IN*(nr+1L)+0.45; H_IN <- H_MAN+H_LEG
  draw <- function() {
    layout(matrix(1:2, nrow=2), heights=c(H_MAN, H_LEG))
    par(mar=c(4.5,5.2,3,1.2), las=1, cex.axis=1.15, cex.lab=1.2, cex.main=1.05)
    plot(NA, xlim=XL, ylim=c(0,ymax), xaxs="i", yaxs="i", axes=FALSE,
         xlab="Chromosome", ylab=expression(-log[10](italic(p))),
         main=sprintf("%s | REVIEW ONLY -- NO GAP RULE | %d clumps: --clump r2>=0.5, +/-2 Mb, span = full min-max of every member absorbed",
                      TRAIT_LAB[[tr]], n))
    rasterImage(img, XL[1], 0, XL[2], ymax, interpolate=FALSE)
    axis(2); axis(1, at=cc$c, labels=paste0(cc$CHR,"H"), tick=FALSE); box()
    abline(h=BONF, lwd=1.4)
    points(M$x, M$nlp, pch=16, cex=1.05, col=M$colour)
    points(K$lead_x, d[K$lead_SNP, nlp], pch=16, cex=1.6, col=K$colour)
    legend("topleft", bty="n", cex=1.0, pch=c(16,16), col="grey20", pt.cex=c(1.2,1.8),
           legend=c("clump member (r2 >= 0.5, NOT severed)","lead SNP"))
    par(mar=c(0.2,0.6,0.2,0.6)); plot.new(); plot.window(c(0,1),c(0,1))
    for (i in seq_len(n)) {
      ci <- (i-1L)%/%nr; ri <- (i-1L)%%nr
      xx <- ci/NC+0.004; yy <- 1-(ri+0.70)/(nr+1L)
      points(xx, yy, pch=16, col=K$colour[i], cex=1.7, xpd=NA)
      text(xx+0.016, yy, K$leg[i], adj=c(0,0.5), cex=CEX_LEG, xpd=NA, family="mono")
    }
    text(0, 1-(nr+0.55)/(nr+1L),
      "REVIEW ONLY. No contiguity rule is applied: a locus spans from its first to its last absorbed member, however large the internal gaps. The analysis set is 02_loci_FINAL/.",
      adj=c(0,0), cex=CEX_FOOT, col="grey35", xpd=NA)
  }
  tmp <- tempfile(fileext=".png", tmpdir=Sys.getenv("TMPDIR"))
  png(tmp, width=W_IN, height=H_MAN, units="in", res=RDPI, bg="transparent", type="cairo-png")
  par(mar=c(0,0,0,0)); plot(d$x, d$nlp, xlim=XL, ylim=c(0,ymax), xaxs="i", yaxs="i",
       axes=FALSE, xlab="", ylab="", pch=16, cex=0.85, col=d$col); dev.off()
  img <- readPNG(tmp); unlink(tmp)
  b <- sprintf("manhattan_%s__raw_nogap", tr)
  cairo_pdf(file.path(FIG, paste0(b,".pdf")), width=W_IN, height=H_IN, pointsize=10); draw(); dev.off()
  png(file.path(FIG, paste0(b,".png")), width=W_IN, height=H_IN, units="in", res=PDPI,
      pointsize=10, type="cairo-png"); draw(); dev.off()
  cat(sprintf("[37] %-11s %2d clumps | %5d SNPs painted | median span %8.1f kb | max %8.1f kb\n",
              tr, n, nrow(M), median(K$span_kb), max(K$span_kb)))
  rm(gw, d, img); gc(verbose=FALSE)
}
S <- rbindlist(allS); setorder(S, trait, chr, lead_bp)
fwrite(S, file.path(TAB,"loci_raw_summary.tsv"), sep="\t")
cat(sprintf("\n[37] %d clumps | median span %.1f kb | mean %.1f kb | total %.1f Mb\n",
            nrow(S), median(S$span_kb), mean(S$span_kb), sum(S$span_kb)/1000))
cat(sprintf("[37] OK -> %s\n", OUT))
