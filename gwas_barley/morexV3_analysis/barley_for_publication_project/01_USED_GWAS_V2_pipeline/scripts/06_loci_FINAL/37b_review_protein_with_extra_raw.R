#!/usr/bin/env Rscript
# 37b_review_protein_with_extra_raw.R
# REVIEW ONLY -- companion to 37_review_raw_no_gap.R.
#
# Protein Manhattan with NO gap rule, including the two sub-threshold 3H peaks as
# their own index SNPs. Each entry spans the full min->max of every SNP at
# r^2 >= 0.5 within +/-2 Mb of its anchor, however large the internal gaps.
#
# The two protein clumps are read from 02_loci_FINAL/plink/pass1/protein.clumped and
# the two extra peaks from 02_loci_FINAL/extra/*.ld, so this view uses exactly the
# same underlying LD as the analysis set -- only the severing is omitted.
#
# The extras keep their established colours (3H:173351655 green #00A878,
# 3H:198076308 amber #E8A33D) and are drawn as diamonds.
#
# Output -> 02b_REVIEW_raw_no_gap/figures/manhattan_protein__raw_nogap_with_extra.{pdf,png}
# Created 2026-09-09.
Sys.setenv(TMPDIR="/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(png) })
PIPE  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
INTER <- file.path(PIPE,"intermediates"); PS <- file.path(PIPE,"results","emmax_ps")
SET   <- file.path(PIPE,"results","00_FINAL_BLUP_3PC")
LOC   <- file.path(SET,"02_loci_FINAL")
OUT   <- file.path(SET,"02b_REVIEW_raw_no_gap"); FIG <- file.path(OUT,"figures"); TAB <- file.path(OUT,"tables")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

BONF <- -log10(0.10/length(readLines(file.path(INTER,"morexV3_pruned_for_covs.prune.in"))))
BG <- c("grey72","grey55"); H_MAN <- 7; RDPI <- 900; PDPI <- 400
W_IN <- 18; CEX_LEG <- 1.05; CEX_FOOT <- 0.85; ROW_IN <- 0.42
EXTRA_COL <- c("3H:173351655"="#00A878","3H:198076308"="#E8A33D")
LOC_COL   <- c("#E41A1C","#377EB8")

sm <- fread(file.path(INTER,"morexV3_290.bim"), header=FALSE,
            col.names=c("CHR_raw","SNP","cm","BP","A1","A2"))
sm[, CHR := as.integer(sub("H$","",CHR_raw))]; setorder(sm, CHR, BP)
cmx <- sm[, .(mx=max(BP)), by=CHR][order(CHR)]; cmx[, off := cumsum(as.numeric(shift(mx, fill=0)))]
sm <- merge(sm, cmx[, .(CHR,off)], by="CHR"); setorder(sm, CHR, BP)
sm[, x := as.numeric(BP)+off][, col := ifelse(CHR%%2L==0L, BG[2], BG[1])]
cc <- sm[, .(c=(min(x)+max(x))/2), by=CHR][order(CHR)]
XL <- c(min(sm$x), max(sm$x)); setkey(sm, SNP)

gw <- fread(file.path(PS,"morexV3__protein__BLUP__pc3.ps"), header=FALSE,
            col.names=c("SNP","beta","SE","P"), showProgress=FALSE)
d <- data.table(SNP=sm$SNP, x=sm$x, col=sm$col,
                nlp=-log10(gw$P[match(sm$SNP, gw$SNP)]))[is.finite(nlp)]
ymax <- ceiling(max(max(d$nlp), BONF, 6))+1; setkey(d, SNP)

E <- list(); MM <- list()
# the 2 significant protein clumps, raw
cl <- fread(file.path(LOC,"plink","pass1","protein.clumped"), fill=TRUE); cl <- cl[!is.na(BP) & SNP!=""]
setorder(cl, CHR, BP)
for (i in seq_len(nrow(cl))) {
  sp <- cl$SP2[i]
  ms <- if (is.na(sp)||sp %chin% c("NONE","")) character(0) else trimws(gsub("\\(1\\)","",strsplit(sp,",")[[1]]))
  ms <- ms[grepl(":", ms, fixed=TRUE)]
  bp <- suppressWarnings(as.numeric(sub(".*:","", c(cl$SNP[i], ms)))); bp <- sort(unique(bp[is.finite(bp)]))
  id <- sprintf("protein_R%02d", i)
  E[[id]] <- data.table(locus_id=id, kind="locus", lead_SNP=cl$SNP[i], chr=as.character(cl$CHR[i]),
    lead_bp=as.numeric(cl$BP[i]), start=min(bp), end=max(bp), span_kb=(max(bp)-min(bp))/1000,
    n_members=length(bp), max_internal_gap_kb=if(length(bp)>1) max(diff(bp))/1000 else 0,
    colour=LOC_COL[i])
  MM[[id]] <- data.table(locus_id=id, SNP=sprintf("%s:%d", as.character(cl$CHR[i]), bp), colour=LOC_COL[i])
}
# the 2 sub-threshold peaks, raw
for (s in names(EXTRA_COL)) {
  f <- file.path(LOC,"extra",paste0(sub(":","_",s),".ld")); stopifnot(file.exists(f))
  bp <- sort(unique(fread(f)$BP_B)); a <- as.numeric(sub(".*:","",s))
  E[[s]] <- data.table(locus_id=s, kind="extra", lead_SNP=s, chr="3H", lead_bp=a,
    start=min(bp), end=max(bp), span_kb=(max(bp)-min(bp))/1000, n_members=length(bp),
    max_internal_gap_kb=if(length(bp)>1) max(diff(bp))/1000 else 0, colour=EXTRA_COL[[s]])
  MM[[s]] <- data.table(locus_id=s, SNP=sprintf("3H:%d", bp), colour=EXTRA_COL[[s]])
}
S <- rbindlist(E); setorder(S, chr, lead_bp)
S[, lead_x := sm[lead_SNP, x]]
S[, leg := sprintf("%-13s %s:%s-%s  |  %9s kb, %5d SNPs  |  lead %s  %s",
     locus_id, chr, format(start, big.mark=","), format(end, big.mark=","),
     formatC(span_kb, format="f", digits=1, big.mark=","), n_members,
     format(lead_bp, big.mark=","), fifelse(kind=="locus","[significant]","[SUB-THRESHOLD]"))]
M <- rbindlist(MM); M[, `:=`(x=sm[SNP,x], nlp=d[SNP,nlp])]; M <- M[is.finite(nlp) & is.finite(x)]
fwrite(S[, .(locus_id, kind, lead_SNP, chr, lead_bp, start, end, span_kb, n_members, max_internal_gap_kb)],
       file.path(TAB,"protein_with_extra_raw_summary.tsv"), sep="\t")

n <- nrow(S); H_LEG <- ROW_IN*(n+2L)+0.45; H_IN <- H_MAN+H_LEG
draw <- function() {
  layout(matrix(1:2, nrow=2), heights=c(H_MAN, H_LEG))
  par(mar=c(4.5,5.2,3,1.2), las=1, cex.axis=1.15, cex.lab=1.2, cex.main=1.0)
  plot(NA, xlim=XL, ylim=c(0,ymax), xaxs="i", yaxs="i", axes=FALSE,
       xlab="Chromosome", ylab=expression(-log[10](italic(p))),
       main="Protein | REVIEW ONLY -- NO GAP RULE | r2>=0.5, +/-2 Mb, span = full min-max of every member  +  the 2 sub-threshold 3H peaks")
  rasterImage(img, XL[1], 0, XL[2], ymax, interpolate=FALSE)
  axis(2); axis(1, at=cc$c, labels=paste0(cc$CHR,"H"), tick=FALSE); box()
  abline(h=BONF, lwd=1.4)
  points(M$x, M$nlp, pch=16, cex=1.05, col=M$colour)
  L <- S[kind=="locus"]; X <- S[kind=="extra"]
  points(L$lead_x, d[L$lead_SNP,nlp], pch=16, cex=1.6, col=L$colour)
  points(X$lead_x, d[X$lead_SNP,nlp], pch=18, cex=2.2, col=X$colour)
  legend("topleft", bty="n", cex=1.0, pch=c(16,16,18), col="grey20", pt.cex=c(1.2,1.8,2.5),
         legend=c("member (r2 >= 0.5, NOT severed)","lead SNP (significant)",
                  "sub-threshold peak used as its own index SNP"))
  par(mar=c(0.2,0.6,0.2,0.6)); plot.new(); plot.window(c(0,1),c(0,1))
  for (i in seq_len(n)) {
    yy <- 1-(i-0.7)/(n+2L)
    points(0.004, yy, pch=if (S$kind[i]=="locus") 16 else 18, col=S$colour[i],
           cex=if (S$kind[i]=="locus") 1.7 else 2.3, xpd=NA)
    text(0.020, yy, S$leg[i], adj=c(0,0.5), cex=CEX_LEG, xpd=NA, family="mono")
  }
  text(0, 1-(n+1.4)/(n+2L),
    "REVIEW ONLY. No contiguity rule: each entry spans from its first to its last absorbed member, however large the internal gaps. The analysis set is 02_loci_FINAL/.",
    adj=c(0,0), cex=CEX_FOOT, col="grey35", xpd=NA)
}
tmp <- tempfile(fileext=".png", tmpdir=Sys.getenv("TMPDIR"))
png(tmp, width=W_IN, height=H_MAN, units="in", res=RDPI, bg="transparent", type="cairo-png")
par(mar=c(0,0,0,0)); plot(d$x, d$nlp, xlim=XL, ylim=c(0,ymax), xaxs="i", yaxs="i",
     axes=FALSE, xlab="", ylab="", pch=16, cex=0.85, col=d$col); dev.off()
img <- readPNG(tmp); unlink(tmp)
b <- "manhattan_protein__raw_nogap_with_extra"
cairo_pdf(file.path(FIG, paste0(b,".pdf")), width=W_IN, height=H_IN, pointsize=10); draw(); dev.off()
png(file.path(FIG, paste0(b,".png")), width=W_IN, height=H_IN, units="in", res=PDPI,
    pointsize=10, type="cairo-png"); draw(); dev.off()
cat(sprintf("[37b] %d entries (%d loci + %d extras) | %d SNPs painted\n",
            n, sum(S$kind=="locus"), sum(S$kind=="extra"), nrow(M)))
print(S[, .(locus_id, span_kb=round(span_kb,1), n_members, maxgap_kb=round(max_internal_gap_kb,1))])
