#!/usr/bin/env Rscript
# 33_manhattan_loci_painted.R
# Per-trait Manhattan plots with every locus member painted in that locus's colour,
# so each locus reads as one uniform colour and its extent is visible directly. The
# lead SNP is the same colour, slightly larger; a locus whose lead is its only
# surviving member is drawn as an open circle.
#
# Loci come from 31_loci_clump_iterative.R. All parameter strings shown in the title
# and footnote are read from tables/locus_definition_params.tsv, so a caption can
# never state a value the run did not use.
#
# The legend lists, per locus: id, span coordinates (chr:start-end), span in kb,
# member count, and the lead SNP position.
#
# Output -> 02_loci_FINAL/figures/manhattan_<trait>__loci_r05.{pdf,png}
#           02_loci_FINAL/figures/locus_colour_key__<trait>.tsv
Sys.setenv(TMPDIR="/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(png) })
PIPE  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
INTER <- file.path(PIPE,"intermediates"); PS <- file.path(PIPE,"results","emmax_ps")
SET   <- file.path(PIPE,"results","00_FINAL_BLUP_3PC")
LOC   <- file.path(SET,"02_loci_FINAL")
OUT   <- file.path(SET,"02_loci_FINAL","figures")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

BONF <- -log10(0.10 / length(readLines(file.path(INTER,"morexV3_pruned_for_covs.prune.in"))))
TRAITS <- c("betaglucan","fiber","protein","starch")
TRAIT_LAB <- c(betaglucan="beta-glucan", fiber="Fiber", protein="Protein", starch="Starch")
BG_COLS <- c("grey72","grey55"); H_MAN <- 7; RASTER_DPI <- 900; PNG_DPI <- 400
# Legend text sizing (raised 2026-09-09: the previous cex 0.62 / 0.55 was unreadable
# at normal viewing size). Row height and canvas width scale with it so the lines
# still fit without wrapping.
W_IN <- 18; CEX_LEG <- 1.05; CEX_FOOT <- 0.85; ROW_IN <- 0.42
PALETTE <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF",
             "#1B9E77","#666666","#66A61E","#E6AB02","#7570B3","#00CED1","#B15928","#8DD3C7")

snp_map <- fread(file.path(INTER,"morexV3_290.bim"), header=FALSE,
                 col.names=c("CHR_raw","SNP","cm","BP","A1","A2"))
snp_map[, CHR := as.integer(sub("H$","",CHR_raw))]; setorder(snp_map, CHR, BP)
cm <- snp_map[, .(mx=max(BP)), by=CHR][order(CHR)]
cm[, off := cumsum(as.numeric(shift(mx, fill=0)))]
snp_map <- merge(snp_map, cm[, .(CHR,off)], by="CHR"); setorder(snp_map, CHR, BP)
snp_map[, x := as.numeric(BP)+off]
snp_map[, col := ifelse(CHR %% 2L == 0L, BG_COLS[2], BG_COLS[1])]
cc <- snp_map[, .(c=(min(x)+max(x))/2), by=CHR][order(CHR)]
XLIM <- c(min(snp_map$x), max(snp_map$x)); setkey(snp_map, SNP)

TAB_PARAM <- file.path(LOC,"tables")
PAR <- fread(file.path(TAB_PARAM,"locus_definition_params.tsv"))
gp  <- function(k) PAR[parameter==k, value][1]
P_KB <- as.numeric(gp("clump_kb")); P_GAP <- as.numeric(gp("gap_rule_kb"))
P_R2 <- as.numeric(gp("clump_r2")); P_P1 <- gp("clump_p1")
S <- fread(file.path(LOC,"tables","loci_summary.tsv")); M <- fread(file.path(LOC,"tables","loci_members.tsv"))

for (tr in TRAITS) {
  gw <- fread(file.path(PS, sprintf("morexV3__%s__BLUP__pc3.ps", tr)), header=FALSE,
              col.names=c("SNP","beta","SE","P"), showProgress=FALSE)
  d <- data.table(SNP=snp_map$SNP, x=snp_map$x, col=snp_map$col,
                  nlp=-log10(gw$P[match(snp_map$SNP, gw$SNP)]))[is.finite(nlp)]
  ymax <- ceiling(max(max(d$nlp), BONF, 6)) + 1; setkey(d, SNP)

  K <- S[trait == tr]; setorder(K, chr, lead_bp)
  K[, colour := PALETTE[((seq_len(.N)-1L) %% length(PALETTE)) + 1L]]
  K[, lead_x := snp_map[lead_SNP, x]]
  # legend line carries the SPAN COORDINATES, not just the lead position
  K[, leg := sprintf("%-14s %s:%s-%s  |  %9s kb, %4d SNPs  |  lead %s",
       locus_id, chr, format(start, big.mark=","), format(end, big.mark=","),
       formatC(span_kb, format="f", digits=1, big.mark=","),
       n_members, format(lead_bp, big.mark=","))]
  MM <- merge(M[trait == tr], K[, .(locus_id, colour)], by="locus_id")
  MM[, `:=`(x = snp_map[SNP, x], nlp = d[SNP, nlp])]; MM <- MM[is.finite(nlp)]

  n <- nrow(K); NC <- if (n > 8) 2L else 1L; nr <- ceiling(n/NC)
  H_LEG <- ROW_IN*(nr+1L) + 0.45; H_IN <- H_MAN + H_LEG
  draw <- function() {
    layout(matrix(1:2, nrow=2), heights=c(H_MAN, H_LEG))
    par(mar=c(4.5,5.2,3,1.2), las=1, cex.axis=1.15, cex.lab=1.2, cex.main=1.05)
    plot(NA, xlim=XLIM, ylim=c(0,ymax), xaxs="i", yaxs="i", axes=FALSE,
         xlab="Chromosome", ylab=expression(-log[10](italic(p))),
         main=sprintf("%s | BLUP x 3 PCs | %d loci: --clump r2>=%.1f, +/-%.0f Mb, no p-value on members, severed at any gap > %.0f kb",
                      TRAIT_LAB[[tr]], n, P_R2, P_KB/1000, P_GAP))
    rasterImage(img, XLIM[1], 0, XLIM[2], ymax, interpolate=FALSE)
    axis(2); axis(1, at=cc$c, labels=paste0(cc$CHR,"H"), tick=FALSE); box()
    abline(h=BONF, lwd=1.4)
    if (nrow(MM)) points(MM$x, MM$nlp, pch=16, cex=1.05, col=MM$colour)
    W <- K[n_members > 1]; N <- K[n_members == 1]
    if (nrow(W)) points(W$lead_x, d[W$lead_SNP, nlp], pch=16, cex=1.6, col=W$colour)
    if (nrow(N)) points(N$lead_x, d[N$lead_SNP, nlp], pch=1, cex=1.7, lwd=2.2, col=N$colour)
    legend("topleft", bty="n", cex=1.0, pch=c(16,16,1), col="grey20",
           pt.cex=c(1.2,1.8,1.9), pt.lwd=c(1,1,2.4),
           legend=c("locus member (r2 >= 0.5, connected)","lead SNP","lead SNP, no member survived the gap rule"))
    par(mar=c(0.2,0.6,0.2,0.6)); plot.new(); plot.window(c(0,1),c(0,1))
    for (i in seq_len(n)) {
      col_i <- (i-1L)%/%nr; row_i <- (i-1L)%%nr
      xx <- col_i/NC + 0.004; yy <- 1 - (row_i+0.70)/(nr+1L)
      points(xx, yy, pch=if (K$n_members[i] > 1) 16 else 1, col=K$colour[i], cex=1.7, lwd=2.2, xpd=NA)
      text(xx+0.016, yy, K$leg[i], adj=c(0,0.5), cex=CEX_LEG, xpd=NA, family="mono")
    }
    text(0, 1-(nr+0.55)/(nr+1L),
      sprintf("Loci: --clump p1=%s (leads), p2=1 (membership on LD alone), r2>=%.1f, +/-%.0f Mb; then severed at the first gap > %.0f kb from the lead outward.", P_P1, P_R2, P_KB/1000, P_GAP),
      adj=c(0,0), cex=CEX_FOOT, col="grey35", xpd=NA)
  }
  tmp <- tempfile(fileext=".png", tmpdir=Sys.getenv("TMPDIR"))
  png(tmp, width=W_IN, height=H_MAN, units="in", res=RASTER_DPI, bg="transparent", type="cairo-png")
  par(mar=c(0,0,0,0))
  plot(d$x, d$nlp, xlim=XLIM, ylim=c(0,ymax), xaxs="i", yaxs="i", axes=FALSE,
       xlab="", ylab="", pch=16, cex=0.85, col=d$col)
  dev.off(); img <- readPNG(tmp); unlink(tmp)
  b <- sprintf("manhattan_%s__loci_r05", tr)
  cairo_pdf(file.path(OUT, paste0(b,".pdf")), width=W_IN, height=H_IN, pointsize=10); draw(); dev.off()
  png(file.path(OUT, paste0(b,".png")), width=W_IN, height=H_IN, units="in", res=PNG_DPI,
      pointsize=10, type="cairo-png"); draw(); dev.off()
  fwrite(K[, .(locus_id, lead_SNP, chr, lead_bp, span_kb, n_members, colour)],
         file.path(OUT, sprintf("locus_colour_key__%s.tsv", tr)), sep="\t")
  cat(sprintf("[33] %-11s %2d loci | %5d member SNPs painted | %d lead-only\n",
              tr, n, nrow(MM), sum(K$n_members==1)))
  rm(gw, d, img); gc(verbose=FALSE)
}
cat(sprintf("[33] OK -> %s\n", OUT))
