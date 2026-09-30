#!/usr/bin/env Rscript
# 03_violin_plus_barcode.R -- EXPLORATORY: the step-07 "shared_sites" figure for
# 7HG0729030, at MGmin = 3.  Layout, palette, statistics and captions are ported
# from 07_.../scripts/04_figures.R so these read identically to the published ones:
#   TOP     violin + boxplot per wild haplotype group, x label "A\n(n=..)"
#           brackets = Wilcoxon of every group vs the LARGEST group, Holm-corrected
#   BOTTOM  aligned barcodes: wild group consensus rows, a gap, then elite lines
#   ADDED   red triangles under the 4 GWAS signal SNP columns (the crosshap marker
#           group MG1) -- the one deliberate departure from step 07, which draws no
#           marker-group annotation.
suppressPackageStartupMessages({library(vcfR); library(data.table); library(ggplot2); library(patchwork)})
T <- Sys.getenv("TEMP_ROOT"); P <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
E <- file.path(P,"07_USED_elite_lines_compariosn_to_wild_lines")
GENE <- "HORVU.MOREX.r3.7HG0729030"; SHORT <- "GDSL"
TITLE <- "Secreted GDSL esterase/lipase (SGNH hydrolase)"
CHR <- "7H"; GS <- 573604051; GE <- 573605728; WS <- 573603051; WE <- 573606728; WIN <- 1000
SIG <- c(573606282,573606306,573606460,573606491)
STATE_COLS <- c(Reference="#FFFACD", Alternate="#2F4F4F", Missing="grey70", Heterozygous="#C46210")
elite_cfg <- fread(cmd=paste("grep -v '^#'", shQuote(file.path(E,"config/elite_lines.tsv"))))
gt_state <- function(gt){o <- rep("Missing", length(gt))
  o[gt %in% c("0/0","0|0")] <- "Reference"; o[gt %in% c("1/1","1|1")] <- "Alternate"
  o[gt %in% c("0/1","1/0","0|1","1|0")] <- "Heterozygous"; o}
read_gt <- function(p){ v <- suppressMessages(read.vcfR(p, verbose=FALSE))
  fx <- as.data.frame(getFIX(v), stringsAsFactors=FALSE); g <- extract.gt(v,"GT",as.numeric=FALSE)
  k <- nchar(fx$REF)==1 & nchar(fx$ALT)==1
  list(key=paste(sub("^chr","",fx$CHROM[k]),fx$POS[k],fx$REF[k],fx$ALT[k],sep=":"),
       pk=paste(sub("^chr","",fx$CHROM[k]),fx$POS[k],sep=":"),
       pos=as.integer(fx$POS[k]), gt=g[k,,drop=FALSE]) }
wild <- read_gt(file.path(T,"work/raw/fiber",paste0(GENE,".vcf.gz")))
el   <- read_gt(file.path(T,"work/elite/7HG0729030.elite.vcf.gz"))
shared <- intersect(wild$key, el$key)
shared <- shared[order(as.integer(sub("^[^:]+:([0-9]+):.*","\\1",shared)))]
triallelic <- length(intersect(wild$pk, el$pk)) - length(shared)
n_wild_only <- length(wild$key) - length(shared) - triallelic
n_elite_only <- length(el$key) - length(shared) - triallelic
spos <- as.integer(sub("^[^:]+:([0-9]+):.*","\\1", shared))
theme_pub <- function(base=11) theme_bw(base_size=base, base_family="sans") +
  theme(panel.grid=element_blank(), panel.border=element_rect(colour="grey40",fill=NA,linewidth=0.4),
        axis.text=element_text(colour="black"), plot.title=element_blank(), plot.subtitle=element_blank())
p_sym <- function(p) if(is.na(p)) "ns" else if(p<=1e-4) "****" else if(p<=1e-3) "***" else if(p<=0.01) "**" else if(p<=0.05) "*" else "ns"
pairwise_vs_largest <- function(d){ gsv <- sort(unique(d$hap)); if(length(gsv)<=1) return(NULL)
  sz <- table(d$hap); ref <- names(sz)[which.max(sz)]; oth <- setdiff(gsv, ref); if(!length(oth)) return(NULL)
  raw <- vapply(oth, function(g) tryCatch(wilcox.test(d$Pheno[d$hap==ref], d$Pheno[d$hap==g], exact=FALSE)$p.value,
         error=function(e) NA_real_), numeric(1))
  adj <- p.adjust(raw,"holm"); ymax <- max(d$Pheno,na.rm=TRUE); ysp <- max(1e-8, ymax-min(d$Pheno,na.rm=TRUE))
  data.frame(g1=ref, g2=oth, p=raw, p.adj=adj, sym=vapply(adj,p_sym,character(1)),
             y=ymax+ysp*(0.10+seq_along(oth)*0.08), stringsAsFactors=FALSE) }
pw_out <- list()
for (tr in c("fiber","starch")) for (EPS in c(0.9)) {   # epsilon 0.6 dropped 2026-09-24 (user)
  ho <- readRDS(file.path(T,"Cache",tr,GENE,"MGmin_3",paste0("eps_",EPS),"HapObject.rds"))
  ind <- as.data.table(ho$HapObject[[paste0("Haplotypes_MGmin3_E",EPS)]]$Indfile)
  ind[, `:=`(Ind=as.character(Ind), hap=as.character(hap), Pheno=as.numeric(Pheno))]
  ind <- ind[hap!="0" & !is.na(Pheno)]
  lv <- sort(unique(ind$hap)); ind[, hap := factor(hap, levels=lv)]
  gstat <- ind[, .(n=.N), by=hap][order(hap)]
  # ---- TOP: violin
  fillv <- setNames(grDevices::hcl.colors(length(lv), palette="Dark 3")[seq_along(lv)], lv)
  p_v <- ggplot(ind, aes(hap, Pheno, fill=hap)) +
    geom_violin(trim=FALSE, colour="black", linewidth=0.35, width=0.9) +
    geom_boxplot(width=0.14, fill="white", colour="black", outlier.size=0.6, linewidth=0.35) +
    scale_fill_manual(values=fillv, guide="none") +
    scale_x_discrete(labels=setNames(sprintf("%s\n(n=%d)", gstat$hap, gstat$n), gstat$hap)) +
    labs(x="Wild haplotype group", y=paste0(tr," BLUP")) + theme_pub()
  pw <- pairwise_vs_largest(ind)
  if (!is.null(pw)) { tip <- diff(range(ind$Pheno,na.rm=TRUE))*0.012
    br <- data.frame(x=match(pw$g1,lv), xend=match(pw$g2,lv), y=pw$y, lab=pw$sym)
    p_v <- p_v +
      geom_segment(data=br, aes(x=x,xend=xend,y=y,yend=y), inherit.aes=FALSE, linewidth=0.3) +
      geom_segment(data=br, aes(x=x,xend=x,y=y-tip,yend=y), inherit.aes=FALSE, linewidth=0.3) +
      geom_segment(data=br, aes(x=xend,xend=xend,y=y-tip,yend=y), inherit.aes=FALSE, linewidth=0.3) +
      geom_text(data=br, aes(x=(x+xend)/2, y=y+tip*0.6, label=lab), inherit.aes=FALSE, size=3.2, vjust=0) +
      expand_limits(y=max(pw$y)+tip*4)
    pw_out[[length(pw_out)+1]] <- cbind(trait=tr, eps=EPS, test="Wilcoxon vs largest group, Holm", pw) }
  # ---- BOTTOM: barcodes
  wm <- wild$gt[match(shared, wild$key), , drop=FALSE]
  cons <- sapply(lv, function(h){ sub <- wm[, colnames(wm) %in% ind$Ind[ind$hap==h], drop=FALSE]
    apply(sub, 1, function(r){ tb <- table(gt_state(r)[gt_state(r)!="Missing"])
      if(!length(tb)) return("Missing"); top <- names(tb)[tb==max(tb)]
      # exact tie -> Reference if it is one of the tied states, else Missing: step 07's rule
      # since 2026-09-24 (aligned here 2026-09-27, user approval; was: all ties -> Missing;
      # this gene has no tie, so no output changed)
      if(length(top)>1) return(if("Reference" %in% top) "Reference" else "Missing"); top }) })
  em <- el$gt[match(shared, el$key), , drop=FALSE]
  em <- em[, match(elite_cfg$sample_id, colnames(em)), drop=FALSE]; colnames(em) <- elite_cfg$line_name
  bar <- rbind(
    data.table(site=rep(shared, times=length(lv)), label=rep(paste0("Group ",lv), each=length(shared)),
               state=as.vector(cons), block="wild"),
    data.table(site=rep(shared, times=ncol(em)), label=rep(colnames(em), each=length(shared)),
               state=as.vector(apply(em,2,gt_state)), block="elite"))
  SPACER <- "​"
  bar[, site := factor(site, levels=shared)]
  bar[, label := factor(label, levels=c(rev(colnames(em)), SPACER, rev(paste0("Group ",lv))))]
  sigdf <- data.frame(site=factor(shared[spos %in% SIG], levels=shared))
  p_b <- ggplot(bar, aes(site, label, fill=state)) +
    geom_tile(colour="grey78", linewidth=0.25, height=0.72) +
    geom_point(data=sigdf, aes(x=site, y=0.34), inherit.aes=FALSE, shape=17, colour="red", size=2.3) +
    scale_fill_manual(values=STATE_COLS, drop=TRUE, name=NULL) +
    guides(fill=guide_legend(override.aes=list(colour="grey40"))) +
    scale_y_discrete(drop=FALSE, expand=expansion(add=c(1.0,0.6))) +
    labs(x=sprintf("%d SNPs in %s:%s-%s (gene %s-%s, ±%d bp), genomic order\n▲ = crosshap marker group MG1 — the 4 GWAS signal SNPs",
         length(shared), CHR, format(WS,big.mark=","), format(WE,big.mark=","),
         format(GS,big.mark=","), format(GE,big.mark=","), WIN), y=NULL) +
    theme_pub() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.position="bottom", legend.key.size=unit(0.4,"cm"), axis.text.y=element_text(size=9))
  header <- sprintf("%s  |  %s  |  %s  (%s)", SHORT, TITLE, GENE, tr)
  sub <- sprintf("wild %d SNPs, elite %d records in window; shared %d, triallelic %d, wild-only %d, elite-only %d (not drawn); %d columns drawn\nMGmin = 3, eps = %s   [EXPLORATORY \u2014 the pipeline uses MGmin = 2, under which this gene has no usable grouping]",
    length(wild$key), length(el$key), length(shared), triallelic, n_wild_only, n_elite_only, length(shared), EPS)
  n_rows <- length(lv) + ncol(em) + 1
  fig <- (p_v / p_b) + plot_layout(heights=c(2.1, max(1, n_rows*0.22))) +
    plot_annotation(title=header, subtitle=sub,
      theme=theme(plot.title=element_text(face="bold", size=12),
                  plot.subtitle=element_text(size=8.2, colour="grey25")))
  base <- file.path(T,"results/figures/shared_sites", sprintf("%s__%s__%s__MGmin3_eps%s", tr, SHORT, GENE, EPS))
  h <- 8.9 + max(0, n_rows-8)*0.16
  ggsave(paste0(base,".pdf"), fig, width=12, height=h, device=grDevices::cairo_pdf)
  ggsave(paste0(base,".png"), fig, width=12, height=h, dpi=400)
  cat(sprintf("  wrote %s (eps %s)\n", tr, EPS))
}
fwrite(rbindlist(pw_out), file.path(T,"results/tables/pairwise_group_tests.tsv"), sep="\t")
cat(sprintf("\nshared %d | triallelic %d | wild-only %d | elite-only %d\n", length(shared), triallelic, n_wild_only, n_elite_only))
