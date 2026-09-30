#!/usr/bin/env Rscript
# 01_run_mgmin3.R  -- 7HG0729030 at MGmin = 3, epsilon 0.9.
# 2026-09-24 (user decision): epsilon 0.6 dropped from the run; only the chosen epsilon 0.9 is
# kept. The comparison that chose 0.9 over 0.6 (larger assignment) stays in
# results/tables/param_grid_eps0.05-1.5_MGmin2-3.tsv and README section 2.
#
# WHY: at the pipeline's MGmin = 2 this gene is null, because DBSCAN assigns all four
# GWAS signal SNPs to marker group 0 (noise) and builds haplotypes from two NULL SNPs.
# At MGmin = 3 the four correlated signal SNPs form a proper marker group, 85% of the
# panel is assigned, and the minority haplotype is EXACTLY the 34 accessions that are
# ALT at both lead SNPs in the direct genotype split.
#
# Uses step 04's run_crosshap.R and its plot renderers UNMODIFIED. Writes HapObject
# caches in step 04's cache layout so the elite-line scripts can consume them.
# NOTHING here is part of the pipeline.
suppressPackageStartupMessages({library(data.table); library(dplyr); library(crosshap)})
P <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
B <- file.path(P,"04_USED_haplotype_analysis_crosshap")
for (f in c("utils.R","run_crosshap.R","plot_combined_pdf.R","plot_heatmaps.R")) source(file.path(B,"01_scripts/R",f))
T <- Sys.getenv("TEMP_ROOT"); stopifnot(nzchar(T))
W <- file.path(T,"work"); GENE <- "HORVU.MOREX.r3.7HG0729030"; G <- paste0(GENE,".vcf.gz")
MG <- 3L; EPSL <- c(0.9);   # was c(0.6, 0.9) until 2026-09-24
 SIG <- paste0("7H:", c(573606282,573606306,573606460,573606491))
stats <- list(); groups <- list(); assign <- list()
for (tr in c("fiber","starch")) for (EPS in EPSL) {
  cfg <- list(raw_vcf_dir=file.path(W,"raw"), imputed_vcf_dir=file.path(W,"imp"),
              pheno_root="/mnt/data/shahar/gwas_barley/data/inputs", pheno_suffix="_corrected_V3.pheno",
              plink_bin="/usr/local/bin/plink", epsilon_vector=EPS, minHap=9,
              hetmiss_as="allele", keep_outliers=FALSE)
  res <- run_crosshap(cfg, tr, G, MGmin=MG, tmp_dir=file.path(W,"tmp"))
  # cache in step-04 layout: Cache/<trait>/<gene>/MGmin_<N>/HapObject.rds
  cdir <- file.path(T,"Cache",tr,GENE,paste0("MGmin_",MG),paste0("eps_",EPS))
  dir.create(cdir, recursive=TRUE, showWarnings=FALSE); saveRDS(res, file.path(cdir,"HapObject.rds"))
  lab <- paste0("Haplotypes_MGmin",MG,"_E",EPS); HO <- res$HapObject[[lab]]
  v <- as.data.table(HO$Varfile); ind <- as.data.table(HO$Indfile)
  ind[, `:=`(Ind=as.character(Ind), hap=as.character(hap))]
  assign[[length(assign)+1]] <- ind[, .(trait=tr, epsilon=EPS, MGmin=MG, Ind, hap, Pheno)]
  d <- ind[hap!="0" & !is.na(as.numeric(Pheno))][, Pheno := as.numeric(Pheno)]
  k <- uniqueN(d$hap); kt <- kruskal.test(Pheno~hap, data=d); n <- nrow(d); H <- unname(kt$statistic)
  g <- d[, .(n=.N, mean=mean(Pheno), median=median(Pheno), sd=sd(Pheno)), by=hap][order(-mean)]
  groups[[length(groups)+1]] <- cbind(trait=tr, epsilon=EPS, MGmin=MG, g)
  stats[[length(stats)+1]] <- data.table(trait=tr, MGmin=MG, epsilon=EPS,
    marker_groups=uniqueN(v[MGs!="0"]$MGs), gwas_signal_snps_kept=sum(v[ID %in% SIG]$MGs!="0"),
    n_snps_window=nrow(v), haplotype_groups=k, n_assigned=n, n_unassigned=sum(ind$hap=="0"),
    group_sizes=paste(sort(g$n,decreasing=TRUE),collapse="|"), kw_H=H, kw_df=unname(kt$parameter),
    kw_p_raw=kt$p.value, eta_squared=max(0,min(1,(H-k+1)/(n-k))),
    delta_top_bottom_sd=(g$mean[1]-g$mean[nrow(g)])/sd(d$Pheno))
  ttl <- sprintf("7HG0729030 | %s | MGmin=%d eps=%s [EXPLORATORY - pipeline uses MGmin=2 eps=0.6]", tr, MG, EPS)
  try(write_combined_pdf(HapObject=res$HapObject,
      out_pdf=file.path(T,"results/figures/crosshap",sprintf("CrosshapTree+Violin__7HG0729030__%s__MGmin%d_eps%s.pdf",tr,MG,EPS)),
      title=ttl, trait=tr, gene_file=G, gene_name="7HG0729030", MGmin=MG, epsilon_vector=EPS,
      mgmin_test_stats=NULL, gene_summary_row=NULL))
  try(write_heatmaps_for_eps(HapObject=res$HapObject, label=lab, gene_file=G, title=ttl, trait=tr,
      eps=EPS, MGmin=MG, raw_path=res$raw_path, common_ids=res$common_ids, vcf_raw_ids=res$vcf_raw_ids,
      out_pdf=file.path(T,"results/figures/heatmaps",sprintf("Heatmaps__7HG0729030__%s__MGmin%d_eps%s.pdf",tr,MG,EPS))))
  cat(sprintf("  done %s eps=%s\n", tr, EPS))
}
O <- file.path(T,"results","tables")
fwrite(rbindlist(stats),  file.path(O,"mgmin3_gene_results.tsv"), sep="\t")
fwrite(rbindlist(groups), file.path(O,"mgmin3_haplotype_groups.tsv"), sep="\t")
fwrite(rbindlist(assign), file.path(O,"mgmin3_haplotype_assignment.tsv"), sep="\t")
print(rbindlist(stats)[, .(trait,epsilon,marker_groups,gwas_signal_snps_kept,haplotype_groups,
     n_assigned,group_sizes,kw_p_raw=signif(kw_p_raw,3),eta2=round(eta_squared,3),delta=round(delta_top_bottom_sd,2))])
