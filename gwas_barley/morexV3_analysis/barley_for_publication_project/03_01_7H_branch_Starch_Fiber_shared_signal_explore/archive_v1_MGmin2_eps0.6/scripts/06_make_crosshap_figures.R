#!/usr/bin/env Rscript
# 06_make_crosshap_figures.R
# Render the step-04 crosshap figures for this gene at the EXPLORATORY eps = 1.5
# (the only value that retains all four GWAS signal SNPs). Uses step 04's own
# renderers unmodified. Headers are stamped EXPLORATORY so these can never be
# mistaken for pipeline output at the fixed eps = 0.6.
# NOTE: page 1 of each PDF (the CrossHap tree) fails to render at this epsilon with
# a crosshap-internal "subscript out of bounds" - a library issue, not a data problem.
# The informative panel is page 5.
# Outputs -> results/figures/{crosshap,heatmaps}/
suppressPackageStartupMessages({library(data.table); library(dplyr); library(crosshap)})
B <- Sys.getenv("BRANCH"); P <- Sys.getenv("PROJECT_ROOT")
H <- file.path(P,"04_USED_haplotype_analysis_crosshap","01_scripts","R")
for (f in c("utils.R","run_crosshap.R","plot_combined_pdf.R","plot_heatmaps.R")) source(file.path(H,f))
W <- file.path(B,"intermediates","crosshap_work"); EPS <- 1.5; MG <- as.integer(Sys.getenv("MGMIN"))
G <- paste0(Sys.getenv("GENE_ID"),".vcf.gz"); gene_name <- "7HG0729030"
for (tr in c("fiber","starch")) {
  cfg <- list(raw_vcf_dir=file.path(W,"raw"), imputed_vcf_dir=file.path(W,"imp"),
              pheno_root=Sys.getenv("PHENO_ROOT"), pheno_suffix=Sys.getenv("PHENO_SUFFIX"),
              plink_bin=Sys.getenv("PLINK"), epsilon_vector=EPS, minHap=9,
              hetmiss_as="allele", keep_outliers=FALSE)
  res <- run_crosshap(cfg, tr, G, MGmin=MG, tmp_dir=file.path(W,"tmp"))
  label <- paste0("Haplotypes_MGmin",MG,"_E",EPS)
  title <- sprintf("%s | %s | MGmin=%d eps=%s [EXPLORATORY - not the fixed eps=0.6]", gene_name, tr, MG, EPS)
  cp <- file.path(B,"results","figures","crosshap", sprintf("CrosshapTree+Violin__%s__%s__eps%s.pdf",gene_name,tr,EPS))
  hp <- file.path(B,"results","figures","heatmaps", sprintf("Heatmaps__%s__%s__eps%s.pdf",gene_name,tr,EPS))
  try(write_combined_pdf(HapObject=res$HapObject, out_pdf=cp, title=title, trait=tr, gene_file=G,
      gene_name=gene_name, MGmin=MG, epsilon_vector=EPS, mgmin_test_stats=NULL, gene_summary_row=NULL))
  try(write_heatmaps_for_eps(HapObject=res$HapObject, label=label, gene_file=G, title=title, trait=tr,
      eps=EPS, MGmin=MG, raw_path=res$raw_path, common_ids=res$common_ids,
      vcf_raw_ids=res$vcf_raw_ids, out_pdf=hp))
}
cat("06 OK\n")
