#!/usr/bin/env Rscript
# 05_crosshap_and_epsilon_sweep.R
# Run crosshap on 7HG0729030 (gene +/- 1000 bp) exactly as step 04 would, then sweep
# epsilon to diagnose WHY the gene is not recovered.
#
# RESULT AT THE PIPELINE'S FIXED eps = 0.6: NOT significant (fiber p=0.16, starch p=0.36).
# REASON: DBSCAN assigns ALL FOUR GWAS signal SNPs to marker group 0 (noise). The only
# marker group that forms is built from two NULL SNPs at the far end of the window.
# The KW test therefore tests a partition unrelated to the signal.
# This is the documented epsilon trap: eps is a Euclidean radius over each SNP's full
# r2 PROFILE, not an r2 cutoff, so its stringency depends on SNP density in the window.
#
# *** THE SWEEP IS DIAGNOSTIC ONLY - NOT A RESULT. ***
# Choosing epsilon by outcome is the forking path the main pipeline deliberately closed
# by fixing eps = 0.6. Nothing from the sweep may be reported as a finding.
# Outputs -> results/tables/{crosshap_epsilon_sweep.tsv, crosshap_marker_groups.tsv}
suppressPackageStartupMessages({library(data.table); library(dplyr)})
B <- Sys.getenv("BRANCH"); P <- Sys.getenv("PROJECT_ROOT"); OUT <- file.path(B,"results","tables")
H <- file.path(P,"04_USED_haplotype_analysis_crosshap","01_scripts","R")
for (f in c("utils.R","run_crosshap.R")) source(file.path(H,f))   # read-only reuse of step-04 code
W <- file.path(B,"intermediates","crosshap_work"); dir.create(file.path(W,"tmp"), recursive=TRUE, showWarnings=FALSE)
G <- paste0(Sys.getenv("GENE_ID"),".vcf.gz")
for (tr in c("fiber","starch")) for (k in c("raw","imp")) {
  dir.create(file.path(W,k,tr), recursive=TRUE, showWarnings=FALSE)
  file.copy(file.path(B,"intermediates", if(k=="raw") "gene_window_raw.vcf.gz" else "gene_window_imputed.vcf.gz"),
            file.path(W,k,tr,G), overwrite=TRUE)
  file.copy(file.path(B,"intermediates", if(k=="raw") "gene_window_raw.vcf.gz.csi" else "gene_window_imputed.vcf.gz.csi"),
            file.path(W,k,tr,paste0(G,".csi")), overwrite=TRUE)
}
SIG <- paste0("7H:", c(573606282,573606306,573606460,573606491))
MG <- as.integer(Sys.getenv("MGMIN")); EPSV <- as.numeric(strsplit(Sys.getenv("EPS_EXPLORE")," ")[[1]])
sweep <- list(); MGCOLLECT <- list()
for (tr in c("fiber","starch")) for (e in EPSV) {
  cfg <- list(raw_vcf_dir=file.path(W,"raw"), imputed_vcf_dir=file.path(W,"imp"),
              pheno_root=Sys.getenv("PHENO_ROOT"), pheno_suffix=Sys.getenv("PHENO_SUFFIX"),
              plink_bin=Sys.getenv("PLINK"), epsilon_vector=e, minHap=9,
              hetmiss_as="allele", keep_outliers=FALSE)
  out <- tryCatch({
    r <- run_crosshap(cfg, tr, G, MGmin=MG, tmp_dir=file.path(W,"tmp"))
    lab <- names(r$HapObject)[1]; HO <- r$HapObject[[lab]]
    v <- as.data.table(HO$Varfile); ind <- HO$Indfile
    if (e %in% c(0.6, 1.5)) { vv <- copy(v); vv[, `:=`(trait=tr, eps=e,
        is_gwas_signal = ID %in% SIG)]
      assign("MGCOLLECT", c(get("MGCOLLECT", envir=.GlobalEnv),
             list(vv[, .(trait,eps,ID,MGs,is_gwas_signal)])), envir=.GlobalEnv) }
    d <- data.table(hap=as.character(ind$hap), Pheno=suppressWarnings(as.numeric(ind$Pheno)))
    nun <- sum(d$hap=="0", na.rm=TRUE); d <- d[hap!="0" & !is.na(Pheno)]; k <- uniqueN(d$hap)
    kt <- if (k>=2 && nrow(d)>0) kruskal.test(Pheno~hap, data=d) else NULL
    n <- nrow(d); Hs <- if(!is.null(kt)) unname(kt$statistic) else NA_real_
    data.table(trait=tr, epsilon=e, marker_groups=uniqueN(v[MGs!="0"]$MGs),
      gwas_signal_snps_kept=sum(v[ID %in% SIG]$MGs!="0"), haplotype_groups=k,
      n_assigned=n, n_unassigned=nun, kw_p=if(!is.null(kt)) kt$p.value else NA_real_,
      eta_squared=if(!is.null(kt) && n>k) max(0,min(1,(Hs-k+1)/(n-k))) else NA_real_, status="ok")
  }, error=function(err) data.table(trait=tr, epsilon=e, marker_groups=NA_integer_,
      gwas_signal_snps_kept=NA_integer_, haplotype_groups=NA_integer_, n_assigned=NA_integer_,
      n_unassigned=NA_integer_, kw_p=NA_real_, eta_squared=NA_real_,
      status=paste("crosshap_error:", substr(conditionMessage(err),1,60))))
  sweep[[length(sweep)+1]] <- out
}
S <- rbindlist(sweep); fwrite(S, file.path(OUT,"crosshap_epsilon_sweep.tsv"), sep="\t")
fwrite(rbindlist(MGCOLLECT), file.path(OUT,"crosshap_marker_groups.tsv"), sep="\t")
print(S); cat("05 OK\n")
