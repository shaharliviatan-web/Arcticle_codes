#!/usr/bin/env Rscript
# eps_coverage_scan.R -- DIAGNOSTIC ONLY, never part of the pipeline.
# For every gene in gene_windows.tsv, run crosshap across a range of epsilon and
# record, per epsilon, whether a valid haplotype grouping exists (genotype-only:
# no p-values are used or reported here). Answers "which fixed epsilon maximises
# the number of testable genes on THIS gene set".
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table) })
sd_ <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/04_USED_haplotype_analysis_crosshap/01_scripts"
source(file.path(sd_,"R","utils.R")); source(file.path(sd_,"R","run_crosshap.R"))
cfg <- yaml::read_yaml(file.path(dirname(sd_),"00_config","config.yaml"))
cfg$epsilon_vector <- c(0.05,0.2,0.4,0.5,0.6,0.8,0.85,1.0,1.5,2.0,3.0)
gw <- as.data.table(read_gene_windows(file.path(dirname(sd_),"00_config","gene_windows.tsv")))
tmp <- file.path(cfg$output_root,"04_runs","_diagnostic_eps_scan","tmp"); ensure_dir(tmp)
out <- list()
for (i in seq_len(nrow(gw))) {
  tr <- gw$trait[i]; gf <- gw$gene_file[i]; gn <- make_gene_name(gf)
  cat(sprintf("[%2d/%d] %s %s\n", i, nrow(gw), tr, gn))
  r <- tryCatch(run_crosshap(cfg=cfg, trait=tr, gene_file=gf, MGmin=2L, tmp_dir=tmp),
                error=function(e) NULL)
  if (is.null(r)) { out[[length(out)+1]] <- data.table(trait=tr, gene_id=gn, n_snps=NA_integer_,
      eps=NA_real_, ok=FALSE, n_groups=NA_integer_); next }
  for (e in cfg$epsilon_vector) {
    lab <- paste0("Haplotypes_MGmin2_E", e); H <- r$HapObject[[lab]]
    ng <- NA_integer_; ok <- FALSE
    if (!is.null(H) && !is.null(H$Indfile) && nrow(H$Indfile)) {
      hp <- as.character(H$Indfile$hap); ph <- suppressWarnings(as.numeric(H$Indfile$Pheno))
      k <- hp!="0" & !is.na(ph); ng <- length(unique(hp[k])); ok <- ng>=2
    }
    out[[length(out)+1]] <- data.table(trait=tr, gene_id=gn, n_snps=r$n_common,
                                       eps=e, ok=ok, n_groups=ng)
  }
}
d <- rbindlist(out)
fwrite(d, file.path(cfg$output_root,"04_runs","_diagnostic_eps_scan","eps_coverage.tsv"), sep="\t")
cat("\n=== testable genes per fixed epsilon (of", uniqueN(d$gene_id), "genes) ===\n")
print(d[!is.na(eps), .(testable=sum(ok)), by=eps][order(eps)])
