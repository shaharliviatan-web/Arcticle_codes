#!/usr/bin/env Rscript
# =============================================================================
# param_grid_scan.R -- DIAGNOSTIC ONLY. Never part of the pipeline.
#
# WHAT   Scan a grid of MGmin x epsilon and report, per config, how many genes
#        yield a USABLE haplotype grouping (>= 2 groups) and how well accessions
#        are assigned. These are GENOTYPE-ONLY measures -- crosshap never uses the
#        phenotype to build groups, so choosing a config on them is legitimate.
#
#        p-values are NOT computed here, deliberately. Choosing a config because it
#        maximises significance is the forking path that step 04 was rewritten to
#        remove; this script cannot be used that way.
#
# SPEED  The expensive work per gene -- reading the raw+imputed VCFs and running
#        PLINK for the LD matrix -- does not depend on MGmin or epsilon, so it is
#        done ONCE per gene and every config is then evaluated against the cached
#        LD. Each config is wrapped in its own tryCatch, because crosshap aborts the
#        whole call if any single epsilon in a vector fails (the bug that made the
#        earlier vector-sweep under-count coverage).
#
# WRITES 04_runs/_diagnostic_param_grid/grid_by_config.tsv   one row per config
#        04_runs/_diagnostic_param_grid/grid_by_gene.tsv     one row per gene x config
#        04_runs/_diagnostic_param_grid/grid_by_trait.tsv    one row per trait x config
# =============================================================================
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table); library(crosshap); library(dplyr) })

sd_ <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/04_USED_haplotype_analysis_crosshap/01_scripts"
source(file.path(sd_,"R","utils.R")); source(file.path(sd_,"R","run_crosshap.R"))
cfg <- yaml::read_yaml(file.path(dirname(sd_),"00_config","config.yaml"))

MGMINS <- c(2L, 3L, 4L)
EPSS   <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0)

gw  <- as.data.table(read_gene_windows(file.path(dirname(sd_),"00_config","gene_windows.tsv")))
out <- file.path(cfg$output_root, "04_runs", "_diagnostic_param_grid"); ensure_dir(out)
tmp <- file.path(out, "tmp"); ensure_dir(tmp)

rows <- list()
for (i in seq_len(nrow(gw))) {
  tr <- gw$trait[i]; gf <- gw$gene_file[i]; gn <- make_gene_name(gf)
  cat(sprintf("[%2d/%d] %-11s %s\n", i, nrow(gw), tr, gn)); flush.console()

  # ---- expensive setup, once per gene: VCFs + PLINK LD -----------------------
  # Reuse run_crosshap with a single trivial epsilon purely to obtain vcf/LD/pheno.
  cfg1 <- cfg; cfg1$epsilon_vector <- 0.6
  base <- tryCatch(run_crosshap(cfg=cfg1, trait=tr, gene_file=gf, MGmin=2L, tmp_dir=tmp),
                   error=function(e) NULL)
  if (is.null(base)) {
    for (M in MGMINS) for (E in EPSS)
      rows[[length(rows)+1]] <- data.table(trait=tr, gene_id=gn, MGmin=M, eps=E,
        n_snps=NA_integer_, ok=FALSE, n_groups=NA_integer_, assign_rate=NA_real_,
        note="setup_failed_no_variants_or_LD")
    next
  }
  n_snps <- base$n_common

  # Rebuild the inputs run_haplotyping needs. Cheap: no PLINK, no file IO.
  raw_pheno <- data.table::fread(pheno_path(cfg, tr), header=FALSE)
  colnames(raw_pheno) <- c("FID","IID","trait")
  pheno <- raw_pheno %>% dplyr::mutate(Ind=as.character(IID)) %>%
    dplyr::select(Ind, Pheno=trait) %>% dplyr::distinct(Ind, .keep_all=TRUE)
  vcf_raw <- read_vcf_robust(base$raw_path)
  vcf_raw <- vcf_raw[paste0(vcf_raw[["#CHROM"]],":",vcf_raw[["POS"]]) %in% base$common_ids, ]
  vcf_raw$ID <- base$vcf_raw_ids

  for (M in MGMINS) {
    for (E in EPSS) {
      r <- tryCatch(suppressWarnings(suppressMessages(run_haplotyping(
             vcf=vcf_raw, LD=base$LD, pheno=pheno, epsilon=E, MGmin=M,
             minHap=cfg$minHap, hetmiss_as=cfg$hetmiss_as,
             keep_outliers=cfg$keep_outliers))), error=function(e) NULL)
      ng <- NA_integer_; ar <- NA_real_; ok <- FALSE; note <- "no_marker_groups"
      if (!is.null(r)) {
        H <- r[[paste0("Haplotypes_MGmin",M,"_E",E)]]
        if (!is.null(H) && !is.null(H$Indfile) && nrow(H$Indfile)) {
          hp <- as.character(H$Indfile$hap)
          ph <- suppressWarnings(as.numeric(H$Indfile$Pheno))
          k  <- hp!="0" & !is.na(ph)
          ng <- length(unique(hp[k])); ar <- mean(hp!="0")
          ok <- ng >= 2; note <- if (ok) "ok" else "fewer_than_2_groups"
        }
      }
      rows[[length(rows)+1]] <- data.table(trait=tr, gene_id=gn, MGmin=M, eps=E,
        n_snps=n_snps, ok=ok, n_groups=ng, assign_rate=ar, note=note)
    }
  }
}

d <- rbindlist(rows)
fwrite(d, file.path(out,"grid_by_gene.tsv"), sep="\t", na="NA")

by_cfg <- d[, .(genes_testable=sum(ok), pct_testable=round(100*mean(ok)),
                median_groups=as.numeric(median(n_groups[ok], na.rm=TRUE)),
                median_assign_pct=round(100*median(assign_rate[ok], na.rm=TRUE))),
            by=.(MGmin, eps)][order(-genes_testable)]
fwrite(by_cfg, file.path(out,"grid_by_config.tsv"), sep="\t", na="NA")

by_tr <- d[, .(genes=uniqueN(gene_id), testable=sum(ok)), by=.(trait, MGmin, eps)]
fwrite(by_tr, file.path(out,"grid_by_trait.tsv"), sep="\t", na="NA")

cat("\n=== testable genes per config (of", uniqueN(d$gene_id), ") ===\n")
print(head(by_cfg, 15))
cat("\n=== PROTEIN only ===\n")
print(by_tr[trait=="protein"][order(-testable)][1:8])
