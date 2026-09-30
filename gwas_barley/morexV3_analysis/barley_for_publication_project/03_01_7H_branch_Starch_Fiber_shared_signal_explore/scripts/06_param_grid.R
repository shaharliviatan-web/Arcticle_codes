#!/usr/bin/env Rscript
# 06_param_grid.R -- the MGmin x epsilon grid behind the choice of MGmin = 3, epsilon = 0.9
# for 7HG0729030 (README section 2). Added 2026-09-30 (user-approved): the grid table
# results/tables/param_grid_eps0.05-1.5_MGmin2-3.tsv existed since 2026-09-24 but no script in
# the project wrote it; it is an Online Resource of the manuscript (M&M, the 7H region), so it
# must be reproducible. This script regenerates it with the same crosshap wrapper and settings
# as 01_run_mgmin3.R (step 04's run_crosshap.R, unmodified; minHap 9, hetmiss "allele",
# keep_outliers FALSE; gene +/- 1 kb; LD from the imputed window, haplotypes from the raw one).
#
# GRID   MGmin 2 and 3 x epsilon 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
#        1.2, 1.5, for fiber and starch. One crosshap run per cell (a vector of epsilons would
#        abort the whole run if one value fails; see 04_.../01_scripts/README.md).
#
# SELECTION (README section 2) reads only the genotype columns: sig_kept (how many of the four
#        GWAS signal SNPs are in a marker group) and assigned (accessions in a haplotype group).
#        kw_p, eta2 and delta are recorded for completeness and were NOT used to choose the
#        parameters.
#
# COLUMNS trait, MGmin, eps, mgroups (marker groups), sig_kept, hgroups (haplotype groups),
#        assigned, unassigned, sizes (group sizes, largest first), kw_p (Kruskal-Wallis, raw),
#        eta2 ((H - k + 1) / (n - k), bounded to [0, 1]), delta ((highest - lowest group mean) /
#        SD of the tested values), status ("ok" or "ERR:<crosshap message>").
#
# RUN    source scripts/config.sh; Rscript scripts/06_param_grid.R        (~1-2 min)
#        OUT_FILE=<path> writes elsewhere (used to compare with the existing table before
#        replacing it). Not part of run_all.sh: it changes no result.
suppressPackageStartupMessages({library(data.table); library(dplyr); library(crosshap)})
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
P <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
B <- file.path(P, "04_USED_haplotype_analysis_crosshap")
for (f in c("utils.R", "run_crosshap.R")) source(file.path(B, "01_scripts/R", f))   # read-only reuse
T <- Sys.getenv("TEMP_ROOT"); stopifnot(nzchar(T))
W <- file.path(T, "work"); GENE <- Sys.getenv("GENE_ID", "HORVU.MOREX.r3.7HG0729030")
G <- paste0(GENE, ".vcf.gz")
SIG <- paste0("7H:", strsplit(Sys.getenv("SIGNAL_SNPS", "573606282 573606306 573606460 573606491"), " ")[[1]])
OUT_FILE <- Sys.getenv("OUT_FILE", file.path(T, "results", "tables", "param_grid_eps0.05-1.5_MGmin2-3.tsv"))
TMP <- file.path(W, "tmp", paste0("grid_", Sys.getpid())); dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(TMP, recursive = TRUE), add = TRUE)

MGMINS <- c(2L, 3L)
EPSV   <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5)

one_cell <- function(tr, mg, e) {
  cfg <- list(raw_vcf_dir = file.path(W, "raw"), imputed_vcf_dir = file.path(W, "imp"),
              pheno_root = "/mnt/data/shahar/gwas_barley/data/inputs", pheno_suffix = "_corrected_V3.pheno",
              plink_bin = "/usr/local/bin/plink", epsilon_vector = e, minHap = 9,
              hetmiss_as = "allele", keep_outliers = FALSE)
  tryCatch({
    r  <- run_crosshap(cfg, tr, G, MGmin = mg, tmp_dir = TMP)
    HO <- r$HapObject[[names(r$HapObject)[1]]]
    v  <- as.data.table(HO$Varfile); ind <- as.data.table(HO$Indfile)
    ind[, `:=`(hap = as.character(hap), Pheno = suppressWarnings(as.numeric(Pheno)))]
    d  <- ind[hap != "0" & !is.na(Pheno)]; k <- uniqueN(d$hap); n <- nrow(d)
    g  <- d[, .(n = .N, mean = mean(Pheno)), by = hap][order(-mean)]
    kt <- if (k >= 2) kruskal.test(Pheno ~ hap, data = d) else NULL
    H  <- if (!is.null(kt)) unname(kt$statistic) else NA_real_
    data.table(trait = tr, MGmin = mg, eps = e,
      mgroups = uniqueN(v[MGs != "0"]$MGs), sig_kept = sum(v[ID %in% SIG]$MGs != "0"),
      hgroups = k, assigned = n, unassigned = sum(ind$hap == "0"),
      sizes = paste(sort(g$n, decreasing = TRUE), collapse = "|"),
      kw_p = if (!is.null(kt)) kt$p.value else NA_real_,
      eta2 = if (!is.null(kt) && n > k) max(0, min(1, (H - k + 1) / (n - k))) else NA_real_,
      delta = if (k >= 2) (g$mean[1] - g$mean[nrow(g)]) / sd(d$Pheno) else NA_real_,
      status = "ok")
  }, error = function(err) data.table(trait = tr, MGmin = mg, eps = e,
      mgroups = NA_integer_, sig_kept = NA_integer_, hgroups = NA_integer_, assigned = NA_integer_,
      unassigned = NA_integer_, sizes = NA_character_, kw_p = NA_real_, eta2 = NA_real_,
      delta = NA_real_, status = paste0("ERR:", substr(conditionMessage(err), 1, 26))))   # 26 chars, as in the 2026-09-24 table
}

res <- rbindlist(lapply(c("fiber", "starch"), function(tr)
  rbindlist(lapply(MGMINS, function(mg) rbindlist(lapply(EPSV, function(e) {
    x <- one_cell(tr, mg, e)
    cat(sprintf("  %-6s MGmin %d eps %-4s -> %s\n", tr, mg, e, x$status)); x }))))))
fwrite(res, OUT_FILE, sep = "\t", na = "")
cat(sprintf("06 OK: %d cells -> %s\n", nrow(res), OUT_FILE))
