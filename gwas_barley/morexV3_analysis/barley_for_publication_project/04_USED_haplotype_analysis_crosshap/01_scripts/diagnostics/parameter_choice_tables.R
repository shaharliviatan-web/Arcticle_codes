#!/usr/bin/env Rscript
# =============================================================================
# parameter_choice_tables.R
#
# WHAT   Regenerates the epsilon-choice table that justifies the crosshap epsilon in
#        the manuscript's Materials & Methods, and the per-gene data behind it. This is
#        Online Resource / supplementary material, so it must be reproducible from the
#        project -- this script is that reproducibility.
#
#        MGmin is NOT scanned here. It is stated as a parameter, not justified by a
#        table (MGmin = 1 yields no haplotypes at all, so 2 is simply the smallest
#        workable value).
#
# GENOTYPE-ONLY BY CONSTRUCTION. crosshap builds haplotype groups from genotypes
#        alone; the phenotype is never used. No p-value is computed anywhere in this
#        script, so the parameter choice cannot be, and was not, tuned toward
#        significance. Do not add a p-value column.
#
# READS  00_config/config.yaml           (paths, minHap, hetmiss_as, keep_outliers)
#        00_config/gene_windows.tsv      the candidate genes of the current run
#        03_per_gene_vcfs/{raw,imputed}_1000bp/  per-gene VCFs
#        03_per_gene_vcfs/raw_1000bp_manifest.tsv  SNP counts per window
#
# WRITES 04_runs/<run_id>/Diagnostics/
#          epsilon_choice_supplementary.tsv    summary, one row per epsilon
#          epsilon_choice_per_gene.tsv         per gene x epsilon
#
# GRID   epsilon 0.2/0.4/0.6/0.8/1.0, with MGmin held at its production value.
#
# RUNTIME ~5 min for 55 genes. Run under screen.
#
# WHY IT IS NOT PART OF run_all: it re-runs crosshap once per gene per epsilon and
#        changes no result.
#        Re-run it only when the candidate-gene set changes, or before submission.
# =============================================================================
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(crosshap); library(dplyr); library(tibble)
})

argv <- commandArgs(trailingOnly = FALSE)
fa <- argv[grepl("^--file=", argv)]
sdir <- if (length(fa)) dirname(dirname(normalizePath(sub("^--file=", "", fa[1])))) else getwd()
source(file.path(sdir, "R", "utils.R")); source(file.path(sdir, "R", "run_crosshap.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg <- yaml::read_yaml(if (length(args) >= 1) args[1] else
        file.path(dirname(sdir), "00_config", "config.yaml"))

EPS_PROD <- as.numeric(cfg$epsilon_vector[[1]])
MG_PROD  <- as.integer(cfg$mgmin_values[[1]])
EPS_GRID <- c(0.2, 0.4, 0.6, 0.8, 1.0)

gw  <- as.data.table(read_gene_windows(file.path(dirname(sdir), "00_config", "gene_windows.tsv")))
man <- fread(file.path(dirname(sdir), "03_per_gene_vcfs", "raw_1000bp_manifest.tsv"))
OUT <- file.path(cfg$output_root, "04_runs", cfg$run_id, "Diagnostics"); ensure_dir(OUT)
tmp <- file.path(Sys.getenv("TMPDIR"), paste0("paramtab_", Sys.getpid())); ensure_dir(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
N_CAND <- nrow(gw)
cat(sprintf("Candidate genes: %d | production eps=%.2f MGmin=%d\n", N_CAND, EPS_PROD, MG_PROD))

## ---- per gene: the expensive setup once, then every configuration ------------
rows <- list()
for (i in seq_len(nrow(gw))) {
  gf <- gw$gene_file[i]; tr <- gw$trait[i]; gn <- make_gene_name(gf)
  nsnp <- man[gene_id == gn]$n_snps; if (!length(nsnp)) nsnp <- NA_integer_
  cat(sprintf("[%2d/%d] %-11s %s\n", i, nrow(gw), tr, gn)); flush.console()

  fail_all <- function(reason) {
    for (E in EPS_GRID) rows[[length(rows)+1]] <<- data.table(gene=gn, trait=tr, n_snps=nsnp,
      grid="epsilon", epsilon=E, MGmin=MG_PROD, testable=FALSE, n_groups=NA_integer_,
      n_assigned=NA_integer_, assign_rate=NA_real_, note=reason)
  }

  raw <- raw_vcf_path(cfg, tr, gf); imp <- imputed_vcf_path(cfg, tr, gf)
  vr <- tryCatch(read_vcf_robust(raw), error = function(e) NULL)
  vi <- tryCatch(read_vcf_robust(imp), error = function(e) NULL)
  if (is.null(vr) || is.null(vi) || nrow(vr) < 2 || nrow(vi) < 2) { fail_all("no_or_too_few_variants"); next }
  ir <- paste0(vr[["#CHROM"]], ":", vr[["POS"]]); ii <- paste0(vi[["#CHROM"]], ":", vi[["POS"]])
  cm <- intersect(ir, ii)
  if (length(cm) < 2) { fail_all("fewer_than_2_common_variants"); next }
  vr <- vr[ir %in% cm, ]; vi <- vi[ii %in% cm, ]
  vr$ID <- make.unique(paste0(vr[["#CHROM"]], ":", vr[["POS"]]))

  pre <- file.path(tmp, paste0(gn, "_"))
  h <- paste0(pre,"h.vcf"); b <- paste0(pre,"b.vcf"); r <- paste0(pre,"r.vcf")
  writeLines(system(paste0("zgrep \"^#\" ", shQuote(imp)), intern = TRUE), h)
  fwrite(vi, b, sep = "\t", col.names = FALSE, quote = FALSE)
  system(paste("cat", shQuote(h), shQuote(b), ">", shQuote(r)))
  system(paste(shQuote(cfg$plink_bin), "--vcf", shQuote(r), "--r2 square --keep-allele-order",
               "--allow-extra-chr --double-id --silent --out", shQuote(paste0(pre,"ld"))))
  ldf <- paste0(pre, "ld.ld")
  if (!file.exists(ldf)) { fail_all("plink_ld_failed"); next }
  LD <- read_LD(ldf, vcf = vr)

  ph <- fread(pheno_path(cfg, tr), header = FALSE); colnames(ph) <- c("FID","IID","trait")
  pheno <- ph %>% mutate(Ind = as.character(IID)) %>%
           select(Ind, Pheno = trait) %>% distinct(Ind, .keep_all = TRUE)

  # One crosshap call per configuration. NOT a vector sweep: crosshap aborts the whole
  # call if any single epsilon in a vector fails, which silently under-counts coverage.
  one <- function(E, M, grid) {
    res <- tryCatch(suppressWarnings(suppressMessages(run_haplotyping(
             vcf = vr, LD = LD, pheno = pheno, epsilon = E, MGmin = M,
             minHap = cfg$minHap, hetmiss_as = cfg$hetmiss_as,
             keep_outliers = cfg$keep_outliers))), error = function(e) NULL)
    k <- NA_integer_; asg <- NA_integer_; ok <- FALSE; note <- "no_marker_groups"
    if (!is.null(res)) {
      H <- res[[paste0("Haplotypes_MGmin", M, "_E", E)]]
      if (!is.null(H) && !is.null(H$Indfile) && nrow(H$Indfile)) {
        hp <- as.character(H$Indfile$hap); g <- hp[hp != "0"]
        k <- length(unique(g)); asg <- length(g)
        ok <- k >= 2; note <- if (ok) "ok" else "fewer_than_2_groups"
      }
    }
    data.table(gene=gn, trait=tr, n_snps=nsnp, grid=grid, epsilon=E, MGmin=M,
               testable=ok, n_groups=k, n_assigned=asg,
               assign_rate=if (is.na(asg)) NA_real_ else asg/nrow(pheno), note=note)
  }
  for (E in EPS_GRID) rows[[length(rows)+1]] <- one(E, MG_PROD, "epsilon")
}
d <- rbindlist(rows)
fwrite(d[grid=="epsilon"], file.path(OUT, "epsilon_choice_per_gene.tsv"), sep="\t", na="NA")

## ---- summaries ----------------------------------------------------------------
summarise <- function(x, by) {
  t <- x[testable == TRUE]
  s <- t[, .(genes_testable = .N,
             total_accessions_assigned = sum(n_assigned),
             median_assignment_pct = round(100*median(assign_rate), 1),
             mean_assignment_pct   = round(100*mean(assign_rate), 1),
             median_n_haplotype_groups = as.numeric(median(n_groups))), by = by]
  all <- data.table(V1 = sort(unique(x[[by]]))); setnames(all, "V1", by)
  s <- merge(all, s, by = by, all.x = TRUE)
  for (cc in setdiff(names(s), by)) s[is.na(get(cc)), (cc) := 0]
  s[, pct_of_candidates := round(100*genes_testable/N_CAND)]
  setcolorder(s, c(by, "genes_testable", "pct_of_candidates", "total_accessions_assigned",
                   "median_assignment_pct", "mean_assignment_pct", "median_n_haplotype_groups"))
  s[order(get(by))]
}
se <- summarise(d[grid=="epsilon"], "epsilon"); se[, setting_used := epsilon == EPS_PROD]

hdr_common <- c(
"# GENOTYPE-ONLY: crosshap builds haplotype groups from genotypes alone -- the phenotype",
"# is never used -- and no p-values enter this table. The parameter choice therefore",
"# cannot be, and was not, tuned toward significance.",
"#",
"# total_accessions_assigned sums over ALL candidate genes (a gene with no grouping",
"# counts 0), so unlike the per-gene percentages it is not conditioned on which genes",
"# happened to cluster at that setting.",
"#",
"# Regenerate: Rscript 01_scripts/diagnostics/parameter_choice_tables.R")

writeLines(c(
sprintf("# Supplementary table: choice of the DBSCAN epsilon for crosshap haplotyping."),
sprintf("# Measured %s on the %d candidate genes of run %s, MGmin = %d fixed.",
        format(Sys.Date(), "%Y-%m-%d"), N_CAND, cfg$run_id, MG_PROD),
"#",
"# THE TRADE-OFF. Assignment rate and gene coverage move monotonically against each",
"# other. A low epsilon assigns a high proportion of accessions but only in the minority",
"# of genes whose variants cluster tightly, and with fewer haplotype groups. A high",
"# epsilon makes more genes testable but leaves more accessions unassigned in each.",
"# Neither extreme is preferable: an unassigned accession contributes nothing to a",
"# gene's test, and an untestable gene contributes nothing at all.",
sprintf("# eps = %.2f was fixed in advance, applied to every gene, and sits at the midpoint of",
        EPS_PROD),
"# that trade-off. It is NOT the maximum of either metric and is not claimed to be.",
"#",
"# Superseded claim, recorded so it is not reintroduced: an earlier table measured on a",
"# retired 64-gene set reported eps 0.6 as having the highest assignment rate. On the",
"# real gene set that is not so -- an artefact of the wrong gene set and of a grid that",
"# omitted eps < 0.6.", "#", hdr_common),
  file.path(OUT, "epsilon_choice_supplementary.tsv"))
suppressWarnings(write.table(se, file.path(OUT, "epsilon_choice_supplementary.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE, append=TRUE))

cat("\n=== epsilon ===\n"); print(se)
cat(sprintf("\nWrote 2 files to %s\n", OUT))
