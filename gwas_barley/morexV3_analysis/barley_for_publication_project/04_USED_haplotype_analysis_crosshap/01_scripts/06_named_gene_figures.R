#!/usr/bin/env Rscript
# =============================================================================
# 06_named_gene_figures.R
#
# WHAT   Re-render each significant gene's crosshap figures with its FUNCTIONAL NAME
#        in the title, and merge the violin/tree page and the heatmap page into ONE
#        PDF per gene, filed under its trait.
#
#        From this step a gene is referred to by its name, not its accession: the
#        name appears in the filename and in the header inside the PDF. The serial
#        number (rank by BH q) is kept as a prefix so files still sort by significance.
#
# READS  05_.../08_USED_annotation_master/results/tables/Table_significant_genes_paper.tsv
#        04_runs/<run>/Cache/<trait>/<gene>/MGmin_2/HapObject.rds
# WRITES 04_runs/<run>/Significant_genes/by_trait/<trait>/NN__<Name>__<gene>.pdf
#
# The figures are RE-RENDERED rather than re-titled after the fact, because the title
# is drawn into the page by the plotting code. The renderers themselves are untouched.
# =============================================================================
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table) })

argv <- commandArgs(trailingOnly = FALSE)
fa <- argv[grepl("^--file=", argv)]
sdir <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
source(file.path(sdir,"R","utils.R")); source(file.path(sdir,"R","plot_combined_pdf.R"))
source(file.path(sdir,"R","plot_heatmaps.R")); source(file.path(sdir,"R","run_crosshap.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg <- yaml::read_yaml(if (length(args)>=1) args[1] else file.path(dirname(sdir),"00_config","config.yaml"))
MG <- as.integer(cfg$mgmin_values[[1]]); EPS <- as.numeric(cfg$epsilon_vector[[1]])
label <- paste0("Haplotypes_MGmin",MG,"_E",EPS)

run  <- file.path(cfg$output_root, "04_runs", cfg$run_id)
SG   <- file.path(run, "Significant_genes")
TBL  <- file.path(dirname(cfg$output_root),
                  "05_USED_gene_annotation_analysis/08_USED_annotation_master",
                  "results/tables/Table_significant_genes_paper.tsv")
stopifnot(file.exists(TBL))
d <- fread(TBL)

outroot <- file.path(SG, "by_trait")
unlink(outroot, recursive = TRUE); dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
tmpd <- file.path(cfg$tmpdir %||% "/mnt/data/shahar/.tmp", paste0("namedfig_", Sys.getpid()))
dir.create(tmpd, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmpd, recursive = TRUE), add = TRUE)

# Filesystem-safe short display name from the functional call.
safe_name <- function(x) {
  if (is.na(x) || !nzchar(x)) return("Uncharacterized")
  x <- sub(";.*$", "", x)                    # first clause of a multi-domain list
  x <- sub("^IPR[0-9]+:", "", x)             # drop a leading InterPro accession
  x <- gsub("[^A-Za-z0-9 ._+-]", " ", x)
  x <- gsub("\\s+", "_", trimws(x))
  substr(x, 1, 48)
}

ok <- 0; failed <- character(0)
for (i in seq_len(nrow(d))) {
  tr <- d$trait[i]; gid <- d$gene_id[i]; gshort <- d$gene_short[i]; sn <- d$serial_no[i]
  nm_raw <- d$final_call[i]
  nm <- safe_name(nm_raw)
  disp <- if (is.na(nm_raw) || !nzchar(nm_raw)) "Uncharacterized (conserved, unnamed)" else nm_raw

  rds <- file.path(run, "Cache", tr, gid, paste0("MGmin_", MG), "HapObject.rds")
  if (!file.exists(rds)) { failed <- c(failed, paste(gid, "no cache")); next }
  res <- tryCatch(readRDS(rds), error = function(e) NULL)
  if (is.null(res) || is.null(res$HapObject[[label]])) { failed <- c(failed, paste(gid,"no HapObject")); next }

  # Title carried into BOTH pages: name first, then identity and statistics.
  ttl <- sprintf("#%02d  %s  |  %s  |  %s\nq=%.2g   eta2=%.3f   top-vs-bottom %+.2f SD   lead %s",
                 sn, disp, tr, gshort, d$fdr_q[i], d$eta_squared[i],
                 d$delta_top_bottom_sd[i], d$lead_SNP[i])

  p1 <- file.path(tmpd, sprintf("%02d_a.pdf", sn)); p2 <- file.path(tmpd, sprintf("%02d_b.pdf", sn))
  gene_file <- paste0(gid, ".vcf.gz")
  okc <- tryCatch({ write_combined_pdf(HapObject=res$HapObject, out_pdf=p1, title=ttl, trait=tr,
            gene_file=gene_file, gene_name=gid, MGmin=MG, epsilon_vector=EPS,
            mgmin_test_stats=NULL, gene_summary_row=NULL); TRUE }, error=function(e){message(e$message); FALSE})
  okh <- tryCatch({ write_heatmaps_for_eps(HapObject=res$HapObject, label=label, gene_file=gene_file,
            title=ttl, trait=tr, eps=EPS, MGmin=MG, raw_path=res$raw_path,
            common_ids=res$common_ids, vcf_raw_ids=res$vcf_raw_ids, out_pdf=p2); TRUE },
            error=function(e){message(e$message); FALSE})
  if (!okc) { failed <- c(failed, paste(gid,"combined render")); next }

  tdir <- file.path(outroot, tr); dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(tdir, sprintf("%02d__%s__%s.pdf", sn, nm, gshort))
  parts <- c(p1, if (okh && file.exists(p2)) p2 else NULL)
  if (length(parts) == 2) {
    rc <- system2("pdfunite", c(shQuote(parts), shQuote(out)), stdout=NULL, stderr=NULL)
    if (rc != 0 || !file.exists(out)) file.copy(p1, out, overwrite = TRUE)
  } else file.copy(p1, out, overwrite = TRUE)
  if (file.exists(out)) { ok <- ok + 1
    cat(sprintf("[%2d] %-11s %s\n", sn, tr, basename(out)))
  } else failed <- c(failed, paste(gid,"merge"))
}

cat(sprintf("\nWrote %d merged PDFs (violin/tree + heatmap in one file) to\n  %s\n", ok, outroot))
for (tr in sort(unique(d$trait)))
  cat(sprintf("  %-11s %d\n", tr, length(list.files(file.path(outroot,tr), pattern="\\.pdf$"))))
if (length(failed)) { cat("\nFAILED:\n"); cat(paste0("  ",failed,collapse="\n"),"\n") }
