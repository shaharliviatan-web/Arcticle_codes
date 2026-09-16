# =============================================================================
# 02_screen_elite_lines.R -- call rate of every candidate elite line at our sites
# =============================================================================
# Added 2026-09-14. The first selection of five lines (2026-09-10) was made on
# fame alone and two of them turned out to be poorly covered at our sites (RGT
# Planet 53%, Propino 63% called), drawn as grey tiles. The DivBrowse panel is
# unimputed low-coverage WGS, so a line's missingness is a property of its
# sequencing depth and can be measured before choosing.
#
# For every line in the candidate pool (inputs/elite_pool_samples.tsv, from step
# 01) and every target gene:
#
#   shared sites = wild SNPs that also carry an elite record with IDENTICAL REF
#                  and ALT (key CHROM:POS:REF:ALT) -- exactly the columns on which
#                  an elite row can be compared with the wild group rows
#   called       = REF, ALT or HET genotype at a shared site (./. is not called)
#   call rate    = called / shared sites, per gene and pooled over all genes
#
# Selection on call rate is a genotype-QUALITY criterion only: it never looks at
# which allele or haplotype a line carries, so it cannot steer the comparison.
# The choice of WHICH passing lines to show ("famous, currently grown") stays a
# manual, documented decision in config/elite_lines.tsv.
#
# HARD CHECK: every line configured in config/elite_lines.tsv must reach
# CALL_RATE_MIN pooled. The table is written first, then the run stops if any
# configured line fails -- so the table can be used to pick a replacement.
#
# Output -> results/tables/Table_elite_line_screen.tsv   one row per pool line
#
# Author : Shahar Liviatan
# Created: 2026-09-14
# =============================================================================

suppressPackageStartupMessages({ library(vcfR) })

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(script_dir, "_load_params.R"))
P <- load_params(file.path(script_dir, "..", "config", "params.sh"))
Sys.setenv(TMPDIR = P$TMPDIR)

genes <- target_genes_df(P)
elite <- read_elite_lines(P)
pool  <- read.delim(file.path(P$DIR_INPUTS, "elite_pool_samples.tsv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
gw    <- read.delim(P$STEP04_GENE_WINDOWS, stringsAsFactors = FALSE)
thr   <- as.numeric(P$CALL_RATE_MIN)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

site_keys <- function(v) {
  fx <- as.data.frame(vcfR::getFIX(v), stringsAsFactors = FALSE)
  paste(fx$CHROM, fx$POS, fx$REF, fx$ALT, sep = ":")
}

called_mat <- NULL    # rows = pool lines, one column per gene (called counts)
shared_n   <- c()

for (i in seq_len(nrow(genes))) {
  gid <- genes$gene_id[i]; short <- genes$short_name[i]
  trait <- gw$trait[gw$gene_id == gid][1]

  wv <- suppressMessages(vcfR::read.vcfR(
          file.path(P$STEP04_RAW_VCF_DIR, trait, paste0(gid, ".vcf.gz")), verbose = FALSE))
  ev <- suppressMessages(vcfR::read.vcfR(
          file.path(P$DIR_ELITE_TRIMMED, paste0(gid, ".pool.vcf.gz")), verbose = FALSE))

  shared <- intersect(site_keys(wv), site_keys(ev))
  gt <- vcfR::extract.gt(ev, element = "GT", as.numeric = FALSE)
  rownames(gt) <- site_keys(ev)
  gt <- gt[shared, , drop = FALSE]

  called <- apply(gt, 2, function(x) sum(!is.na(x) & !(x %in% c("./.", ".|.", "."))))

  if (is.null(called_mat)) {
    called_mat <- matrix(NA_real_, nrow = length(called), ncol = nrow(genes),
                         dimnames = list(names(called), genes$short_name))
  }
  stopifnot(setequal(rownames(called_mat), names(called)))
  called_mat[names(called), short] <- called
  shared_n[short] <- length(shared)
  msg(short, ": ", length(shared), " shared sites, ", ncol(gt), " pool lines")
}

stopifnot(setequal(rownames(called_mat), pool$sample_id))

out <- data.frame(sample_id = rownames(called_mat), stringsAsFactors = FALSE)
for (s in genes$short_name) {
  out[[paste0("called_", s)]] <- called_mat[, s]
  out[[paste0("call_rate_", s)]] <- round(called_mat[, s] / shared_n[[s]], 3)
}
out$called_all   <- rowSums(called_mat)
out$shared_all   <- sum(shared_n)
out$call_rate_all <- round(out$called_all / out$shared_all, 3)
out$min_gene_call_rate <- apply(out[, paste0("call_rate_", genes$short_name), drop = FALSE], 1, min)
out$passes_threshold <- out$call_rate_all >= thr
out$configured <- out$sample_id %in% elite$sample_id

out <- merge(pool, out, by = "sample_id")
out <- out[order(-out$call_rate_all, -out$min_gene_call_rate, out$accession_name), ]

write.table(out, file.path(P$DIR_TABLES, "Table_elite_line_screen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

msg("pool lines passing call rate >= ", thr, ": ", sum(out$passes_threshold), " of ", nrow(out))
cfg <- out[out$configured, c("accession_name", "sample_id", "call_rate_all",
                             "min_gene_call_rate", "passes_threshold")]
print(cfg, row.names = FALSE)

bad <- cfg[!cfg$passes_threshold, ]
if (nrow(bad))
  stop("configured elite line(s) below CALL_RATE_MIN=", thr, ": ",
       paste(bad$accession_name, collapse = ", "),
       ". Pick replacements from results/tables/Table_elite_line_screen.tsv.")

msg("DONE -> ", file.path(P$DIR_TABLES, "Table_elite_line_screen.tsv"))
