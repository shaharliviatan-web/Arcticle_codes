#!/usr/bin/env Rscript
# 08b_top_snps_per_chr.R
# For each of the 24 EMMAX runs, write the top N most significant SNPs PER
# CHROMOSOME to its own TSV. With 7 barley chromosomes (1H..7H) and TOP_PER_CHR=15,
# each file has 15 * 7 = 105 rows.
#
# Outputs (24 files, one per run):
#   results/tables/top_snps_per_chr/top15perChr__{trait}__{pheno_type}__pc{N}.tsv
#
# Columns: chrom, rank_in_chrom, SNP, pos, neg_log10_p, p_value, beta, SE
# Rows are sorted by chromosome (1H..7H), then by ascending p-value within chromosome.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
})

PIPE       <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/01_USED_GWAS_V2_pipeline"
PS_DIR     <- file.path(PIPE, "results", "emmax_ps")
TABLES_DIR <- file.path(PIPE, "results", "tables")
TOP_DIR    <- file.path(TABLES_DIR, "top_snps_per_chr")
dir.create(TOP_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_PER_CHR <- 15L

ps_files <- list.files(PS_DIR, pattern = "^morexV3__.+__.+__pc\\d+\\.ps$", full.names = TRUE)
stopifnot(length(ps_files) == 24)
parts <- as.data.table(do.call(rbind, strsplit(basename(ps_files), "__|\\.", perl = FALSE)))
setnames(parts, c("prefix", "trait", "pheno_type", "pcs_str", "ext"))
parts[, pcs := as.integer(sub("^pc", "", pcs_str))]
parts[, path := ps_files]

extract_top_per_chr <- function(i) {
  row <- parts[i, ]
  # .ps columns: SNP, beta, SE, P  (tab-separated)
  d <- fread(row$path, header = FALSE, col.names = c("SNP", "beta", "SE", "P"),
             showProgress = FALSE)
  d <- d[!is.na(P) & P > 0 & P <= 1]
  # Parse "1H:48374" -> chrom="1H", chrom_num=1, pos=48374 (chrom_num for ordering)
  parts_id <- tstrsplit(d$SNP, ":", fixed = TRUE)
  d[, chrom := parts_id[[1]]]
  d[, pos   := as.integer(parts_id[[2]])]
  d[, chrom_num := as.integer(sub("H$", "", chrom))]
  # Top N per chromosome by ascending p-value
  setorder(d, chrom_num, P)
  top <- d[, head(.SD, TOP_PER_CHR), by = chrom_num]
  top[, rank_in_chrom := seq_len(.N), by = chrom_num]
  top[, neg_log10_p := round(-log10(P), 3)]
  top[, p_value := signif(P, 4)]
  top[, beta := signif(as.numeric(beta), 5)]
  top[, SE   := signif(as.numeric(SE),   5)]
  out_cols <- c("chrom", "rank_in_chrom", "SNP", "pos",
                "neg_log10_p", "p_value", "beta", "SE")
  top <- top[, ..out_cols]
  # Per-cell file
  fname <- file.path(TOP_DIR,
                     sprintf("top%dperChr__%s__%s__pc%d.tsv",
                             TOP_PER_CHR, row$trait, row$pheno_type, row$pcs))
  fwrite(top, fname, sep = "\t")
  fname
}

cat(sprintf("[08b] Extracting top-%d SNPs per chromosome for %d runs...\n",
            TOP_PER_CHR, nrow(parts)))
written <- vapply(seq_len(nrow(parts)), extract_top_per_chr, character(1))
cat(sprintf("[08b] Per-cell file count: %d (expected 24)\n",
            length(list.files(TOP_DIR,
                              pattern = sprintf("^top%dperChr__.+\\.tsv$", TOP_PER_CHR)))))
cat(sprintf("[08b] Per-cell files in: %s/  (each has %d x 7 = %d rows)\n",
            TOP_DIR, TOP_PER_CHR, TOP_PER_CHR * 7L))
cat("\n[08b] OK\n")
