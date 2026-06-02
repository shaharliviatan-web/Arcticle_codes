#!/usr/bin/env Rscript
# 00_lock_fdr_genes.R
# Step 00 of 05_USED_gene_annotation.
# Lock the input set: select significant_fdr == TRUE from the step-04 review CSV.
# Outputs:
#   inputs/fdr_genes.txt        -> unique gene_id, one per line (expected 45)
#   inputs/fdr_genes_table.tsv  -> one row per (trait, gene_id) (expected 48)
# Carries the legacy annotation as legacy_annotation (comparison only, never the call).

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

base_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/05_USED_gene_annotation"
src_csv  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/05_results/tables/gene_annotation_review.csv"

inputs_dir <- file.path(base_dir, "inputs")
dir.create(inputs_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(src_csv))
df <- read.csv(src_csv, stringsAsFactors = FALSE, check.names = FALSE)

# Robust TRUE detection (handles logical or character "TRUE")
is_true <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  toupper(trimws(as.character(x))) == "TRUE"
}

fdr <- df[is_true(df$significant_fdr), , drop = FALSE]

# Normalize Unicode punctuation in the legacy notes to ASCII so all outputs are
# pure-ASCII (the legacy text uses en/em dashes and curly quotes). Comparison-only
# column, so meaning is preserved while the file stays portable.
ascii_punct <- function(x) {
  # code points kept as intToUtf8() so this source file stays pure ASCII
  x <- gsub(intToUtf8(0x2013), "-", x)   # en dash  -> hyphen
  x <- gsub(intToUtf8(0x2014), "-", x)   # em dash  -> hyphen
  x <- gsub(intToUtf8(0x2018), "'", x)   # left single quote  -> '
  x <- gsub(intToUtf8(0x2019), "'", x)   # right single quote -> '
  x <- gsub(intToUtf8(0x201C), "\"", x)  # left double quote  -> "
  x <- gsub(intToUtf8(0x201D), "\"", x)  # right double quote -> "
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")  # drop any rest
  x
}
fdr$annotation <- ascii_punct(fdr$annotation)

# Deliverable: one row per (trait, gene_id), legacy annotation renamed.
out_cols <- data.frame(
  trait                  = fdr$trait,
  gene_id                = fdr$gene_id,
  lead_SNP               = fdr$lead_SNP,
  trait_fdr_p            = fdr$trait_fdr_p,
  significant_bonferroni = fdr$significant_bonferroni,
  legacy_annotation      = fdr$annotation,
  stringsAsFactors       = FALSE
)
out_cols <- out_cols[order(out_cols$trait, out_cols$gene_id), , drop = FALSE]

tsv_path <- file.path(inputs_dir, "fdr_genes_table.tsv")
write.table(out_cols, tsv_path, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

# Unique gene ids
uniq_genes <- sort(unique(out_cols$gene_id))
txt_path <- file.path(inputs_dir, "fdr_genes.txt")
writeLines(uniq_genes, txt_path)

# ---- Report ----
cat("=== Step 00: lock FDR gene set ===\n")
cat("Source CSV:", src_csv, "\n")
cat("Total (trait, gene) FDR rows:", nrow(out_cols), "\n")
cat("Unique genes:", length(uniq_genes), "\n\n")

cat("Per-trait (trait, gene) rows:\n")
print(table(out_cols$trait))
cat("\nPer-trait unique genes:\n")
print(tapply(out_cols$gene_id, out_cols$trait, function(x) length(unique(x))))

# Genes shared across >1 trait
shared <- names(which(tapply(out_cols$trait, out_cols$gene_id,
                             function(x) length(unique(x))) > 1))
cat("\nGenes shared across >1 trait:", length(shared), "\n")
if (length(shared)) for (g in shared) cat("  ", g, "\n")

cat("\nWrote:\n  ", txt_path, "\n  ", tsv_path, "\n")
