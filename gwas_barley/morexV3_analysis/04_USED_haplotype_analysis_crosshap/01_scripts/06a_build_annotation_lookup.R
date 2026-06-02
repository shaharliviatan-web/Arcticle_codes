#!/usr/bin/env Rscript
# 06a_build_annotation_lookup.R
# Build the gene_id -> annotation lookup from the previous tool's hand-curated
# gene_explantions.csv. Each gene data row (17 CSV fields) is followed by free-text
# annotation lines; "Relevance to ..." lines are stripped.
#
# Some genes appear in more than one trait in the source (e.g. the 7HG0729xxx cluster
# in fiber + starch). For those, keep the MOST INFORMATIVE annotation (longest non-empty)
# and flag the gene needs_verification when the per-trait notes are not identical.
#
# Outputs (05_results/tables/):
#   prev_annotations_lookup.tsv     gene_id <tab> annotation   (one row per unique gene_id)
#   prev_annotations_dualtrait.tsv  genes seen in >1 trait: each trait note + chosen + identical?

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

args <- commandArgs(trailingOnly = TRUE)
src_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_Haplotype_Analysis/Haplotype_Analysis_All_Genes/Runs/All_Genes_All_Traits_V2/Stats/gene_explantions.csv"
out_dir <- if (length(args) >= 2) args[2] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/05_results/tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

src <- sub("\r$", "", readLines(src_path, warn = FALSE))
is_gr <- function(l) {
  f <- strsplit(l, ",", fixed = TRUE)[[1]]
  length(f) >= 17 && f[1] %in% c("betaglucan","fiber","protein","starch") && grepl("\\.vcf\\.gz$", f[2])
}
gl <- which(vapply(src, is_gr, logical(1)))
trait <- vapply(src[gl], function(l) strsplit(l, ",", fixed = TRUE)[[1]][1], character(1))
gene  <- vapply(src[gl], function(l) strsplit(l, ",", fixed = TRUE)[[1]][3], character(1))
ends  <- c(gl[-1] - 1, length(src))
ann <- vapply(seq_along(gl), function(i) {
  if (ends[i] < gl[i] + 1) return("")                       # adjacent gene rows -> no annotation
  blk <- trimws(src[(gl[i] + 1):ends[i]])
  blk <- blk[nzchar(blk) & !grepl("^Relevance to", blk) & !grepl("^trait,", blk)]
  paste(blk, collapse = " | ")
}, character(1))

df <- data.frame(gene_id = gene, trait = trait, annotation = ann, stringsAsFactors = FALSE)

# ---- Collapse to one row per gene_id: most informative (longest non-empty) ----
pick_best <- function(a) { a <- a[order(-nchar(a))]; a[1] }   # longest first; "" sinks to bottom
genes <- unique(df$gene_id)
best <- vapply(genes, function(g) pick_best(df$annotation[df$gene_id == g]), character(1))
lookup <- data.frame(gene_id = genes, annotation = best, stringsAsFactors = FALSE)
lookup <- lookup[order(lookup$gene_id), ]
write.table(lookup, file.path(out_dir, "prev_annotations_lookup.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- Dual-trait report (genes seen in >1 trait) ----
dups <- names(which(table(df$gene_id) > 1))
rep_rows <- list(); ver_ids <- character(0)
for (g in dups) {
  sub <- df[df$gene_id == g, ]
  identical_all <- length(unique(sub$annotation)) == 1
  chosen <- pick_best(sub$annotation)
  if (!identical_all) ver_ids <- c(ver_ids, g)
  for (j in seq_len(nrow(sub))) {
    rep_rows[[length(rep_rows) + 1]] <- data.frame(
      gene_id = g, trait = sub$trait[j], annotation = sub$annotation[j],
      all_traits_identical = identical_all, chosen_best = chosen,
      stringsAsFactors = FALSE)
  }
}
dt <- if (length(rep_rows)) do.call(rbind, rep_rows) else
  data.frame(gene_id=character(), trait=character(), annotation=character(),
             all_traits_identical=logical(), chosen_best=character())
write.table(dt, file.path(out_dir, "prev_annotations_dualtrait.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(sort(unique(ver_ids)), file.path(out_dir, "prev_annotations_needs_verification.txt"))

cat(sprintf("Lookup: %d unique genes (from %d source rows).\n", nrow(lookup), nrow(df)))
cat(sprintf("Genes in >1 trait: %d | with NON-identical notes (flagged needs_verification): %d\n",
            length(dups), length(unique(ver_ids))))
if (length(ver_ids)) { cat("Flagged genes:\n"); print(unique(ver_ids)) }
