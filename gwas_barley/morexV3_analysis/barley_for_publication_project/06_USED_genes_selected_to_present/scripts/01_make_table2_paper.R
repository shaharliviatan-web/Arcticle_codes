#!/usr/bin/env Rscript
# 01_make_table2_paper.R -- Table 2 of the TAG manuscript: the three candidate genes carried forward
# (GPAT6, GH17, PHT4;3), the association signal each lies in, and its haplotype-analysis result.
#
# Created 2026-09-27 (user approval: review item 4 of 10_USED_Paper_writing/new_publishing_paper/build/
# REVIEW_2026-09-27_Results_Discussion.md -- the manuscript table must come from a script, not be typed).
# Read-only on every input; writes only into ../results/tables/.
#
# Inputs
#   04_.../04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv           position, lead SNP, distance, n, groups, q, eta2, delta
#   05_.../08_USED_annotation_master/results/tables/Table_significant_genes_annotated.tsv   sp_pident, sp_qcovhsp
# Outputs
#   ../results/tables/Table_2_genes_carried_forward.tsv   one row per gene, raw values + the formatted cells
#   ../results/tables/Table_2_genes_carried_forward.md    the pipe table exactly as it goes into Results_Discussion.md
#
# Run:  Rscript 06_USED_genes_selected_to_present/scripts/01_make_table2_paper.R   (from the project root, or anywhere)

suppressPackageStartupMessages(library(utils))
P   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
GR  <- file.path(P, "04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv")
AN  <- file.path(P, "05_USED_gene_annotation_analysis/08_USED_annotation_master/results/tables/Table_significant_genes_annotated.tsv")
OUT <- file.path(P, "06_USED_genes_selected_to_present/results/tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# the three genes, in manuscript order, with the short names used in the text and in Fig. 4
genes <- data.frame(
  gene_id    = c("HORVU.MOREX.r3.3HG0301300", "HORVU.MOREX.r3.5HG0487060", "HORVU.MOREX.r3.3HG0301710"),
  short_name = c("GPAT6", "GH17", "PHT4;3"),
  stringsAsFactors = FALSE)

gr <- read.delim(GR, stringsAsFactors = FALSE)
an <- read.delim(AN, stringsAsFactors = FALSE)
g  <- merge(genes, gr, by = "gene_id", sort = FALSE)
g  <- merge(g, an[, c("gene_id", "sp_pident", "sp_qcovhsp")], by = "gene_id", sort = FALSE)
g  <- g[match(genes$gene_id, g$gene_id), ]
stopifnot(nrow(g) == 3, !anyNA(g$fdr_q), !anyNA(g$sp_pident), all(g$significant_fdr))

# ---- formatting helpers ---------------------------------------------------
sup_digits <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴","5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹","-"="⁻")
sci <- function(x) {                      # 5.26e-06 -> "5.3 × 10⁻⁶"
  e <- floor(log10(x)); m <- round(x / 10^e, 1)
  if (m >= 10) { m <- m / 10; e <- e + 1 }
  paste0(formatC(m, format = "f", digits = 1), " × 10",
         paste(sup_digits[strsplit(as.character(e), "")[[1]]], collapse = ""))
}
dist_fmt <- function(bp) ifelse(bp >= 1000, paste0(round(bp / 1000), " kb"), paste0(bp, " bp"))
mb3 <- function(bp) formatC(bp / 1e6, format = "f", digits = 3)
cap <- function(s) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))

g$trait_cell    <- cap(g$trait)
g$position_cell <- paste0(mb3(g$gene_start), "–", mb3(g$gene_end))
g$lead_cell     <- mb3(g$lead_pos)
g$dist_cell     <- dist_fmt(g$dist_to_lead_bp)
g$grouped_cell  <- paste0(g$n_ind, " (", g$n_groups, ")")
g$q_cell        <- vapply(g$fdr_q, sci, "")
g$eta_cell      <- formatC(g$eta_squared, format = "f", digits = 3)
g$delta_cell    <- formatC(g$delta_top_bottom_sd, format = "f", digits = 2)
g$ident_cell    <- paste0(formatC(g$sp_pident, format = "f", digits = 1), " / ", formatC(g$sp_qcovhsp, format = "f", digits = 1))
g$id_cell       <- sub("^HORVU\\.MOREX\\.r3\\.", "", g$gene_id)

# ---- TSV: raw values next to the formatted cells, for checking --------------
tsv <- g[, c("short_name", "gene_id", "trait", "chr", "gene_start", "gene_end", "lead_SNP", "lead_pos", "dist_to_lead_bp",
             "n_ind", "n_groups", "fdr_q", "eta_squared", "delta_top_bottom_sd", "sp_pident", "sp_qcovhsp",
             "position_cell", "lead_cell", "dist_cell", "grouped_cell", "q_cell", "eta_cell", "delta_cell", "ident_cell")]
write.table(tsv, file.path(OUT, "Table_2_genes_carried_forward.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- markdown: exactly the block inserted into the manuscript ----------------
md <- c(
  "**Table 2** The three candidate genes carried forward from the haplotype analysis: the association signal they lie in and the result of their haplotype analysis",
  "",
  "| Gene | Gene ID^a^ | Trait | Chr | Gene position (Mb) | Lead SNP (Mb) | Distance to lead SNP | Accessions grouped (groups) | *q* | η² | Difference (SD)^b^ | Identity / coverage (%)^c^ |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| *%s* | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
          g$short_name, g$id_cell, g$trait_cell, g$chr, g$position_cell, g$lead_cell, g$dist_cell,
          g$grouped_cell, g$q_cell, g$eta_cell, g$delta_cell, g$ident_cell),
  "",
  paste0("^a^ MorexV3 gene IDs, prefix HORVU.MOREX.r3. ^b^ Difference between the highest and the lowest haplotype group, ",
         "in phenotypic standard deviations. ^c^ Best UniProt Swiss-Prot match"))
writeLines(md, file.path(OUT, "Table_2_genes_carried_forward.md"), useBytes = TRUE)

cat("wrote", file.path(OUT, "Table_2_genes_carried_forward.{tsv,md}"), "\n")
cat(md[5:7], sep = "\n")
