#!/usr/bin/env Rscript
# 03_build_annotation_table.R
# Step 03 of 05_USED_gene_annotation.
# Join the per-gene best Swiss-Prot hit onto the (trait, gene) FDR table and
# classify each gene as CONFIDENT_SWISSPROT_HIT or WEAK_OR_NO_HIT.
#
# Best hit per gene = lowest evalue, ties broken by highest bitscore.
# The protein is annotated once per unique gene and replicated across its traits
# (3 genes are shared between fiber and starch).
#
# Outputs:
#   results/tables/fdr_gene_annotation_swissprot.tsv  (one row per trait,gene; 48 rows)
#   results/tables/residual_weak_or_no_hit.txt        (unique WEAK_OR_NO_HIT genes)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

base_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/05_USED_gene_annotation"
tbl_in   <- file.path(base_dir, "inputs/fdr_genes_table.tsv")
blast_in <- file.path(base_dir, "intermediates/fdr_vs_swissprot.diamond.tsv")
out_tsv  <- file.path(base_dir, "results/tables/fdr_gene_annotation_swissprot.tsv")
out_res  <- file.path(base_dir, "results/tables/residual_weak_or_no_hit.txt")

# ---- Thresholds (explicit, adjustable; documented in README) ----
TH_EVALUE   <- 1e-5
TH_QCOVHSP  <- 50
TH_PIDENT   <- 30
UNCHAR_RE   <- "predicted|uncharacterized|hypothetical|DUF[0-9]|unknown function|putative uncharacterized"

stopifnot(file.exists(tbl_in))
gtab <- read.delim(tbl_in, stringsAsFactors = FALSE, check.names = FALSE)

# ---- Read DIAMOND output ----
bcols <- c("qseqid", "sseqid", "stitle", "pident", "qcovhsp",
           "length", "evalue", "bitscore")
if (file.exists(blast_in) && file.info(blast_in)$size > 0) {
  bl <- read.delim(blast_in, header = FALSE, quote = "", comment.char = "",
                   stringsAsFactors = FALSE)
  names(bl) <- bcols
} else {
  bl <- setNames(data.frame(matrix(ncol = length(bcols), nrow = 0)), bcols)
}

bl$gene_id  <- sub("\\.[0-9]+$", "", bl$qseqid)          # strip isoform suffix
bl$pident   <- as.numeric(bl$pident)
bl$qcovhsp  <- as.numeric(bl$qcovhsp)
bl$evalue   <- as.numeric(bl$evalue)
bl$bitscore <- as.numeric(bl$bitscore)

# ---- Best hit per gene: min evalue, then max bitscore ----
bl_ord <- bl[order(bl$gene_id, bl$evalue, -bl$bitscore), , drop = FALSE]
best   <- bl_ord[!duplicated(bl_ord$gene_id), , drop = FALSE]

# ---- Parse Swiss-Prot fields from sseqid / stitle ----
# sseqid: sp|ACCESSION|ENTRYNAME  ;  stitle: "<sseqid> <description> OS=<org> OX=.."
best$sp_accession <- vapply(strsplit(best$sseqid, "\\|"),
                            function(x) if (length(x) >= 2) x[2] else NA_character_,
                            character(1))
nm <- sub("^\\S+\\s+", "", best$stitle)     # drop leading sseqid token
nm <- sub(" OS=.*$", "", nm)                # cut at OS=
best$sp_name <- trimws(nm)
org <- sub("^.*OS=", "", best$stitle)
org <- sub(" OX=.*$", "", org)
best$sp_organism <- ifelse(grepl("OS=", best$stitle), trimws(org), NA_character_)

# characterized_flag: best hit TITLE does NOT match the uncharacterized regex
best$characterized_flag <- !grepl(UNCHAR_RE, best$stitle, ignore.case = TRUE)

# ---- Merge onto the (trait, gene) table (annotate once per gene, replicate) ----
keep <- c("gene_id", "sp_name", "sp_accession", "sp_organism",
          "pident", "qcovhsp", "evalue", "bitscore", "characterized_flag")
m <- merge(gtab, best[, keep], by = "gene_id", all.x = TRUE, sort = FALSE)

# ---- status ----
has_hit <- !is.na(m$evalue)
m$status <- ifelse(
  has_hit &
    m$evalue   <= TH_EVALUE &
    m$qcovhsp  >= TH_QCOVHSP &
    m$pident   >= TH_PIDENT &
    m$characterized_flag %in% TRUE,
  "CONFIDENT_SWISSPROT_HIT", "WEAK_OR_NO_HIT")

# characterized_flag for no-hit genes is undefined -> NA
m$characterized_flag[!has_hit] <- NA

# ---- Final column order ----
final_cols <- c("trait", "gene_id", "lead_SNP", "sp_name", "sp_accession",
                "sp_organism", "pident", "qcovhsp", "evalue", "bitscore",
                "characterized_flag", "status", "legacy_annotation")
m <- m[, final_cols]
m <- m[order(m$trait, m$gene_id), , drop = FALSE]

write.table(m, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# residual = unique genes that are WEAK_OR_NO_HIT (gene-level, not trait-level)
gene_status <- tapply(m$status, m$gene_id,
                      function(s) if (all(s == "WEAK_OR_NO_HIT")) "WEAK_OR_NO_HIT" else "CONFIDENT_SWISSPROT_HIT")
residual <- sort(names(gene_status[gene_status == "WEAK_OR_NO_HIT"]))
writeLines(residual, out_res)

# ---- Report ----
cat("=== Step 03: build annotation table ===\n")
cat("Thresholds: evalue<=", TH_EVALUE, " qcovhsp>=", TH_QCOVHSP,
    " pident>=", TH_PIDENT, " AND characterized\n", sep = "")
cat("Rows (trait, gene):", nrow(m), " | unique genes:",
    length(unique(m$gene_id)), "\n\n")

cat("status counts (per trait,gene row):\n")
print(table(m$trait, m$status))

cat("\nUnique-gene status:\n")
print(table(gene_status))

cat("\nGenes with no Swiss-Prot hit at all:",
    length(setdiff(unique(m$gene_id), best$gene_id)), "\n")
cat("Residual WEAK_OR_NO_HIT unique genes:", length(residual), "\n")

cat("\nWrote:\n  ", out_tsv, "\n  ", out_res, "\n")
