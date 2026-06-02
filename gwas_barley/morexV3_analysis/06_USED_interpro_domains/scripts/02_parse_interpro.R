#!/usr/bin/env Rscript
# 02_parse_interpro.R
# Step 02 of 06_USED_interpro_domains.
# Parse the per-sequence InterProScan TSV results into ONE row per gene:
#   gene_id -> matched InterPro + Pfam accessions/descriptions + GO terms.
# Every one of the 45 genes gets a row (NO_MATCH if InterProScan found nothing).
#
# InterProScan TSV columns (no header):
#  1 protein_acc  2 md5  3 length  4 analysis  5 signature_acc  6 signature_desc
#  7 start  8 stop  9 score  10 status  11 date  12 interpro_acc  13 interpro_desc
#  14 GO (with --goterms)  15 pathways (with --pathways)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

base_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/06_USED_interpro_domains"
raw_dir  <- file.path(base_dir, "intermediates/raw")
faa      <- file.path(base_dir, "inputs/fdr_proteins.faa")
out_tsv  <- file.path(base_dir, "results/tables/fdr_gene_interpro.tsv")
dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)

strip_iso <- function(x) sub("\\.[0-9]+$", "", x)   # peptide id -> gene_id

# ---- gene list + protein lengths from the query FASTA (independent of results) ----
fl <- readLines(faa)
hdr <- grep("^>", fl)
ids <- sub("^>", "", sub("\\s.*$", "", fl[hdr]))
starts <- hdr + 1L; ends <- c(hdr[-1] - 1L, length(fl))
lens <- mapply(function(s, e) sum(nchar(fl[s:e])), starts, ends)
genes <- data.frame(seqid = ids, gene_id = strip_iso(ids),
                    protein_length = lens, stringsAsFactors = FALSE)
genes <- genes[order(genes$gene_id), ]

uniq_join <- function(v) paste(unique(v[nzchar(v)]), collapse = "; ")
uniq_semi <- function(v) paste(sort(unique(v[nzchar(v)])), collapse = ";")

# ---- parse each TSV ----
parse_one <- function(seqid) {
  f <- file.path(raw_dir, paste0(seqid, ".tsv.tsv"))
  empty <- list(n_signatures = 0L, member_dbs = "", n_interpro = 0L,
                interpro_ids = "", interpro_descs = "", n_pfam = 0L,
                pfam_ids = "", pfam_descs = "", n_go = 0L, go_terms = "")
  if (!file.exists(f) || file.info(f)$size == 0) return(empty)
  d <- tryCatch(
    read.delim(f, header = FALSE, sep = "\t", quote = "", comment.char = "",
               stringsAsFactors = FALSE, fill = TRUE, colClasses = "character"),
    error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(empty)
  # pad to >=15 columns
  while (ncol(d) < 15) d[[paste0("V", ncol(d) + 1)]] <- ""
  analysis <- d[[4]]; sig_acc <- d[[5]]; sig_desc <- d[[6]]
  ipr_acc  <- d[[12]]; ipr_desc <- d[[13]]; go_raw <- d[[14]]

  ipr_ok <- grepl("^IPR[0-9]+", ipr_acc)
  pfam_ok <- analysis == "Pfam"

  go_ids <- unlist(strsplit(go_raw, "\\|"))
  go_ids <- regmatches(go_ids, regexpr("GO:[0-9]+", go_ids))

  list(
    n_signatures = nrow(d),
    member_dbs   = uniq_semi(analysis),
    n_interpro   = length(unique(ipr_acc[ipr_ok])),
    interpro_ids = uniq_semi(ipr_acc[ipr_ok]),
    interpro_descs = uniq_join(ifelse(ipr_ok, paste0(ipr_acc, ":", ipr_desc), "")),
    n_pfam       = length(unique(sig_acc[pfam_ok])),
    pfam_ids     = uniq_semi(sig_acc[pfam_ok]),
    pfam_descs   = uniq_join(ifelse(pfam_ok, paste0(sig_acc, ":", sig_desc), "")),
    n_go         = length(unique(go_ids)),
    go_terms     = uniq_semi(go_ids)
  )
}

rows <- lapply(genes$seqid, parse_one)

res <- data.frame(
  gene_id        = genes$gene_id,
  protein_length = genes$protein_length,
  n_signatures   = vapply(rows, function(x) x$n_signatures, integer(1)),
  member_dbs     = vapply(rows, function(x) x$member_dbs, character(1)),
  n_interpro     = vapply(rows, function(x) x$n_interpro, integer(1)),
  interpro_ids   = vapply(rows, function(x) x$interpro_ids, character(1)),
  interpro_descs = vapply(rows, function(x) x$interpro_descs, character(1)),
  n_pfam         = vapply(rows, function(x) x$n_pfam, integer(1)),
  pfam_ids       = vapply(rows, function(x) x$pfam_ids, character(1)),
  pfam_descs     = vapply(rows, function(x) x$pfam_descs, character(1)),
  n_go           = vapply(rows, function(x) x$n_go, integer(1)),
  go_terms       = vapply(rows, function(x) x$go_terms, character(1)),
  stringsAsFactors = FALSE
)

res$status <- ifelse(res$n_interpro > 0, "INTERPRO_MATCH",
               ifelse(res$n_signatures > 0, "SIGNATURE_ONLY", "NO_MATCH"))

write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ---- report ----
cat("=== Step 02: parse InterProScan -> per-gene table ===\n")
cat("Genes:", nrow(res), " | result TSVs found:",
    length(list.files(raw_dir, pattern = "\\.tsv\\.tsv$")), "\n\n")
cat("status counts:\n"); print(table(res$status))
cat("\ngenes with >=1 InterPro:", sum(res$n_interpro > 0),
    " | with >=1 Pfam:", sum(res$n_pfam > 0),
    " | with >=1 GO:", sum(res$n_go > 0), "\n")
cat("genes with NO match at all:", sum(res$status == "NO_MATCH"), "\n")
cat("\nWrote:", out_tsv, "\n")
