#!/usr/bin/env Rscript
# 01_make_protein_fasta.R
# Step 01 of 05_USED_gene_annotation.
# Copy the local PGSB Morex V3 r3 HC proteome into inputs/ (provenance), then
# subset it to one representative peptide per FDR gene.
#
# Representative rule: among a gene's isoforms take the LONGEST sequence; ties are
# broken by the lowest isoform index. For single-isoform genes this yields the
# canonical ".1"; for the one multi-isoform gene it yields the longest isoform.
#
# Output: inputs/fdr_proteins.faa  (expected 45 sequences)
# Verify: all 45 genes found, count == 45, HORVU.MOREX.r3.1HG0079280 == 206 aa.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

base_dir   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/05_USED_gene_annotation"
inputs_dir <- file.path(base_dir, "inputs")
src_proteome <- "/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.aa.fa"

genes_txt  <- file.path(inputs_dir, "fdr_genes.txt")
local_copy <- file.path(inputs_dir, "Hv_Morex.pgsb.Jul2020.HC.aa.fa")
out_faa    <- file.path(inputs_dir, "fdr_proteins.faa")

stopifnot(file.exists(src_proteome), file.exists(genes_txt))

# Provenance copy of the proteome
if (!file.exists(local_copy)) file.copy(src_proteome, local_copy, overwrite = FALSE)

genes <- readLines(genes_txt)
genes <- genes[nzchar(genes)]

# ---- Minimal FASTA reader (base R, no Biostrings dependency) ----
read_fasta <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  ids  <- sub("^>", "", lines[hdr_idx])
  ids  <- sub("\\s.*$", "", ids)            # header up to first whitespace
  starts <- hdr_idx + 1L
  ends   <- c(hdr_idx[-1] - 1L, length(lines))
  seqs <- mapply(function(s, e) {
    if (s > e) return("")
    paste0(lines[s:e], collapse = "")
  }, starts, ends)
  data.frame(id = ids, seq = seqs, stringsAsFactors = FALSE)
}

fa <- read_fasta(src_proteome)
# Strip a single trailing stop-codon "*" (PGSB peptides carry it; Ensembl reports
# length without it, and DIAMOND/BLASTP should not see stop characters).
n_trailing_stop <- sum(grepl("\\*$", fa$seq))
fa$seq <- sub("\\*$", "", fa$seq)
n_internal_stop <- sum(grepl("\\*", fa$seq))
if (n_internal_stop > 0) cat("WARNING: ", n_internal_stop,
                             " sequence(s) still contain internal '*'\n", sep = "")

# peptide id = gene_id + "." + isoform  ->  derive gene and isoform index
fa$gene    <- sub("\\.[0-9]+$", "", fa$id)
fa$isoform <- suppressWarnings(as.integer(sub("^.*\\.([0-9]+)$", "\\1", fa$id)))
fa$len     <- nchar(fa$seq)

# Restrict to FDR genes and pick representative per gene
sub_fa <- fa[fa$gene %in% genes, , drop = FALSE]

missing <- setdiff(genes, unique(sub_fa$gene))
if (length(missing)) {
  cat("ERROR: genes not found in proteome:\n"); print(missing)
  stop("Missing genes; aborting.")
}

# Order so the representative is first within each gene: longest desc, then isoform asc
sub_fa <- sub_fa[order(sub_fa$gene, -sub_fa$len, sub_fa$isoform), , drop = FALSE]
rep_fa <- sub_fa[!duplicated(sub_fa$gene), , drop = FALSE]
rep_fa <- rep_fa[order(rep_fa$gene), , drop = FALSE]

# Report multi-isoform genes
n_iso <- tapply(sub_fa$id, sub_fa$gene, length)
multi <- names(n_iso[n_iso > 1])

# ---- Write FASTA (60-col wrap) ----
wrap60 <- function(s) paste(substring(s, seq(1, nchar(s), 60),
                                      seq(60, nchar(s) + 59, 60)), collapse = "\n")
con <- file(out_faa, "w")
for (i in seq_len(nrow(rep_fa))) {
  cat(">", rep_fa$id[i], "\n", wrap60(rep_fa$seq[i]), "\n", sep = "", file = con)
}
close(con)

# ---- Report ----
cat("=== Step 01: build FDR-gene protein FASTA ===\n")
cat("Source proteome:", src_proteome, "\n")
cat("Local copy:     ", local_copy, "\n")
cat("FDR genes requested:", length(genes), "\n")
cat("Representative sequences written:", nrow(rep_fa), "\n")
cat("Multi-isoform genes:", length(multi), "\n")
if (length(multi)) for (g in multi) {
  iso <- sub_fa[sub_fa$gene == g, c("id", "len")]
  chosen <- rep_fa$id[rep_fa$gene == g]
  cat("  ", g, "-> chosen", chosen,
      "(", nrow(iso), "isoforms; lengths",
      paste(iso$len, collapse = "/"), ")\n")
}

chk_gene <- "HORVU.MOREX.r3.1HG0079280"
chk_len  <- rep_fa$len[rep_fa$gene == chk_gene]
cat("\nCheck", chk_gene, "length:", chk_len, "aa (expected 206)\n")
stopifnot(nrow(rep_fa) == length(genes))
stopifnot(length(chk_len) == 1 && chk_len == 206)

cat("\nWrote:", out_faa, "\n")
cat("All checks passed.\n")
