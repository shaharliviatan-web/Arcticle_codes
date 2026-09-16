#!/usr/bin/env Rscript
# ============================================================================
# 05_check_publication_genes.R
#
# WHAT   A regression check, not an analysis. Step 04 carries a hand-curated
#        list of the genes the paper's haplotype chapter is actually built on
#        (04_.../00_config/final_genes.tsv). Those genes were selected under the
#        OLD search rule. This script asks, every run: does the current search
#        rule still recover them, and if not, by how much does it miss?
#
#        Run it after any change to FLANK_BP or to the step-01 locus set. If a
#        publication gene silently drops out of the candidate list, you find out
#        here rather than three steps downstream.
#
# READS  04_.../00_config/final_genes.tsv   the curated list (may be absent)
#        intermediates/loci.tsv, genes_1to7H.gff
#        results/tables/candidate_genes.tsv
# WRITES results/tables/publication_gene_recovery.tsv
#
# min_flank_to_recover_bp = the smallest FLANK_BP that would pull the gene into
#        its nearest same-trait locus. 0 means it is already recovered.
# ============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_load_params.R"))
p <- load_params()

base   <- p$STEP03_BASE
outdir <- file.path(base, "results", "tables")
final_tsv <- file.path(p$PROJECT_ROOT, "04_USED_haplotype_analysis_crosshap",
                       "00_config", "final_genes.tsv")

cat("=== 05_check_publication_genes.R ===\n")
if (!file.exists(final_tsv)) {
  cat("No curated gene list at", final_tsv, "-- nothing to check. Skipping.\n"); quit(save = "no")
}

fg   <- read.delim(final_tsv, stringsAsFactors = FALSE)
loci <- read.delim(file.path(base, "intermediates", "loci.tsv"), stringsAsFactors = FALSE)
cand <- read.delim(file.path(outdir, "candidate_genes.tsv"), stringsAsFactors = FALSE)

## Gene coordinates straight from the filtered annotation.
gff <- read.delim(file.path(base, "intermediates", "genes_1to7H.gff"),
                  header = FALSE, stringsAsFactors = FALSE)
gff$gene_id <- sub("^.*?gene_id=([^;]*).*$", "\\1", gff$V9)
gi <- match(fg$gene_id, gff$gene_id)

out <- do.call(rbind, lapply(seq_len(nrow(fg)), function(i) {
  if (is.na(gi[i])) return(data.frame(
    trait = fg$trait[i], gene_id = fg$gene_id[i],
    gene_symbol = if ("gene_symbol" %in% names(fg)) fg$gene_symbol[i] else NA,
    chr = NA, gene_start = NA, gene_end = NA, recovered = FALSE,
    matched_locus = NA, gap_to_nearest_locus_bp = NA,
    min_flank_to_recover_bp = NA, note = "gene_id not found in the annotation",
    stringsAsFactors = FALSE))
  chr <- gff$V1[gi[i]]; gs <- gff$V4[gi[i]]; ge <- gff$V5[gi[i]]
  L <- loci[loci$trait == fg$trait[i] & loci$chr == chr, ]
  rec <- any(cand$gene_id == fg$gene_id[i] & cand$trait == fg$trait[i])
  if (!nrow(L)) return(data.frame(
    trait = fg$trait[i], gene_id = fg$gene_id[i],
    gene_symbol = if ("gene_symbol" %in% names(fg)) fg$gene_symbol[i] else NA,
    chr = chr, gene_start = gs, gene_end = ge, recovered = rec,
    matched_locus = NA, gap_to_nearest_locus_bp = NA, min_flank_to_recover_bp = NA,
    note = sprintf("no %s locus on %s at all", fg$trait[i], chr), stringsAsFactors = FALSE))
  gap <- pmax(0L, pmax(L$locus_start - ge, gs - L$locus_end))
  j <- which.min(gap)
  data.frame(
    trait = fg$trait[i], gene_id = fg$gene_id[i],
    gene_symbol = if ("gene_symbol" %in% names(fg)) fg$gene_symbol[i] else NA,
    chr = chr, gene_start = gs, gene_end = ge, recovered = rec,
    matched_locus = L$locus_id[j], gap_to_nearest_locus_bp = as.integer(gap[j]),
    min_flank_to_recover_bp = as.integer(gap[j]),
    note = if (rec) "inside the current search interval"
           else sprintf("outside; needs FLANK_BP >= %s", format(gap[j], big.mark = ",")),
    stringsAsFactors = FALSE)
}))

write.table(out, file.path(outdir, "publication_gene_recovery.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

n_ok <- sum(out$recovered)
cat(sprintf("FLANK_BP = %s\n", format(p$FLANK_BP, big.mark = ",")))
cat(sprintf("Curated publication genes recovered: %d of %d\n\n", n_ok, nrow(out)))
print(out[, c("trait","gene_symbol","gene_id","recovered","matched_locus","min_flank_to_recover_bp")],
      row.names = FALSE)
if (n_ok < nrow(out)) {
  need <- max(out$min_flank_to_recover_bp[!out$recovered], na.rm = TRUE)
  cat(sprintf("\n*** %d curated gene(s) are NOT in the current candidate list.\n", nrow(out) - n_ok))
  cat(sprintf("*** Smallest FLANK_BP that recovers all of them: %s bp\n", format(need, big.mark = ",")))
  cat("*** Set FLANK_BP in config/params.sh and re-run scripts/run_all.sh.\n")
} else cat("\nAll curated publication genes are present in the candidate list.\n")
cat("\nWrote results/tables/publication_gene_recovery.tsv\n")
