#!/usr/bin/env Rscript
# ============================================================================
# 04_flank_sensitivity.R
#
# WHAT   Re-run the intersection at several flank sizes and tabulate what each
#        would have returned. Nothing here changes the main result -- it exists
#        so the FLANK_BP choice is defensible with numbers rather than asserted,
#        and so the cost of the choice is visible instead of hidden.
#
# READS  intermediates/loci.tsv, intermediates/genes_1to7H.gff
#        config/params.sh -> FLANK_SENSITIVITY_SET
# WRITES results/tables/flank_sensitivity.tsv        genome-wide, per flank
#        results/tables/flank_sensitivity_by_trait.tsv  per trait x flank
#
# The row matching FLANK_BP is flagged is_chosen=TRUE.
# ============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_load_params.R"))
p <- load_params()

base   <- p$STEP03_BASE
inter  <- file.path(base, "intermediates")
outdir <- file.path(base, "results", "tables")
genes_gff <- file.path(inter, "genes_1to7H.gff")
stopifnot(file.exists(genes_gff))

loci <- read.delim(file.path(inter, "loci.tsv"), stringsAsFactors = FALSE)

cat("=== 04_flank_sensitivity.R ===\n")
cat("Flanks tested (bp): ", paste(p$FLANK_SENSITIVITY_SET, collapse = ", "), "\n\n", sep = "")

tmp <- file.path(p$TMPDIR, sprintf("step03_flank_%d", Sys.getpid()))
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

gw <- list(); bt <- list()

for (f in p$FLANK_SENSITIVITY_SET) {
  bed <- data.frame(chrom = loci$chr,
                    start = pmax(1L, loci$locus_start - f) - 1L,
                    end   = loci$locus_end + f,
                    name  = loci$locus_id, stringsAsFactors = FALSE)
  bed <- bed[order(bed$chrom, bed$start), ]
  bf <- file.path(tmp, sprintf("f%d.bed", f))
  write.table(bed, bf, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # -wa -wb for the gene rows; -c for per-locus counts including zeros.
  hits <- system2(p$BEDTOOLS, c("intersect","-a",shQuote(bf),"-b",shQuote(genes_gff),"-wa","-wb"),
                  stdout = TRUE)
  cnts <- system2(p$BEDTOOLS, c("intersect","-a",shQuote(bf),"-b",shQuote(genes_gff),"-c"),
                  stdout = TRUE)

  n_by_locus <- setNames(
    as.integer(vapply(strsplit(cnts, "\t"), `[`, character(1), 5L)),
    vapply(strsplit(cnts, "\t"), `[`, character(1), 4L))

  if (length(hits)) {
    parts   <- strsplit(hits, "\t")
    lid     <- vapply(parts, `[`, character(1), 4L)
    attr9   <- vapply(parts, `[`, character(1), 13L)
    gid     <- sub("^.*?gene_id=([^;]*).*$", "\\1", attr9)
    gid[!grepl("gene_id=", attr9)] <- NA_character_
    has_des <- grepl("description=", attr9)
  } else { lid <- character(0); gid <- character(0); has_des <- logical(0) }

  trait_of <- setNames(loci$trait, loci$locus_id)
  searched_mb <- sum(bed$end - bed$start) / 1e6

  gw[[as.character(f)]] <- data.frame(
    flank_bp            = f,
    flank_kb            = f/1000,
    total_search_Mb     = round(searched_mb, 2),
    n_loci              = nrow(loci),
    n_loci_with_genes   = sum(n_by_locus > 0),
    n_loci_without_genes= sum(n_by_locus == 0),
    n_gene_locus_rows   = length(lid),
    n_unique_genes      = length(unique(gid)),
    n_with_annotation   = sum(has_des),
    genes_per_Mb        = round(length(lid)/max(searched_mb, 1e-9), 1),
    is_chosen           = (f == p$FLANK_BP),
    stringsAsFactors = FALSE)

  bt[[as.character(f)]] <- do.call(rbind, lapply(sort(unique(loci$trait)), function(tr) {
    lt <- loci$locus_id[loci$trait == tr]; sel <- lid %in% lt
    data.frame(flank_bp = f, flank_kb = f/1000, trait = tr,
               n_loci = length(lt),
               n_loci_with_genes    = sum(n_by_locus[lt] > 0),
               n_loci_without_genes = sum(n_by_locus[lt] == 0),
               n_genes              = sum(sel),
               n_with_annotation    = sum(has_des[sel]),
               is_chosen            = (f == p$FLANK_BP),
               stringsAsFactors = FALSE)
  }))
  cat(sprintf("  flank %6d bp -> %3d genes, %2d/%d loci empty, %5.2f Mb searched%s\n",
              f, length(lid), sum(n_by_locus == 0), nrow(loci), searched_mb,
              if (f == p$FLANK_BP) "   <- CHOSEN" else ""))
}

gw_df <- do.call(rbind, gw); bt_df <- do.call(rbind, bt)
rownames(gw_df) <- NULL; rownames(bt_df) <- NULL
write.table(gw_df, file.path(outdir, "flank_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bt_df, file.path(outdir, "flank_sensitivity_by_trait.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nflank_sensitivity.tsv:\n"); print(gw_df, row.names = FALSE)
cat("\nWrote flank_sensitivity.tsv and flank_sensitivity_by_trait.tsv\n")
