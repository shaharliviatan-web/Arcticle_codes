#!/usr/bin/env Rscript
# ============================================================================
# 01_build_intervals.R
#
# WHAT   Turn the step-01 locus handoff table into the gene-search intervals.
#        This script does NOT define loci -- step 01 already did that by LD
#        clumping. All it does is: freeze a copy of the input, assign clean
#        locus IDs, apply FLANK_BP, and write a BED.
#
# READS  config/params.sh -> LOCI_HANDOFF_TSV
#          01_.../results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/
#          Table_loci_for_gene_search.tsv   (34 rows: 32 loci + 2 sub-threshold peaks)
#
# WRITES inputs/loci_handoff_snapshot.tsv   frozen copy of the input, with the
#                                           source path + mtime, so a run stays
#                                           reproducible if step 01 is re-run
#        intermediates/loci.tsv             one row per searched locus
#        intermediates/loci_intervals.bed   0-based half-open, for bedtools
#
# WHY the locus_id rename: step 01 gives the two sub-threshold peaks a search_id
#        equal to their SNP name ("3H:173351655"). A ":" in an ID is unsafe --
#        it reaches filenames and plot titles downstream in step 04 -- so they
#        become protein_X01 / protein_X02 here. The original value is preserved
#        in the source_search_id column, so nothing is lost.
# ============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_load_params.R"))
p <- load_params()

FLANK <- p$FLANK_BP

cat("=== 01_build_intervals.R ===\n")
cat(sprintf("Input : %s\n", p$LOCI_HANDOFF_TSV))
cat(sprintf("FLANK_BP = %d bp each side\n", FLANK))
cat(sprintf("INCLUDE_SUBTHRESHOLD_PEAKS = %s\n\n", p$INCLUDE_SUBTHRESHOLD_PEAKS))

stopifnot(file.exists(p$LOCI_HANDOFF_TSV))
h <- read.delim(p$LOCI_HANDOFF_TSV, stringsAsFactors = FALSE)

need <- c("search_id","trait","chr","lead_SNP","lead_bp","start","end","span_kb",
          "n_member_SNPs","lead_neg_log10_p","lead_MAF","class","include_in_gene_search")
missing <- setdiff(need, names(h))
if (length(missing)) stop("Handoff table is missing columns: ", paste(missing, collapse = ", "))

## ---- Freeze the input -------------------------------------------------------
snap <- file.path(p$STEP03_BASE, "inputs", "loci_handoff_snapshot.tsv")
writeLines(c(
  sprintf("# frozen copy of the step-01 handoff table, taken by 01_build_intervals.R"),
  sprintf("# source : %s", p$LOCI_HANDOFF_TSV),
  sprintf("# mtime  : %s", format(file.mtime(p$LOCI_HANDOFF_TSV), "%Y-%m-%d %H:%M:%S")),
  sprintf("# copied : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
), snap)
suppressWarnings(write.table(h, snap, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
cat(sprintf("Froze input -> inputs/loci_handoff_snapshot.tsv (%d rows)\n", nrow(h)))

## ---- Honour the include flag -----------------------------------------------
h$include_in_gene_search <- as.logical(h$include_in_gene_search)
keep <- h$include_in_gene_search
if (tolower(p$INCLUDE_SUBTHRESHOLD_PEAKS) != "yes") {
  keep <- keep & h$class == "significant_locus"
}
n_drop <- sum(!keep)
h <- h[keep, ]
cat(sprintf("Loci searched: %d (%d excluded)\n", nrow(h), n_drop))
cat("  by class: "); print(table(h$class))

## ---- Clean locus IDs --------------------------------------------------------
# Significant loci already have safe IDs (betaglucan_L01). Sub-threshold peaks
# are numbered per trait as <trait>_X01, X02, ... in position order.
h$source_search_id <- h$search_id
h$locus_id <- h$search_id
sub_i <- h$class == "subthreshold_peak"
if (any(sub_i)) {
  o <- order(h$trait[sub_i], h$chr[sub_i], h$lead_bp[sub_i])
  idx <- which(sub_i)[o]
  h$locus_id[idx] <- ave(h$trait[idx], h$trait[idx],
                         FUN = function(x) sprintf("%s_X%02d", x, seq_along(x)))
  cat("Renamed sub-threshold peak IDs:\n")
  for (i in idx) cat(sprintf("  %-14s -> %s\n", h$source_search_id[i], h$locus_id[i]))
}
if (any(duplicated(h$locus_id))) stop("Duplicate locus_id after renaming.")
if (any(grepl("[^A-Za-z0-9_]", h$locus_id))) stop("Unsafe character in a locus_id.")

## ---- Apply the flank --------------------------------------------------------
h$locus_start <- as.integer(h$start)
h$locus_end   <- as.integer(h$end)
h$search_start <- pmax(1L, h$locus_start - FLANK)
h$search_end   <- h$locus_end + FLANK
h$flank_bp     <- FLANK
h$search_kb    <- round((h$search_end - h$search_start + 1L) / 1000, 3)

loci <- h[, c("locus_id","source_search_id","trait","chr","class",
              "lead_SNP","lead_bp","lead_neg_log10_p","lead_MAF",
              "n_member_SNPs","locus_start","locus_end","span_kb",
              "flank_bp","search_start","search_end","search_kb")]
loci <- loci[order(loci$trait, loci$chr, loci$lead_bp), ]
rownames(loci) <- NULL

write.table(loci, file.path(p$STEP03_BASE, "intermediates", "loci.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- BED (0-based half-open) ------------------------------------------------
bed <- data.frame(chrom = loci$chr,
                  start = loci$search_start - 1L,   # 1-based inclusive -> 0-based
                  end   = loci$search_end,
                  name  = loci$locus_id,
                  stringsAsFactors = FALSE)
bed <- bed[order(bed$chrom, bed$start), ]
write.table(bed, file.path(p$STEP03_BASE, "intermediates", "loci_intervals.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("\nSearch space: %d intervals, %.2f Mb total\n",
            nrow(loci), sum(loci$search_end - loci$search_start + 1) / 1e6))
cat("Per trait:\n"); print(table(loci$trait, loci$class))
cat("\nWrote intermediates/loci.tsv and intermediates/loci_intervals.bed\n")
