#!/usr/bin/env Rscript
# 01_build_loci.R
# Build per-trait curated SNP tables (significant + marginal) and cluster them
# into loci by single-linkage at 200 kb, per trait per chromosome.
#
# Outputs:
#   inputs/curated_snps.tsv            combined curated input (all traits)
#   intermediates/loci.tsv             one row per locus (full metadata)
#   intermediates/loci_intervals.bed   BED of candidate intervals (span +/- 200 kb)
#
# Window/locus rule (locked, see README):
#   window      = 200 kb each side (genome-wide LD decay to r2=0.2 ~188 kb -> 200 kb)
#   loci        = single-linkage clustering, nearest neighbour <= 200 kb joins
#   interval    = cluster span +/- 200 kb (leftmost member -200kb .. rightmost +200kb)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressWarnings(suppressMessages({}))

WINDOW   <- 200000L  # bp each side
LINK_GAP <- 200000L  # single-linkage join distance (bp)

base_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/03_USED_candidate_genes_around_leading_snps"
lead_tsv <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/01_USED_GWAS_V2_pipeline/results/publication_BonfOnly_BLUP_3PC/tables/lead_snps.tsv"
marg_tsv <- file.path(base_dir, "inputs", "marginal_snps.tsv")

dir.create(file.path(base_dir, "intermediates"), showWarnings = FALSE, recursive = TRUE)

split_snp <- function(snp_id) {
  parts <- strsplit(snp_id, ":", fixed = TRUE)
  chr <- vapply(parts, `[`, character(1), 1L)
  pos <- as.integer(vapply(parts, `[`, character(1), 2L))
  list(chr = chr, pos = pos)
}

## --- Significant SNPs (all above_Bonf=TRUE rows from lead_snps.tsv) ---
lead <- read.delim(lead_tsv, stringsAsFactors = FALSE)
stopifnot(all(lead$above_Bonf == "TRUE" | lead$above_Bonf == TRUE))
sig <- data.frame(
  trait       = lead$trait,
  SNP_id      = lead$SNP_id,
  neg_log10_p = as.numeric(lead$neg_log10_p),
  class       = "significant",
  stringsAsFactors = FALSE
)

## --- Marginal SNPs (curated) ---
marg <- read.delim(marg_tsv, stringsAsFactors = FALSE)
marg$class <- "marginal"
marg <- marg[, c("trait", "SNP_id", "neg_log10_p", "class")]

## --- Combine ---
curated <- rbind(sig, marg)
sp <- split_snp(curated$SNP_id)
curated$chr <- sp$chr
curated$pos <- sp$pos
curated <- curated[order(curated$trait, curated$chr, curated$pos), ]
rownames(curated) <- NULL

write.table(curated[, c("trait", "SNP_id", "chr", "pos", "neg_log10_p", "class")],
            file.path(base_dir, "inputs", "curated_snps.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Curated SNPs: %d total (%d significant, %d marginal)\n",
            nrow(curated), sum(curated$class == "significant"), sum(curated$class == "marginal")))

## --- Single-linkage clustering per trait per chromosome ---
loci_list <- list()
li <- 0L
for (tr in unique(curated$trait)) {
  for (ch in unique(curated$chr[curated$trait == tr])) {
    sub <- curated[curated$trait == tr & curated$chr == ch, ]
    sub <- sub[order(sub$pos), ]
    # assign cluster id: new cluster when gap to previous > LINK_GAP
    cl <- integer(nrow(sub))
    cl[1] <- 1L
    if (nrow(sub) > 1) {
      for (i in 2:nrow(sub)) {
        cl[i] <- if ((sub$pos[i] - sub$pos[i - 1]) <= LINK_GAP) cl[i - 1] else cl[i - 1] + 1L
      }
    }
    for (k in unique(cl)) {
      mem <- sub[cl == k, ]
      li <- li + 1L
      lead_i <- which.max(mem$neg_log10_p)
      span_start <- min(mem$pos)
      span_end   <- max(mem$pos)
      int_start  <- max(1L, span_start - WINDOW)
      int_end    <- span_end + WINDOW
      # locus class: significant if any member is significant, else marginal
      locus_class <- if (any(mem$class == "significant")) "significant" else "marginal"
      loci_list[[li]] <- data.frame(
        trait            = tr,
        locus_id         = sprintf("%s_%s_%d", tr, ch, k),
        chr              = ch,
        lead_SNP         = mem$SNP_id[lead_i],
        lead_pos         = mem$pos[lead_i],
        lead_neg_log10p  = mem$neg_log10_p[lead_i],
        class            = locus_class,
        n_SNPs           = nrow(mem),
        member_SNPs      = paste(mem$SNP_id, collapse = ","),
        span_start       = span_start,
        span_end         = span_end,
        interval_start   = int_start,
        interval_end     = int_end,
        stringsAsFactors = FALSE
      )
    }
  }
}
loci <- do.call(rbind, loci_list)
loci <- loci[order(loci$trait, loci$chr, loci$lead_pos), ]
rownames(loci) <- NULL

write.table(loci, file.path(base_dir, "intermediates", "loci.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- BED of intervals (0-based half-open) for bedtools ---
bed <- data.frame(
  chrom = loci$chr,
  start = loci$interval_start - 1L,  # 1-based inclusive -> 0-based half-open
  end   = loci$interval_end,
  name  = loci$locus_id,
  stringsAsFactors = FALSE
)
write.table(bed, file.path(base_dir, "intermediates", "loci_intervals.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("Loci: %d total (%d significant, %d marginal)\n",
            nrow(loci), sum(loci$class == "significant"), sum(loci$class == "marginal")))
cat("Per trait:\n")
print(table(loci$trait, loci$class))
cat("Wrote inputs/curated_snps.tsv, intermediates/loci.tsv, intermediates/loci_intervals.bed\n")
