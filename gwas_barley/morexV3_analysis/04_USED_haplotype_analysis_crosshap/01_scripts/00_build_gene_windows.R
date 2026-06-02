#!/usr/bin/env Rscript
# 00_build_gene_windows.R
# Build the canonical per-gene input for step 04 from the step-03 candidate genes,
# plus the 290-sample keep-list.
#
# Outputs (00_config/):
#   gene_windows.tsv      one row per gene: coords, +/- window_bp window, and the gene's
#                         ORIGIN (lead_SNP, lead_pos, lead_neg_log10p, class, locus_id,
#                         dist_to_lead_bp, description)
#   samples_keep_290.txt  290 single-ID samples = source VCF (300) minus samples_remove (10)
#
# Window: win_start = max(1, gene_start - window_bp); win_end = gene_end + window_bp.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml) })

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)

config_dir   <- file.path(cfg$output_root, "00_config")
cg_tsv       <- cfg$candidate_genes_tsv
lead_loci_tsv<- file.path(dirname(cg_tsv), "lead_loci.tsv")
window_bp    <- as.integer(cfg$window_bp)
bcftools     <- cfg$bcftools_bin

stopifnot(file.exists(cg_tsv), file.exists(lead_loci_tsv))

## ---- Candidate genes (one row per gene x locus) ----
cg <- read.delim(cg_tsv, stringsAsFactors = FALSE)
need <- c("trait","locus_id","lead_SNP","class","gene_id","chr",
          "gene_start","gene_end","strand","dist_to_lead_bp","description")
stopifnot(all(need %in% names(cg)))

## Defensive: if a gene maps to >1 locus within a trait, keep the nearest lead.
cg <- cg[order(cg$trait, cg$gene_id, abs(cg$dist_to_lead_bp)), ]
key <- paste(cg$trait, cg$gene_id, sep = "\t")
cg1 <- cg[!duplicated(key), ]

## ---- Lead locus metadata (lead_pos, lead_neg_log10p) ----
ll <- read.delim(lead_loci_tsv, stringsAsFactors = FALSE)
ll_keep <- ll[, c("trait","locus_id","lead_pos","lead_neg_log10p")]
m <- merge(cg1, ll_keep, by = c("trait","locus_id"), all.x = TRUE, sort = FALSE)

## ---- Window ----
m$win_start <- pmax(1L, as.integer(m$gene_start) - window_bp)
m$win_end   <- as.integer(m$gene_end) + window_bp

out <- m[, c("trait","gene_id","chr","gene_start","gene_end","strand",
             "win_start","win_end","locus_id","class","lead_SNP","lead_pos",
             "lead_neg_log10p","dist_to_lead_bp","description")]
out <- out[order(out$trait, out$chr, out$gene_start), ]
rownames(out) <- NULL

gw_path <- file.path(config_dir, "gene_windows.tsv")
write.table(out, gw_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat(sprintf("gene_windows.tsv: %d genes\n", nrow(out)))
cat("Per-trait gene counts:\n"); print(table(out$trait))
cat("Per-trait class counts:\n"); print(table(out$trait, out$class))

## ---- 290-sample keep-list ----
all_samples <- system2(bcftools, c("query","-l", shQuote(cfg$source_vcf)), stdout = TRUE)
all_samples <- trimws(all_samples[nzchar(all_samples)])

rm_raw <- readLines(cfg$samples_remove)
rm_raw <- rm_raw[nzchar(trimws(rm_raw))]
rm_ids <- trimws(sub("[[:space:]].*$", "", rm_raw))   # first column (IID == FID here)

keep <- setdiff(all_samples, rm_ids)
keep <- keep[order(keep)]

cat(sprintf("Samples: source=%d, remove=%d, keep=%d\n",
            length(all_samples), length(rm_ids), length(keep)))
if (length(keep) != 290L) {
  stop(sprintf("Expected 290 kept samples, got %d. Check source_vcf / samples_remove.", length(keep)))
}
not_found <- setdiff(rm_ids, all_samples)
if (length(not_found) > 0) {
  cat("WARNING: remove-IDs not present in source VCF: ", paste(not_found, collapse=", "), "\n")
}

keep_path <- file.path(config_dir, "samples_keep_290.txt")
writeLines(keep, keep_path)
cat(sprintf("Wrote %s (%d samples)\n", keep_path, length(keep)))
cat("Wrote ", gw_path, "\n", sep = "")
