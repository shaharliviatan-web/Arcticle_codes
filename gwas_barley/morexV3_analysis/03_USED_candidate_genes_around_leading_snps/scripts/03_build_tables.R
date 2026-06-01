#!/usr/bin/env Rscript
# 03_build_tables.R
# Parse the bedtools intersect output, attach locus/lead metadata, compute the
# signed gene->lead distance, and assemble the output tables (no plots).
#
# Inputs:
#   intermediates/loci.tsv
#   intermediates/intersect_raw.tsv   (bedtools -wa -wb: 4 BED cols + 9 GFF cols)
# Outputs (results/tables/):
#   candidate_genes.tsv  lead_loci.tsv  per_trait_summary.tsv  results_chapter_numbers.txt
#
# dist_to_lead_bp: signed nearest-edge distance from the lead SNP to the gene.
#   0  = lead SNP falls inside the gene body
#   <0 = gene lies upstream of the lead (lower coordinate)
#   >0 = gene lies downstream of the lead (higher coordinate)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

base_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/03_USED_candidate_genes_around_leading_snps"
loci   <- read.delim(file.path(base_dir, "intermediates", "loci.tsv"), stringsAsFactors = FALSE)
raw    <- read.delim(file.path(base_dir, "intermediates", "intersect_raw.tsv"),
                     header = FALSE, stringsAsFactors = FALSE)
out_dir <- file.path(base_dir, "results", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(ncol(raw) == 13)
colnames(raw) <- c("int_chr", "int_start0", "int_end", "locus_id",
                   "g_chr", "g_source", "g_feature", "g_start", "g_end",
                   "g_score", "g_strand", "g_frame", "g_attr")

## --- Parse GFF attributes ---
get_attr <- function(attr, key) {
  # attributes are ';'-separated key=value; encoded ';' inside values is '%3B'.
  # regexpr returns -1 for non-matches; build a full-length vector (NA when absent).
  out <- rep(NA_character_, length(attr))
  m   <- regexpr(paste0("(?:^|;)", key, "=([^;]*)"), attr, perl = TRUE)
  hit <- m != -1L
  vals <- regmatches(attr, m)                      # length == sum(hit)
  vals <- sub(paste0("^.*?", key, "="), "", vals)  # strip leading "...key="
  out[hit] <- vals
  out
}
decode <- function(x) {
  vapply(x, function(v) if (is.na(v)) NA_character_ else utils::URLdecode(v),
         character(1), USE.NAMES = FALSE)
}

raw$gene_id     <- get_attr(raw$g_attr, "gene_id")
raw$biotype     <- get_attr(raw$g_attr, "biotype")
raw$description <- decode(get_attr(raw$g_attr, "description"))

## --- Attach locus / lead metadata ---
idx <- match(raw$locus_id, loci$locus_id)
raw$trait      <- loci$trait[idx]
raw$lead_SNP   <- loci$lead_SNP[idx]
raw$lead_pos   <- loci$lead_pos[idx]
raw$class      <- loci$class[idx]
raw$span_start <- loci$span_start[idx]
raw$span_end   <- loci$span_end[idx]

## --- Signed distance gene -> lead and in-span flag ---
lead_pos <- raw$lead_pos; gs <- raw$g_start; ge <- raw$g_end
raw$dist_to_lead_bp <- ifelse(lead_pos >= gs & lead_pos <= ge, 0L,
                       ifelse(ge < lead_pos, ge - lead_pos, gs - lead_pos))
raw$in_cluster_span <- (ge >= raw$span_start) & (gs <= raw$span_end)

## --- Table 1: candidate_genes.tsv (one row per gene x locus) ---
cand <- data.frame(
  trait           = raw$trait,
  locus_id        = raw$locus_id,
  lead_SNP        = raw$lead_SNP,
  class           = raw$class,
  gene_id         = raw$gene_id,
  biotype         = raw$biotype,
  chr             = raw$g_chr,
  gene_start      = raw$g_start,
  gene_end        = raw$g_end,
  strand          = raw$g_strand,
  dist_to_lead_bp = raw$dist_to_lead_bp,
  in_cluster_span = raw$in_cluster_span,
  description     = raw$description,
  stringsAsFactors = FALSE
)
cand <- cand[order(cand$trait, cand$locus_id, cand$gene_start), ]
rownames(cand) <- NULL
write.table(cand, file.path(out_dir, "candidate_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## --- Per-locus gene counts ---
agg <- function(v, f) tapply(v, cand$locus_id, f)
n_genes <- table(cand$locus_id)
n_pc    <- tapply(cand$biotype == "protein_coding", cand$locus_id, sum)
n_desc  <- tapply(!is.na(cand$description), cand$locus_id, sum)

loci$n_genes          <- as.integer(n_genes[loci$locus_id]);          loci$n_genes[is.na(loci$n_genes)] <- 0L
loci$n_protein_coding <- as.integer(n_pc[loci$locus_id]);             loci$n_protein_coding[is.na(loci$n_protein_coding)] <- 0L
loci$n_with_description<- as.integer(n_desc[loci$locus_id]);          loci$n_with_description[is.na(loci$n_with_description)] <- 0L

## --- Table 2: lead_loci.tsv ---
write.table(loci[, c("trait","locus_id","chr","lead_SNP","lead_pos","lead_neg_log10p",
                     "class","n_SNPs","member_SNPs","span_start","span_end",
                     "interval_start","interval_end","n_genes","n_protein_coding","n_with_description")],
            file.path(out_dir, "lead_loci.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- Table 3: per_trait_summary.tsv ---
traits <- sort(unique(loci$trait))
summ <- do.call(rbind, lapply(traits, function(tr) {
  lt <- loci[loci$trait == tr, ]
  ct <- cand[cand$trait == tr, ]
  data.frame(
    trait              = tr,
    n_loci             = nrow(lt),
    n_loci_significant = sum(lt$class == "significant"),
    n_loci_marginal    = sum(lt$class == "marginal"),
    n_genes            = nrow(ct),
    n_protein_coding   = sum(ct$biotype == "protein_coding"),
    n_with_description = sum(!is.na(ct$description)),
    stringsAsFactors = FALSE
  )
}))
write.table(summ, file.path(out_dir, "per_trait_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- Table 4: results_chapter_numbers.txt (manuscript sentences) ---
con <- file(file.path(out_dir, "results_chapter_numbers.txt"), "w")
wl <- function(...) writeLines(paste0(...), con)
wl("Candidate genes around lead SNPs (Morex V3, gene window = lead/cluster span +/- 200 kb).")
wl("Loci defined by single-linkage clustering at 200 kb per trait per chromosome.")
wl("")
wl(sprintf("Across all traits: %d loci (%d significant, %d marginal) span %d candidate genes (%d protein-coding, %d with a functional description).",
           nrow(loci), sum(loci$class=="significant"), sum(loci$class=="marginal"),
           nrow(cand), sum(cand$biotype=="protein_coding"), sum(!is.na(cand$description))))
wl("")
for (tr in traits) {
  lt <- loci[loci$trait == tr, ]
  ct <- cand[cand$trait == tr, ]
  wl(sprintf("%s: %d loci (%d significant, %d marginal), %d candidate genes (%d with description).",
             tr, nrow(lt), sum(lt$class=="significant"), sum(lt$class=="marginal"),
             nrow(ct), sum(!is.na(ct$description))))
  for (i in seq_len(nrow(lt))) {
    g <- ct[ct$locus_id == lt$locus_id[i], ]
    wl(sprintf("  %s (%s lead %s, -log10p=%.3f): %d SNP(s), %d genes.",
               lt$locus_id[i], lt$class[i], lt$lead_SNP[i], lt$lead_neg_log10p[i],
               lt$n_SNPs[i], nrow(g)))
  }
}
close(con)

cat(sprintf("candidate_genes.tsv: %d rows\n", nrow(cand)))
cat(sprintf("lead_loci.tsv: %d loci\n", nrow(loci)))
cat("per_trait_summary.tsv:\n"); print(summ, row.names = FALSE)
cat("Wrote 4 tables to results/tables/\n")
