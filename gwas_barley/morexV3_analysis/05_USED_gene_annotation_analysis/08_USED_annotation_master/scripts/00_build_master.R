#!/usr/bin/env Rscript
# 00_build_master.R
# Step 08 of the annotation pipeline: consolidate the three evidence sources
#   05 Swiss-Prot BLASTP, 06 InterProScan, 07 nr BLASTP (residual)
# into ONE master table, one row per (trait, gene). The 3 fiber/starch shared
# genes keep BOTH trait rows (48 rows total). For each row we assign a single
# best functional call + source, an evidence tier, and a provisional
# Swiss-Prot-vs-InterPro agreement flag. Nothing is filtered or dropped here.
#
# Decisions (see README, all adjustable):
#   - near-floor Swiss-Prot identity  = pident < 40  -> caps tier at MEDIUM
#   - agreement compares ANY characterized Swiss-Prot name (even if it failed the
#     step-05 coverage/identity cutoff) against the InterPro/Pfam descriptions
#   - concordance is conservative (exact token match, or substring with len>=4):
#     it may flag a true agreement as discordant rather than over-call concordant
#   - reconciliation priority for final_call:
#       confident Swiss-Prot -> InterPro domain family -> confident nr -> none

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis"
f05  <- file.path(ROOT, "05_USED_gene_annotation/results/tables/fdr_gene_annotation_swissprot.tsv")
f06  <- file.path(ROOT, "06_USED_interpro_domains/results/tables/fdr_gene_interpro.tsv")
f07  <- file.path(ROOT, "07_USED_BLASTP_genes_with_no_annotation_left/results/tables/fdr_residual_nr.tsv")
fcsv <- file.path(ROOT, "04_USED_haplotype_analysis_crosshap/05_results/tables/gene_annotation_review.csv")
out  <- file.path(ROOT, "08_USED_annotation_master/results/tables/fdr_annotation_master.tsv")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

NEAR_FLOOR_PIDENT <- 40   # SP identity below this caps the tier at MEDIUM

rd <- function(p) read.delim(p, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
s5 <- rd(f05); s6 <- rd(f06); s7 <- rd(f07)

# lead_SNP class (significant / marginal) per (trait, gene) from the step-04 CSV
cv <- read.csv(fcsv, stringsAsFactors = FALSE, check.names = FALSE)
cv <- cv[is.logical(cv$significant_fdr) & cv$significant_fdr %in% TRUE, ]
cls <- cv[, c("trait", "gene_id", "class")]
names(cls)[3] <- "lead_SNP_class"

# ---- rename per-source columns to avoid collisions ----
names(s5)[match(c("pident","qcovhsp","evalue","bitscore","characterized_flag","status"), names(s5))] <-
  c("sp_pident","sp_qcovhsp","sp_evalue","sp_bitscore","sp_characterized","sp_status")
names(s6)[match("status", names(s6))] <- "interpro_status"
names(s7)[match(c("pident","qcovs","evalue","bitscore","characterized_flag","status"), names(s7))] <-
  c("nr_pident","nr_qcovs","nr_evalue","nr_bitscore","nr_characterized","nr_status")

# ---- merge onto the 48-row (trait, gene) base ----
m <- merge(s5, cls, by = c("trait","gene_id"), all.x = TRUE)
m <- merge(m, s6[, c("gene_id","protein_length","n_signatures","n_interpro",
                     "interpro_ids","interpro_descs","n_pfam","pfam_ids","pfam_descs",
                     "n_go","go_terms","interpro_status")],
           by = "gene_id", all.x = TRUE)
m <- merge(m, s7[, c("gene_id","nr_title","nr_accession","nr_organism",
                     "nr_pident","nr_qcovs","nr_evalue","nr_status")],
           by = "gene_id", all.x = TRUE)
stopifnot(nrow(m) == nrow(s5))   # no row multiplication

# ---- helper: evidence presence ----
TRUEs <- function(x) !is.na(x) & x %in% TRUE
sp_confident <- m$sp_status == "CONFIDENT_SWISSPROT_HIT"
sp_named     <- TRUEs(m$sp_characterized)                       # a real SP name exists
ipr_present  <- !is.na(m$interpro_status) & m$interpro_status == "INTERPRO_MATCH"
nr_confident <- !is.na(m$nr_status) & m$nr_status == "CONFIDENT_NR_HIT"
near_floor   <- sp_confident & !is.na(m$sp_pident) & m$sp_pident < NEAR_FLOOR_PIDENT

# ---- keyword-overlap concordance (best effort, conservative) ----
STOP <- c("protein","proteins","domain","domains","family","superfamily","subfamily",
          "containing","putative","probable","like","type","related","conserved",
          "uncharacterized","predicted","hypothetical","terminal","group","repeat",
          "repeats","homolog","homologue","isozyme","chain","subunit","motif","region",
          "fold","class","similar","partial","isoform","and","the","with","full","length",
          "binding","activity","unknown","function","product","unnamed")
toks <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  x <- gsub("IPR[0-9]+|PF[0-9]+|G3DSA[:0-9.]+|SSF[0-9]+|PTHR[0-9]+", " ", x) # drop accession tokens
  x <- tolower(x); x <- gsub("[^a-z0-9]+", " ", x)
  w <- unlist(strsplit(x, " +")); w <- w[nchar(w) >= 3]
  w <- w[!grepl("^[0-9]+$", w)]; setdiff(unique(w), STOP)
}
concord <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(FALSE)
  for (x in a) for (y in b) {
    if (x == y) return(TRUE)
    if (nchar(x) >= 4 && grepl(x, y, fixed = TRUE)) return(TRUE)
    if (nchar(y) >= 4 && grepl(y, x, fixed = TRUE)) return(TRUE)
  }
  FALSE
}
agreement <- rep(NA_character_, nrow(m))
for (i in seq_len(nrow(m))) {
  sp_t  <- if (sp_named[i])    toks(m$sp_name[i]) else character(0)
  ipr_t <- if (ipr_present[i]) toks(paste(m$interpro_descs[i], m$pfam_descs[i])) else character(0)
  has_sp <- length(sp_t) > 0; has_ipr <- length(ipr_t) > 0
  agreement[i] <- if (has_sp && has_ipr) {
      if (concord(sp_t, ipr_t)) "concordant" else "discordant"
    } else if (has_sp || has_ipr) "single_source" else "none"
}
m$sp_vs_interpro_agreement <- agreement

# ---- final_call + call_source (reconciliation priority) ----
strip_ipr_ids <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NA_character_)
  parts <- trimws(unlist(strsplit(s, ";")))
  parts <- sub("^(IPR[0-9]+|PF[0-9]+):", "", parts)
  paste(unique(parts[nzchar(parts)]), collapse = "; ")
}
final_call <- character(nrow(m)); call_source <- character(nrow(m))
for (i in seq_len(nrow(m))) {
  if (sp_confident[i]) {
    call_source[i] <- "swissprot"; final_call[i] <- m$sp_name[i]
  } else if (ipr_present[i]) {
    call_source[i] <- "interpro"
    fc <- strip_ipr_ids(m$interpro_descs[i])
    if (is.na(fc) || !nzchar(fc)) fc <- strip_ipr_ids(m$pfam_descs[i])
    final_call[i] <- fc
  } else if (nr_confident[i]) {
    call_source[i] <- "nr"; final_call[i] <- m$nr_title[i]
  } else {
    call_source[i] <- "none"
    any_ev <- (!is.na(m$n_signatures[i]) && m$n_signatures[i] > 0) ||
              (!is.na(m$nr_status[i]) && m$nr_status[i] != "NO_NR_HIT") ||
              sp_named[i]
    final_call[i] <- if (any_ev) "uncharacterized (conserved, unnamed)" else "no homolog found"
  }
}
m$final_call <- final_call; m$call_source <- call_source

# ---- evidence_tier ----
tier <- character(nrow(m))
n_conf <- sp_confident + ipr_present + nr_confident
for (i in seq_len(nrow(m))) {
  if (sp_confident[i] && ipr_present[i] && agreement[i] == "concordant") {
    tier[i] <- if (near_floor[i]) "MEDIUM" else "HIGH"
  } else if (n_conf[i] >= 1) {
    tier[i] <- "MEDIUM"
  } else {
    tier[i] <- "LOW"
  }
}
m$evidence_tier <- tier

# ---- assemble final column order ----
cols <- c("trait","gene_id","lead_SNP","lead_SNP_class",
          "final_call","call_source","evidence_tier","sp_vs_interpro_agreement",
          "sp_name","sp_accession","sp_organism","sp_pident","sp_qcovhsp","sp_evalue","sp_status",
          "interpro_status","n_interpro","interpro_ids","interpro_descs",
          "n_pfam","pfam_ids","pfam_descs","n_go","go_terms",
          "nr_title","nr_accession","nr_pident","nr_qcovs","nr_evalue","nr_status",
          "legacy_annotation")
m <- m[, cols]
m <- m[order(m$trait, m$gene_id), ]
write.table(m, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- report ----
cat("=== Step 08: annotation master table ===\n")
cat("rows:", nrow(m), " (expect 48)   unique genes:", length(unique(m$gene_id)), "\n\n")
cat("call_source:\n"); print(table(m$call_source))
cat("\nevidence_tier:\n"); print(table(m$evidence_tier))
cat("\nsp_vs_interpro_agreement:\n"); print(table(m$sp_vs_interpro_agreement))
cat("\nlead_SNP_class:\n"); print(table(m$lead_SNP_class))
cat("\ntier x call_source:\n"); print(table(m$evidence_tier, m$call_source))
cat("\nWrote:", out, "\n")
