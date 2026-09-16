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
#   - evidence_tier was REMOVED on 2026-09-10. It mixed confidence in the protein
#     identification with plausibility as a trait candidate, which are independent.
#     Replaced by plain provenance: annotation_source + annotation_from_check (1/2/3),
#     plus every source's own call in call_check1/2/3, so two differing annotations are
#     both visible and can be judged together.
#   - there is no identity floor and no tier. Swiss-Prot identity/coverage/e-value are
#     REPORTED (sp_pident, sp_qcovhsp, sp_evalue) so quality can be judged directly.
#     (NEAR_FLOOR_PIDENT was removed 2026-09-10: it existed only to cap the retired
#     tier at MEDIUM and had become dead code.)
#   - agreement compares ANY characterized Swiss-Prot name (even if it failed the
#     step-05 coverage/identity cutoff) against the InterPro/Pfam descriptions
#   - concordance is conservative (exact token match, or substring with len>=4):
#     it may flag a true agreement as discordant rather than over-call concordant
#   - reconciliation priority for final_call:
#       confident Swiss-Prot -> InterPro domain family -> confident nr -> none

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
ANNO <- file.path(ROOT, "05_USED_gene_annotation_analysis")
f05  <- file.path(ANNO, "05_USED_gene_annotation/results/tables/fdr_gene_annotation_swissprot.tsv")
f06  <- file.path(ANNO, "06_USED_interpro_domains/results/tables/fdr_gene_interpro.tsv")
f07  <- file.path(ANNO, "07_USED_BLASTP_genes_with_no_annotation_left/results/tables/fdr_residual_nr.tsv")
# 2026-09-09: v1 read the retired step-04 review CSV for the lead-SNP class. That file
# is gone. The locked gene table written by 05_.../scripts/00_lock_fdr_genes.R now supplies
# both the locus class AND the step-04 statistics (q, effect sizes), which belong beside
# an annotation - a functional call is only interpretable next to the effect it explains.
fgenes <- file.path(ANNO, "05_USED_gene_annotation/inputs/fdr_genes_table.tsv")
out  <- file.path(ANNO, "08_USED_annotation_master/results/tables/fdr_annotation_master.tsv")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)


rd <- function(p) read.delim(p, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
s5 <- rd(f05); s6 <- rd(f06)
# Step 07 (nr rescue) can legitimately have nothing to do when Swiss-Prot + InterPro
# already cover every gene. Tolerate an absent/empty table instead of failing.
s7 <- if (file.exists(f07) && file.info(f07)$size > 0) rd(f07) else
      data.frame(gene_id=character(), nr_title=character(), nr_accession=character(),
                 nr_organism=character(), pident=numeric(), qcovs=numeric(),
                 evalue=numeric(), bitscore=numeric(), characterized_flag=logical(),
                 status=character(), stringsAsFactors=FALSE)

# Locus class + step-04 statistics.
# 2026-09-09: these already travel INSIDE the step-05 Swiss-Prot table (00_lock_fdr_genes.R
# puts them there and 03_build_annotation_table.R carries them through), so merging the
# locked gene table in again would collide and produce .x/.y suffixes. Assert their
# presence and just rename locus_class -> lead_SNP_class instead.
need <- c("trait","gene_id","locus_id","lead_SNP","locus_class","dist_to_lead_bp",
          "n_haplotype_groups","kw_p_raw","fdr_q","eta_squared","delta_top_bottom_sd")
miss <- setdiff(need, names(s5))
if (length(miss)) stop("Step-05 table is missing: ", paste(miss, collapse=", "),
                       "\nHas 03_build_annotation_table.R changed?")
names(s5)[names(s5) == "locus_class"] <- "lead_SNP_class"
n_expected <- nrow(s5)

# ---- rename per-source columns to avoid collisions ----
names(s5)[match(c("pident","qcovhsp","evalue","bitscore","characterized_flag","status"), names(s5))] <-
  c("sp_pident","sp_qcovhsp","sp_evalue","sp_bitscore","sp_characterized","sp_status")
names(s6)[match("status", names(s6))] <- "interpro_status"
names(s7)[match(c("pident","qcovs","evalue","bitscore","characterized_flag","status"), names(s7))] <-
  c("nr_pident","nr_qcovs","nr_evalue","nr_bitscore","nr_characterized","nr_status")

# ---- merge onto the 48-row (trait, gene) base ----
m <- s5
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
# ---------------------------------------------------------------------------
# Annotation provenance (replaces the old evidence_tier, removed 2026-09-10).
#
# The tier conflated two different things: how confident we are in the PROTEIN
# IDENTIFICATION, and how plausible the gene is as a candidate FOR THE TRAIT. Those
# are independent judgements, and collapsing them into HIGH/MEDIUM/LOW hid the
# evidence. What matters is simply: which source produced the call, and what did the
# other sources say?
#
# The three sources are checked in a fixed order, and that order is the "check number":
#   check 1 = UniProt Swiss-Prot BLASTP   (curated protein names)
#   check 2 = InterProScan                (domain families + GO)
#   check 3 = NCBI nr BLASTP              (last-resort rescue, residual genes only)
#
# EVERY source's call is now reported in its own column, so when two sources disagree
# both texts are visible side by side and can be judged together -- often they are the
# same protein under different vocabulary.
# ---------------------------------------------------------------------------
call_sp  <- ifelse(sp_confident,  m$sp_name,       NA_character_)
call_ipr <- ifelse(ipr_present,   m$interpro_descs, NA_character_)
call_nr  <- ifelse(nr_confident,  m$nr_title,      NA_character_)

m$call_check1_swissprot <- call_sp
m$call_check2_interpro  <- call_ipr
m$call_check3_nr        <- call_nr

m$annotation_from_check <- ifelse(!is.na(call_sp), 1L,
                           ifelse(!is.na(call_ipr), 2L,
                           ifelse(!is.na(call_nr),  3L, NA_integer_)))
m$annotation_source <- ifelse(!is.na(call_sp), "check1_swissprot",
                       ifelse(!is.na(call_ipr), "check2_interpro",
                       ifelse(!is.na(call_nr),  "check3_nr", "none")))
m$n_sources_with_call <- (!is.na(call_sp)) + (!is.na(call_ipr)) + (!is.na(call_nr))

# TRUE when more than one source produced a call and the texts are not obviously the
# same. Flags a pair for JOINT REVIEW; it does not decide anything and never changes
# the call. The comparison is textual, so it errs toward flagging.
m$needs_review_two_calls <- m$n_sources_with_call >= 2 &
  !is.na(m$sp_vs_interpro_agreement) & m$sp_vs_interpro_agreement == "discordant"

# ---- assemble final column order ----
# 2026-09-09: "legacy_annotation" dropped (v1's non-reproducible manual column, gone
# upstream); step-04 statistics added so the table stands alone as the supplementary
# annotation table.
cols <- c("trait","gene_id","locus_id","lead_SNP","lead_SNP_class","dist_to_lead_bp",
          "n_haplotype_groups","kw_p_raw","fdr_q","eta_squared","delta_top_bottom_sd",
          "final_call","annotation_source","annotation_from_check",
          "n_sources_with_call","needs_review_two_calls",
          "call_check1_swissprot","call_check2_interpro","call_check3_nr",
          "sp_vs_interpro_agreement",
          "sp_name","sp_accession","sp_organism","sp_pident","sp_qcovhsp","sp_evalue","sp_status",
          "interpro_status","n_interpro","interpro_ids","interpro_descs",
          "n_pfam","pfam_ids","pfam_descs","n_go","go_terms",
          "nr_title","nr_accession","nr_pident","nr_qcovs","nr_evalue","nr_status")
miss <- setdiff(cols, names(m))
if (length(miss)) stop("Columns absent after the merge: ", paste(miss, collapse=", "))
m <- m[, cols]
m <- m[order(m$trait, m$gene_id), ]
write.table(m, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- report ----
cat("=== Step 08: annotation master table ===\n")
cat("rows:", nrow(m), " (expect ", n_expected, ")   unique genes:", length(unique(m$gene_id)), "\n\n")
cat("call_source:\n"); print(table(m$call_source))
cat("\nannotation_source:\n"); print(table(m$annotation_source))
cat("\nn_sources_with_call:\n"); print(table(m$n_sources_with_call))
cat("\nflagged for joint review (2 differing calls):", sum(m$needs_review_two_calls), "\n")
cat("\nsp_vs_interpro_agreement:\n"); print(table(m$sp_vs_interpro_agreement))
cat("\nlead_SNP_class:\n"); print(table(m$lead_SNP_class))

cat("\nWrote:", out, "\n")
