#!/usr/bin/env Rscript
# ============================================================================
# 03_build_tables.R
#
# WHAT   Parse the intersect output, attach locus metadata, compute gene->lead
#        distances, and write every output table. No plots (tables only).
#
# READS  intermediates/loci.tsv, intersect_raw.tsv, genes_per_locus.tsv
#        $GWAS_PARAMS_TSV  (step-01 parameters, folded into the methods table)
#
# WRITES results/tables/
#   candidate_genes.tsv      MAIN. one row per gene x locus. Column names and
#                            types are the CONTRACT with step 04 -- see below.
#   lead_loci.tsv            one row per locus + its gene counts. Also part of
#                            the step-04 contract.
#   loci_without_genes.tsv   loci whose interval contains no annotated gene
#   genes_with_annotation.tsv  the subset carrying a functional description
#   per_trait_summary.tsv    per-trait counts
#   analysis_parameters.tsv  every constant used, step 01 + step 03
#   results_chapter_numbers.txt  paste-ready sentences for the manuscript
#
# STEP-04 CONTRACT (04_.../01_scripts/00_build_gene_windows.R):
#   candidate_genes.tsv MUST have: trait, locus_id, lead_SNP, class, gene_id,
#                                  chr, gene_start, gene_end, strand,
#                                  dist_to_lead_bp, description
#   lead_loci.tsv       MUST have: trait, locus_id, lead_pos, lead_neg_log10p
#   Both filenames and both column spellings are load-bearing. lead_pos and
#   lead_neg_log10p keep step-04's spelling, NOT step-01's lead_bp /
#   lead_neg_log10_p -- do not "tidy" them.
#
# dist_to_lead_bp  signed nearest-edge distance from the lead SNP to the gene:
#   0  lead SNP falls inside the gene body
#   <0 gene is upstream of the lead (lower coordinate)
#   >0 gene is downstream of the lead (higher coordinate)
# ============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_load_params.R"))
p <- load_params()

base   <- p$STEP03_BASE
inter  <- file.path(base, "intermediates")
outdir <- file.path(base, "results", "tables")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("=== 03_build_tables.R ===\n")

loci   <- read.delim(file.path(inter, "loci.tsv"), stringsAsFactors = FALSE)
# Loci RECEIVED from step 01, read from the frozen snapshot rather than assumed.
# (Was hard-coded "34" - a leftover from an earlier locus set that made
# analysis_parameters.tsv claim 34 received / 36 searched, implying 2 were invented.)
snap_path <- file.path(p$STEP03_BASE, "inputs", "loci_handoff_snapshot.tsv")
.snap <- if (file.exists(snap_path))
  read.delim(snap_path, comment.char = "#", stringsAsFactors = FALSE) else NULL
n_handoff_rows   <- if (is.null(.snap)) NA_integer_ else nrow(.snap)
n_sub_in_handoff <- if (is.null(.snap) || !"class" %in% names(.snap)) NA_integer_ else
  sum(.snap$class == "subthreshold_peak")
counts <- read.delim(file.path(inter, "genes_per_locus.tsv"), stringsAsFactors = FALSE)
# Genes in the filtered annotation, counted rather than asserted (was a literal 35106L,
# which would go stale silently if the GFF were ever updated).
gff_path <- file.path(inter, "genes_1to7H.gff")
n_annotation_genes <- if (file.exists(gff_path))
  as.integer(length(readLines(gff_path, warn = FALSE))) else NA_integer_

raw_path <- file.path(inter, "intersect_raw.tsv")
has_hits <- file.exists(raw_path) && file.info(raw_path)$size > 0
raw <- if (has_hits) read.delim(raw_path, header = FALSE, stringsAsFactors = FALSE) else
       data.frame(matrix(nrow = 0, ncol = 13))
stopifnot(ncol(raw) == 13)
colnames(raw) <- c("int_chr","int_start0","int_end","locus_id",
                   "g_chr","g_source","g_feature","g_start","g_end",
                   "g_score","g_strand","g_frame","g_attr")

## ---- GFF attribute parsing --------------------------------------------------
# Attributes are ';'-separated key=value. A ';' inside a value is escaped %3B,
# so splitting on ';' is safe, and the value is percent-decoded afterwards.
get_attr <- function(attr, key) {
  out <- rep(NA_character_, length(attr))
  m   <- regexpr(paste0("(?:^|;)", key, "=([^;]*)"), attr, perl = TRUE)
  hit <- m != -1L
  if (any(hit)) {
    vals <- regmatches(attr, m)
    out[hit] <- sub(paste0("^.*?", key, "="), "", vals)
  }
  out
}
decode <- function(x) vapply(x, function(v) if (is.na(v)) NA_character_ else utils::URLdecode(v),
                             character(1), USE.NAMES = FALSE)

raw$gene_id     <- get_attr(raw$g_attr, "gene_id")
raw$biotype     <- get_attr(raw$g_attr, "biotype")
raw$description <- decode(get_attr(raw$g_attr, "description"))

## ---- Attach locus metadata --------------------------------------------------
idx <- match(raw$locus_id, loci$locus_id)
if (anyNA(idx)) stop("intersect_raw.tsv carries a locus_id absent from loci.tsv")
for (f in c("trait","class","lead_SNP","lead_bp","lead_neg_log10_p","lead_MAF",
            "n_member_SNPs","locus_start","locus_end","span_kb","source_search_id")) {
  raw[[f]] <- loci[[f]][idx]
}

## ---- Signed distance and in-span flag ---------------------------------------
lp <- raw$lead_bp; gs <- raw$g_start; ge <- raw$g_end
raw$dist_to_lead_bp <- as.integer(ifelse(lp >= gs & lp <= ge, 0L,
                                  ifelse(ge < lp, ge - lp, gs - lp)))
raw$lead_inside_gene <- (lp >= gs & lp <= ge)
# With FLANK_BP = 0 the interval IS the locus span, so in_locus_span is TRUE for
# every row. It is kept so the column stays meaningful if a flank is switched on.
raw$in_locus_span <- (ge >= raw$locus_start) & (gs <= raw$locus_end)

## ---- Table 1: candidate_genes.tsv (MAIN) ------------------------------------
cand <- data.frame(
  trait            = raw$trait,
  locus_id         = raw$locus_id,
  class            = raw$class,
  lead_SNP         = raw$lead_SNP,
  gene_id          = raw$gene_id,
  biotype          = raw$biotype,
  chr              = raw$g_chr,
  gene_start       = raw$g_start,
  gene_end         = raw$g_end,
  strand           = raw$g_strand,
  gene_length_bp   = raw$g_end - raw$g_start + 1L,
  dist_to_lead_bp  = raw$dist_to_lead_bp,
  lead_inside_gene = raw$lead_inside_gene,
  in_locus_span    = raw$in_locus_span,
  has_annotation   = !is.na(raw$description),
  description      = raw$description,
  stringsAsFactors = FALSE
)
cand <- cand[order(cand$trait, cand$chr, cand$gene_start), ]
rownames(cand) <- NULL
write.table(cand, file.path(outdir, "candidate_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Guard the step-04 contract explicitly rather than discovering a break there.
need04 <- c("trait","locus_id","lead_SNP","class","gene_id","chr",
            "gene_start","gene_end","strand","dist_to_lead_bp","description")
if (!all(need04 %in% names(cand))) stop("candidate_genes.tsv breaks the step-04 contract.")

## ---- Table 2: lead_loci.tsv -------------------------------------------------
cn <- setNames(counts$n_genes, counts$locus_id)
loci$n_genes <- as.integer(cn[loci$locus_id]); loci$n_genes[is.na(loci$n_genes)] <- 0L
tab <- function(v) { t <- tapply(v, cand$locus_id, sum); x <- as.integer(t[loci$locus_id]); x[is.na(x)] <- 0L; x }
loci$n_protein_coding   <- if (nrow(cand)) tab(cand$biotype == "protein_coding") else 0L
loci$n_with_annotation  <- if (nrow(cand)) tab(cand$has_annotation)              else 0L
loci$n_genes_containing_lead <- if (nrow(cand)) tab(cand$lead_inside_gene)       else 0L

# Step-04 spelling: lead_pos / lead_neg_log10p.
loci$lead_pos        <- loci$lead_bp
loci$lead_neg_log10p <- loci$lead_neg_log10_p

lead_cols <- c("trait","locus_id","source_search_id","chr","class",
               "lead_SNP","lead_pos","lead_neg_log10p","lead_MAF",
               "n_member_SNPs","locus_start","locus_end","span_kb",
               "flank_bp","search_start","search_end","search_kb",
               "n_genes","n_protein_coding","n_with_annotation","n_genes_containing_lead")
write.table(loci[, lead_cols], file.path(outdir, "lead_loci.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (!all(c("trait","locus_id","lead_pos","lead_neg_log10p") %in% names(loci)))
  stop("lead_loci.tsv breaks the step-04 contract.")

## ---- Table 3: loci_without_genes.tsv ----------------------------------------
# A deliberate, reportable result under FLANK_BP=0, not missing data.
empty <- loci[loci$n_genes == 0, c("trait","locus_id","chr","class","lead_SNP","lead_pos",
                                   "lead_neg_log10p","n_member_SNPs","locus_start",
                                   "locus_end","span_kb","search_kb")]
empty <- empty[order(empty$trait, empty$chr, empty$lead_pos), ]
write.table(empty, file.path(outdir, "loci_without_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- Table 4: genes_with_annotation.tsv -------------------------------------
ann <- cand[cand$has_annotation, c("trait","locus_id","class","lead_SNP","gene_id","chr",
                                   "gene_start","gene_end","strand","dist_to_lead_bp",
                                   "lead_inside_gene","description")]
ann <- ann[order(ann$trait, abs(ann$dist_to_lead_bp)), ]
write.table(ann, file.path(outdir, "genes_with_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- Table 5: per_trait_summary.tsv -----------------------------------------
traits <- sort(unique(loci$trait))
summ <- do.call(rbind, lapply(traits, function(tr) {
  lt <- loci[loci$trait == tr, ]; ct <- cand[cand$trait == tr, ]
  data.frame(
    trait                  = tr,
    n_loci                 = nrow(lt),
    n_loci_significant     = sum(lt$class == "significant_locus"),
    n_loci_subthreshold    = sum(lt$class == "subthreshold_peak"),
    n_loci_with_genes      = sum(lt$n_genes > 0),
    n_loci_without_genes   = sum(lt$n_genes == 0),
    total_search_kb        = round(sum(lt$search_kb), 1),
    n_genes                = nrow(ct),
    n_unique_genes         = length(unique(ct$gene_id)),
    n_protein_coding       = sum(ct$biotype == "protein_coding"),
    n_with_annotation      = sum(ct$has_annotation),
    n_genes_containing_lead= sum(ct$lead_inside_gene),
    stringsAsFactors = FALSE)
}))
write.table(summ, file.path(outdir, "per_trait_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- Table 6: analysis_parameters.tsv ---------------------------------------
gp <- if (file.exists(p$GWAS_PARAMS_TSV))
        read.delim(p$GWAS_PARAMS_TSV, stringsAsFactors = FALSE) else
        data.frame(parameter = character(), value = character())
gp$stage <- "01_GWAS_and_loci"
s3 <- data.frame(stage = "03_candidate_genes", parameter = c(
  "loci_handoff_table","n_loci_received","n_loci_searched",
  "subthreshold_peaks_in_handoff","subthreshold_peaks_searched",
  "flank_bp_each_side","search_interval_rule",
  "total_search_space_Mb","annotation_source","gene_feature_used",
  "chromosomes_searched","n_genes_in_annotation","bedtools_version",
  "n_gene_locus_rows","n_unique_genes","n_loci_without_genes","date_run"),
  value = c(
  basename(p$LOCI_HANDOFF_TSV), as.character(n_handoff_rows), as.character(nrow(loci)),
  # Report what the handoff actually CONTAINED and what was searched, not the flag
  # setting. INCLUDE_SUBTHRESHOLD_PEAKS="yes" printed as if peaks were included even
  # when step 01 supplied none, which reads as a scientific claim rather than a switch.
  as.character(n_sub_in_handoff), as.character(sum(loci$class == "subthreshold_peak")),
  format(p$FLANK_BP, scientific = FALSE),
  if (p$FLANK_BP == 0) "search interval = step-01 LD locus span, no flank added"
    else sprintf("search interval = step-01 LD locus span +/- %s bp", format(p$FLANK_BP, scientific=FALSE)),
  sprintf("%.2f", sum(loci$search_end - loci$search_start + 1)/1e6),
  "Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3.gz (Ensembl Plants r62)",
  "GFF3 col3 == gene", "1H-7H (unplaced CAJHDD* scaffolds excluded)",
  format(n_annotation_genes, big.mark = ","),
  tryCatch(sub("^bedtools ", "", system2(p$BEDTOOLS, "--version", stdout = TRUE)), error = function(e) "NA"),
  as.character(nrow(cand)), as.character(length(unique(cand$gene_id))),
  as.character(sum(loci$n_genes == 0)), format(Sys.Date(), "%Y-%m-%d")),
  stringsAsFactors = FALSE)
allp <- rbind(gp[, c("stage","parameter","value")], s3)
write.table(allp, file.path(outdir, "analysis_parameters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- Table 7: results_chapter_numbers.txt -----------------------------------
con <- file(file.path(outdir, "results_chapter_numbers.txt"), "w"); wl <- function(...) writeLines(paste0(...), con)
wl("================================================================")
wl("  Candidate genes at the GWAS loci -- numbers for the manuscript")
wl("  generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by 03_build_tables.R")
wl("================================================================")
wl("")
wl("METHOD")
wl(sprintf("  Search interval = the step-01 LD locus span%s.",
   if (p$FLANK_BP == 0) ", with no flanking window added"
   else sprintf(" extended by %s bp on each side", format(p$FLANK_BP, big.mark=","))))
# The step-01 clumping parameters are READ from step 01's own table, never restated
# here. They have changed three times; a literal "+/-1 Mb" survived a change to 2 Mb and
# sat inside this paste-ready manuscript sentence, which is the one place a wrong number
# would reach a reviewer.
gp_lookup <- function(key, default = "?") {
  if (!nrow(gp)) return(default)
  v <- gp$value[gp$parameter == key]
  if (length(v) != 1 || is.na(v)) default else as.character(v)
}
clump_mb <- suppressWarnings(as.numeric(gp_lookup("clump_kb", NA)) / 1000)
wl(sprintf("  Loci were defined in step 01 by PLINK --clump (lead p < %s, members by",
           gp_lookup("clump_p1_lead_threshold", "9.008e-07")))
wl(sprintf("  LD only at r2 >= %s within +/-%s Mb), each locus severed at the first internal",
           gp_lookup("clump_r2", "0.5"),
           if (is.na(clump_mb)) "?" else format(clump_mb, trim = TRUE)))
wl(sprintf("  gap > %s kb. Genes = GFF3 'gene' features on 1H-7H of Morex V3 (Ensembl Plants r62).",
           gp_lookup("gap_rule_kb", "?")))
wl("")
wl("HEADLINE NUMBERS")
wl(sprintf("  Loci searched          : %d (%d significant, %d sub-threshold protein peaks)",
   nrow(loci), sum(loci$class=="significant_locus"), sum(loci$class=="subthreshold_peak")))
wl(sprintf("  Total search space     : %.2f Mb", sum(loci$search_end-loci$search_start+1)/1e6))
wl(sprintf("  Loci containing genes  : %d", sum(loci$n_genes>0)))
wl(sprintf("  Loci with no gene      : %d", sum(loci$n_genes==0)))
wl(sprintf("  Candidate genes        : %d (%d unique gene IDs)", nrow(cand), length(unique(cand$gene_id))))
wl(sprintf("  Protein-coding         : %d", sum(cand$biotype=="protein_coding")))
wl(sprintf("  With a functional description : %d (%.1f%%)",
   sum(cand$has_annotation), 100*sum(cand$has_annotation)/max(1,nrow(cand))))
wl(sprintf("  Lead SNP inside a gene body   : %d locus/gene pair(s)", sum(cand$lead_inside_gene)))
wl("")
wl("PASTE-READY SENTENCE")
wl(sprintf("  The %d loci were searched for annotated genes across %.2f Mb of the Morex V3",
   nrow(loci), sum(loci$search_end-loci$search_start+1)/1e6))
wl(sprintf("  assembly, yielding %d candidate genes at %d loci; the remaining %d loci contained",
   nrow(cand), sum(loci$n_genes>0), sum(loci$n_genes==0)))
wl(sprintf("  no annotated gene. %d of the candidates carry a functional description.",
   sum(cand$has_annotation)))
wl("")
wl("PER TRAIT")
for (i in seq_len(nrow(summ))) with(summ[i,], wl(sprintf(
  "  %-11s %2d loci (%d with genes, %d empty), %3d genes, %2d annotated",
  trait, n_loci, n_loci_with_genes, n_loci_without_genes, n_genes, n_with_annotation)))
wl("")
wl("PER LOCUS")
for (tr in traits) {
  wl(sprintf("  --- %s ---", tr))
  lt <- loci[loci$trait == tr, ]; lt <- lt[order(lt$chr, lt$lead_pos), ]
  for (i in seq_len(nrow(lt))) {
    wl(sprintf("  %-14s %s  lead %-14s -log10p=%.4f  span=%8.1f kb  %2d SNPs  %2d gene(s)%s",
       lt$locus_id[i], lt$chr[i], lt$lead_SNP[i], lt$lead_neg_log10p[i],
       lt$span_kb[i], lt$n_member_SNPs[i], lt$n_genes[i],
       if (lt$class[i]=="subthreshold_peak") "  [sub-threshold]" else ""))
    g <- cand[cand$locus_id == lt$locus_id[i], ]
    if (nrow(g)) for (j in seq_len(nrow(g)))
      wl(sprintf("                   %s  %+8d bp  %s", g$gene_id[j], g$dist_to_lead_bp[j],
                 if (is.na(g$description[j])) "(no annotation)" else substr(g$description[j],1,90)))
  }
}
wl("")
wl("LOCI WITH NO ANNOTATED GENE")
if (nrow(empty)) for (i in seq_len(nrow(empty))) with(empty[i,], wl(sprintf(
  "  %-14s %s  lead %-14s -log10p=%.4f  span=%.1f kb", locus_id, chr, lead_SNP, lead_neg_log10p, span_kb)))
close(con)

## ---- Console report ----------------------------------------------------------
cat(sprintf("\ncandidate_genes.tsv      : %d rows (%d unique genes)\n", nrow(cand), length(unique(cand$gene_id))))
cat(sprintf("lead_loci.tsv            : %d loci\n", nrow(loci)))
cat(sprintf("loci_without_genes.tsv   : %d loci\n", nrow(empty)))
cat(sprintf("genes_with_annotation.tsv: %d genes\n", nrow(ann)))
cat("\nper_trait_summary.tsv:\n"); print(summ[, c("trait","n_loci","n_loci_with_genes",
  "n_loci_without_genes","n_genes","n_with_annotation")], row.names = FALSE)
cat("\nWrote 7 files to results/tables/\n")
