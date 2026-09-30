# =============================================================================
# 06_removed_sites_table.R -- every site the `shared_sites` version drops, one row each
# =============================================================================
# `Table_site_overlap.tsv` says HOW MANY sites were removed and how many of them were
# monomorphic. This script says WHICH ones, and whether each still carried variation,
# so the cost of the shared-sites rule can be inspected site by site rather than trusted.
#
# WHAT COUNTS AS REMOVED
#   The shared set is keyed on CHROM:POS:REF:ALT, so a site is dropped when the other
#   file has no record at that position at all, or has one with a different alternate
#   allele. Three kinds, named in the `removed_from` column:
#     wild_only    in the wild call set, no record at that position in the elite export
#     elite_only   in the elite export, no record at that position in the wild call set
#     triallelic   both files have the position and the same REF, but different ALT.
#                  It is dropped from BOTH files, so it appears twice, once per file.
#
# WHAT "POLYMORPHIC" MEANS HERE  (the same two reference sets as 03_build_matrices.R)
#   wild sites  -- judged only over the accessions crosshap ASSIGNED to a haplotype group
#                  for that gene. Unassigned accessions are dropped from the test and from
#                  the figure, so variation confined to them could not have been drawn.
#                  The set differs per gene (GPAT6 123, GH17 107, PHT4;3 205).
#   elite sites -- judged only over the 5 CONFIGURED lines, not the 136-line pool, because
#                  those five are the only rows drawn.
#   A site is monomorphic if it has at least one call and every call is the same genotype,
#   and `undetermined` if nothing was called at all.
#
# HOW TO READ IT
#   The `cost` column is the point of the table:
#     "wild variation lost"   a polymorphic wild site that could not be drawn
#     "no information lost"   a removed site that was monomorphic anyway
#   Summed per gene, these reproduce the monomorphic counts in Table_site_overlap.tsv.
#
# Output: results/tables/Table_removed_sites.tsv        one row per removed site x file
#         results/tables/Table_removed_sites_summary.tsv one row per gene x removed_from
#
# Created 2026-09-23 (requested for the manuscript: which sites the shared-sites rule
# costs, kept so it can be re-examined later). Reads only; changes no existing output.
# =============================================================================

suppressPackageStartupMessages({ library(vcfR); library(dplyr) })

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(script_dir, "_load_params.R"))
P <- load_params(file.path(script_dir, "..", "config", "params.sh"))
Sys.setenv(TMPDIR = P$TMPDIR)

genes <- target_genes_df(P)
elite <- read_elite_lines(P)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# Genotype coding and VCF reader, identical to 03_build_matrices.R.
gt_code <- function(gt) {
  out <- rep(NA_real_, length(gt))
  out[gt %in% c("0/0", "0|0")] <- 0
  out[gt %in% c("1/1", "1|1")] <- 1
  out[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- 2
  out
}
read_gt <- function(path) {
  v   <- suppressMessages(vcfR::read.vcfR(path, verbose = FALSE))
  gt  <- vcfR::extract.gt(v, element = "GT", as.numeric = FALSE)
  fix <- as.data.frame(vcfR::getFIX(v), stringsAsFactors = FALSE)
  key <- paste(fix$CHROM, fix$POS, fix$REF, fix$ALT, sep = ":")
  m <- matrix(as.numeric(apply(gt, 2, gt_code)), nrow = nrow(gt),
              dimnames = list(key, colnames(gt)))
  list(mat = t(m), fix = fix, key = key)
}
GT_NAME <- c("0" = "REF", "1" = "ALT", "2" = "HET")

# Genotype tally of one site over a chosen set of rows.
tally_site <- function(mat, key, rows) {
  x <- mat[rows, key]
  n_called <- sum(!is.na(x))
  list(n_called   = n_called,
       n_ref      = sum(x == 0, na.rm = TRUE),
       n_alt      = sum(x == 1, na.rm = TRUE),
       n_het      = sum(x == 2, na.rm = TRUE),
       n_missing  = sum(is.na(x)),
       n_genotypes = length(unique(x[!is.na(x)])),
       state      = if (n_called == 0L) "undetermined"
                    else if (length(unique(x[!is.na(x)])) == 1L) "monomorphic" else "polymorphic")
}

rows_out <- list()
for (i in seq_len(nrow(genes))) {
  gid <- genes$gene_id[i]; short <- genes$short_name[i]

  gw  <- read.delim(P$STEP04_GENE_WINDOWS, stringsAsFactors = FALSE)
  gwr <- gw[gw$gene_id == gid, ][1, ]; trait <- gwr$trait

  W <- read_gt(file.path(P$STEP04_RAW_VCF_DIR, trait, paste0(gid, ".vcf.gz")))

  ho  <- readRDS(file.path(P$STEP04_CACHE, trait, gid,
                           paste0("MGmin_", P$CROSSHAP_MGMIN), "HapObject.rds"))
  lab <- paste0("Haplotypes_MGmin", P$CROSSHAP_MGMIN, "_E", P$CROSSHAP_EPSILON)
  ind <- ho$HapObject[[lab]]$Indfile
  ind$hap <- as.character(ind$hap); ind <- ind[ind$hap != "0", ]
  assigned <- intersect(ind$Ind, rownames(W$mat))

  E <- read_gt(file.path(P$DIR_ELITE_TRIMMED, paste0(gid, ".pool.vcf.gz")))
  keep <- elite$sample_id[elite$sample_id %in% rownames(E$mat)]
  E$mat <- E$mat[keep, , drop = FALSE]
  rownames(E$mat) <- elite$line_name[match(keep, elite$sample_id)]

  shared_keys <- intersect(colnames(W$mat), colnames(E$mat))
  wpos <- as.integer(W$fix$POS); epos <- as.integer(E$fix$POS)
  tri_pos <- intersect(wpos, epos)[
    vapply(intersect(wpos, epos), function(pp) {
      w <- W$fix[wpos == pp, ][1, ]; e <- E$fix[epos == pp, ][1, ]
      w$REF == e$REF && w$ALT != e$ALT
    }, logical(1))]

  for (side in c("wild", "elite")) {
    M    <- if (side == "wild") W else E
    rows <- if (side == "wild") assigned else rownames(E$mat)
    pos  <- as.integer(M$fix$POS)
    drop <- setdiff(colnames(M$mat), shared_keys)
    for (k in drop) {
      pp   <- pos[match(k, M$key)]
      kind <- if (pp %in% tri_pos) "triallelic" else paste0(side, "_only")
      fx   <- M$fix[match(k, M$key), ]
      tl   <- tally_site(M$mat, k, rows)
      rows_out[[length(rows_out) + 1L]] <- data.frame(
        gene_id = gid, short_name = short, trait = trait,
        file = side, removed_from = kind,
        chr = fx$CHROM, pos = pp, ref = fx$REF, alt = fx$ALT,
        is_indel = nchar(fx$REF) > 1 || any(nchar(strsplit(fx$ALT, ",")[[1]]) > 1),
        n_rows_judged_on = length(rows),
        n_called = tl$n_called, n_REF = tl$n_ref, n_ALT = tl$n_alt, n_HET = tl$n_het,
        n_missing = tl$n_missing, n_distinct_genotypes = tl$n_genotypes,
        state = tl$state,
        cost = if (side == "wild" && tl$state == "polymorphic") "wild variation lost"
               else if (tl$state == "polymorphic") "elite variation lost"
               else if (tl$state == "monomorphic") "no information lost"
               else "nothing called",
        stringsAsFactors = FALSE)
    }
  }
  msg(short, ": ", sum(vapply(rows_out, function(r) r$gene_id == gid, logical(1))),
      " removed site-rows (", length(assigned), " assigned accessions, ",
      nrow(E$mat), " elite lines)")
}

out <- bind_rows(rows_out) %>% arrange(short_name, file, pos)
write.table(out, file.path(P$DIR_TABLES, "Table_removed_sites.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

summ <- out %>%
  group_by(short_name, trait, file, removed_from) %>%
  summarise(n_sites = dplyr::n(),
            n_polymorphic = sum(state == "polymorphic"),
            n_monomorphic = sum(state == "monomorphic"),
            n_undetermined = sum(state == "undetermined"),
            n_indel = sum(is_indel), .groups = "drop") %>%
  arrange(short_name, file, removed_from)
write.table(summ, file.path(P$DIR_TABLES, "Table_removed_sites_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\n"); print(as.data.frame(summ))
cat(sprintf("\n[06] %d removed site-rows -> %s\n", nrow(out),
            file.path(P$DIR_TABLES, "Table_removed_sites.tsv")))
cat(sprintf("[06] summary                -> %s\n",
            file.path(P$DIR_TABLES, "Table_removed_sites_summary.tsv")))
