# =============================================================================
# 03_build_matrices.R -- wild haplotype consensus + elite genotypes, aligned
# =============================================================================
# Turns the two VCFs into the aligned genotype matrices the figures draw, in the
# THREE SNP-matching treatments defined in config/params.sh. No plotting here.
#
# WILD SIDE -- consumed, never recomputed.
#   Haplotype group membership comes from the cached crosshap HapObject of step
#   04 run `loci_LDspan_eps06_V4` (Indfile: Ind, hap, Pheno). crosshap is NOT
#   re-run and the elite lines are NEVER shown to it. Accessions in hap "0"
#   (unassigned) are dropped, exactly as step 04's Kruskal-Wallis does.
#
#   Each group is reduced to ONE row by PER-SNP MAJORITY CONSENSUS over its
#   members, ignoring NA. Rationale: wild missingness is high (19% of calls at
#   PHT4;3), so a single representative accession would carry grey tiles that are
#   artefacts of that one plant's sequencing rather than features of the
#   haplotype. Ties (exactly 50/50) resolve to REFERENCE when REF is one of the
#   tied states (changed 2026-09-24, user decision; until then ties resolved to NA
#   and were drawn as missing). In the current data this affects exactly two GH17
#   cells: group A at 5H:462728332 (18 REF / 18 ALT) and group E at 5H:462729123
#   (5 / 5). A tie that does not involve REF (ALT vs HET only) still resolves to NA;
#   it does not occur in the current data.
#   The closest real accession to each consensus is still reported, in
#   results/tables/consensus_representatives.tsv, for verification only.
#
# ELITE SIDE -- genotypes only, no group assignment, no phenotype. Read from the
#   candidate-pool VCF of step 01 and subset to the lines in config/elite_lines.tsv
#   (changed 2026-09-14; the call-rate screen in step 02 reads the same file).
#
# ALLELE CONCORDANCE (hard check, requested explicitly).
#   Both call sets are called against Morex V3 and the wild set was produced
#   without PLINK allele-flipping, so at a shared position REF and ALT must be
#   IDENTICAL. If they were ever swapped, 0/0 in one file and 0/0 in the other
#   would denote OPPOSITE alleles and every barcode here would be silently
#   inverted. Every shared position is classified identical / swapped /
#   alt_differs / ref_differs into results/tables/allele_concordance.tsv.
#   `swapped` and `ref_differs` ABORT the run (ALLELE_MISMATCH_ACTION="stop").
#   `alt_differs` is a TRIALLELIC site -- same reference base, different alternate
#   allele segregating in the wild panel and in the elite panel. That is biology,
#   not a coding error, so it does not abort; the site simply never enters the
#   shared set (the key is CHROM:POS:REF:ALT) and is counted in the tables.
#   Measured 2026-09-10: 0 swapped anywhere, 2 triallelic of 90 shared positions.
#
# THE THREE VERSIONS (SNPs differ between the files in BOTH directions):
#   shared_sites   only SNPs present in wild AND elite. No assumption.
#   filled_marked  all wild SNPs; elite tiles with no elite record get their own
#                  "no record" state -- kept visible rather than assumed REF.
#   filled_silent  all wild SNPs; missing elite records become REF, the earlier
#                  behaviour (bcftools merge -0). Kept only for comparison.
#   In BOTH filled versions the TRIALLELIC columns (same REF, different ALT; see
#   ALLELE CONCORDANCE) carry the elite line's REAL genotype, never a fill:
#   0/0 -> REF (same base as the wild REF, comparable); a call containing the
#   elite-only third allele -> Triallelic. Added 2026-09-14. What each configured
#   line carries there is written to triallelic_sites_elite_genotypes.tsv.
#   Elite-only SNPs are never drawn -- the wild groups have no data there, so
#   they would add columns of grey to every group row. They are counted and
#   listed in results/tables/site_overlap_summary.tsv instead.
#
# Output -> intermediates/matrices/<gene>__<version>.rds   (matrix + metadata)
#           results/tables/allele_concordance.tsv
#           results/tables/site_overlap_summary.tsv
#             (2026-09-23: + monomorphic-site counts per file. "Monomorphic" is judged
#              only over the rows that enter the comparison -- for the wild file the
#              crosshap-ASSIGNED accessions of that gene, for the elite file the 5
#              configured lines -- and is reported over all sites, the sites kept in
#              shared_sites, and the sites removed from it.)
#           results/tables/consensus_representatives.tsv
#           results/tables/elite_genotypes_long.tsv
#           results/tables/triallelic_sites_elite_genotypes.tsv
#
# Author : Shahar Liviatan
# Created: 2026-09-10
# =============================================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
})

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(script_dir, "_load_params.R"))
P <- load_params(file.path(script_dir, "..", "config", "params.sh"))
Sys.setenv(TMPDIR = P$TMPDIR)

genes <- target_genes_df(P)
elite <- read_elite_lines(P)

dir.create(P$DIR_MATRICES, recursive = TRUE, showWarnings = FALSE)
dir.create(P$DIR_TABLES,   recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---------------------------------------------------------------------------
# Genotype coding. 0 = REF, 1 = ALT, 2 = HET, NA = no call.
# Step 04's gt_to_bin01() folded het into NA; here het is its own state because
# the elite panel carries het calls and the wild set carries none (verified).
# ---------------------------------------------------------------------------
gt_code <- function(gt) {
  out <- rep(NA_real_, length(gt))
  out[gt %in% c("0/0", "0|0")] <- 0
  out[gt %in% c("1/1", "1|1")] <- 1
  out[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- 2
  out
}
CODE_REF <- 0; CODE_ALT <- 1; CODE_HET <- 2; CODE_NORECORD <- 3; CODE_TRIALLELIC <- 4

read_gt <- function(path) {
  v   <- suppressMessages(vcfR::read.vcfR(path, verbose = FALSE))
  gt  <- vcfR::extract.gt(v, element = "GT", as.numeric = FALSE)
  fix <- as.data.frame(vcfR::getFIX(v), stringsAsFactors = FALSE)
  key <- paste(fix$CHROM, fix$POS, fix$REF, fix$ALT, sep = ":")
  m <- apply(gt, 2, gt_code)
  m <- matrix(as.numeric(m), nrow = nrow(gt), dimnames = list(key, colnames(gt)))
  list(mat = t(m), fix = fix, key = key)   # rows = samples, cols = sites
}

# Per-SNP majority over rows, ignoring NA; exact ties -> REF if REF is among the
# tied states, otherwise NA (tie rule changed 2026-09-24, user decision; was: all ties -> NA).
consensus_row <- function(m) {
  apply(m, 2, function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_real_)
    tb <- table(x)
    top <- names(tb)[tb == max(tb)]
    if (length(top) > 1) return(if (as.character(CODE_REF) %in% top) CODE_REF else NA_real_)
    as.numeric(top)
  })
}

frac_disagree <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (!sum(ok)) return(NA_real_)
  sum(a[ok] != b[ok]) / sum(ok)
}

allele_tbl <- list(); overlap_tbl <- list(); rep_tbl <- list(); elite_long <- list(); tri_tbl <- list()

for (i in seq_len(nrow(genes))) {
  gid <- genes$gene_id[i]; short <- genes$short_name[i]
  msg("=== ", short, " (", gid, ")")

  # ---- wild: raw per-gene VCF of step 04 ----------------------------------
  gw <- read.delim(P$STEP04_GENE_WINDOWS, stringsAsFactors = FALSE)
  gwr <- gw[gw$gene_id == gid, ][1, ]
  trait <- gwr$trait
  wild_vcf <- file.path(P$STEP04_RAW_VCF_DIR, trait, paste0(gid, ".vcf.gz"))
  stopifnot(file.exists(wild_vcf))
  W <- read_gt(wild_vcf)

  # ---- wild: haplotype groups from the cached HapObject --------------------
  cache <- file.path(P$STEP04_CACHE, trait, gid,
                     paste0("MGmin_", P$CROSSHAP_MGMIN), "HapObject.rds")
  stopifnot(file.exists(cache))
  ho <- readRDS(cache)
  lab <- paste0("Haplotypes_MGmin", P$CROSSHAP_MGMIN, "_E", P$CROSSHAP_EPSILON)
  stopifnot(lab %in% names(ho$HapObject))
  ind <- ho$HapObject[[lab]]$Indfile
  ind$hap <- as.character(ind$hap)
  ind <- ind[ind$hap != "0", ]                     # drop unassigned, as step 04 does

  # ---- elite: candidate-pool VCF from 01_fetch_elite_vcfs.sh, subset below to
  #      the configured lines (screen and figures read identical data) ----------
  elite_vcf <- file.path(P$DIR_ELITE_TRIMMED, paste0(gid, ".pool.vcf.gz"))
  stopifnot(file.exists(elite_vcf))
  E <- read_gt(elite_vcf)
  # SAMEA ids -> cultivar names, in the order given by config/elite_lines.tsv
  keep <- elite$sample_id[elite$sample_id %in% rownames(E$mat)]
  missing_lines <- setdiff(elite$sample_id, rownames(E$mat))
  if (length(missing_lines))
    stop(gid, ": elite lines absent from the export: ", paste(missing_lines, collapse = ", "))
  E$mat <- E$mat[keep, , drop = FALSE]
  rownames(E$mat) <- elite$line_name[match(keep, elite$sample_id)]

  # ---- allele concordance on shared POSITIONS ------------------------------
  wf <- W$fix; ef <- E$fix
  wpos <- as.integer(wf$POS); epos <- as.integer(ef$POS)
  shared_pos <- intersect(wpos, epos)

  ac <- lapply(shared_pos, function(p) {
    wi <- which(wpos == p)[1]; ei <- which(epos == p)[1]
    wr <- wf$REF[wi]; wa <- wf$ALT[wi]; er <- ef$REF[ei]; ea <- ef$ALT[ei]
    # Four outcomes, and only two of them are errors:
    #   identical   -- the normal case; codes mean the same thing in both files
    #   swapped     -- REF/ALT inverted. FATAL: 0/0 would denote opposite alleles
    #   ref_differs -- different reference base at the same coordinate. FATAL
    #   alt_differs -- same REF, different ALT: a TRIALLELIC site where the wild
    #                  Levantine panel and the elite panel carry different
    #                  alternate alleles. Not an error -- biology. REF tiles stay
    #                  comparable but ALT tiles would denote different bases, so
    #                  the site is EXCLUDED from the shared set (the shared key is
    #                  CHROM:POS:REF:ALT, so this happens automatically) and counted.
    st <- if (wr == er && wa == ea) "identical"
          else if (wr == ea && wa == er) "swapped"
          else if (wr == er) "alt_differs"
          else "ref_differs"
    data.frame(gene_id = gid, short_name = short, chr = wf$CHROM[wi], pos = p,
               wild_ref = wr, wild_alt = wa, elite_ref = er, elite_alt = ea,
               status = st, stringsAsFactors = FALSE)
  })
  ac <- if (length(ac)) do.call(rbind, ac) else
        data.frame(gene_id = character(), short_name = character(), chr = character(),
                   pos = integer(), wild_ref = character(), wild_alt = character(),
                   elite_ref = character(), elite_alt = character(), status = character())
  allele_tbl[[gid]] <- ac

  n_ident <- sum(ac$status == "identical")
  n_swap  <- sum(ac$status == "swapped")
  n_altd  <- sum(ac$status == "alt_differs")
  n_refd  <- sum(ac$status == "ref_differs")
  msg("  allele concordance on ", nrow(ac), " shared positions: ",
      n_ident, " identical, ", n_swap, " swapped, ", n_altd,
      " triallelic (ALT differs, excluded), ", n_refd, " different REF")
  if (n_altd > 0) print(ac[ac$status == "alt_differs",
                           c("chr","pos","wild_ref","wild_alt","elite_ref","elite_alt")])
  # Only a swap or a differing REF is fatal: both mean a genotype code denotes a
  # different allele in the two files, which would silently invert the barcodes.
  if ((n_swap + n_refd) > 0 && identical(P$ALLELE_MISMATCH_ACTION, "stop")) {
    print(ac[ac$status %in% c("swapped", "ref_differs"), ])
    stop(gid, ": REF/ALT are not consistent between the wild and elite VCFs at ",
         n_swap + n_refd, " shared position(s) (", n_swap, " swapped, ", n_refd,
         " differing REF). Genotype codes are therefore not comparable. Set ",
         "ALLELE_MISMATCH_ACTION=warn in config/params.sh only if you have ",
         "decided how to handle it.")
  }

  # Only alleles-identical sites may be treated as shared. Because the key is
  # CHROM:POS:REF:ALT, triallelic sites drop out here without special handling.
  shared_keys <- intersect(colnames(W$mat), colnames(E$mat))

  # Triallelic sites: same position and REF, different ALT. Wild key -> elite key.
  tri <- ac[ac$status == "alt_differs", , drop = FALSE]
  tri_wild_keys  <- paste(tri$chr, tri$pos, tri$wild_ref,  tri$wild_alt,  sep = ":")
  tri_elite_keys <- paste(tri$chr, tri$pos, tri$elite_ref, tri$elite_alt, sep = ":")
  stopifnot(all(tri_wild_keys %in% colnames(W$mat)), all(tri_elite_keys %in% colnames(E$mat)))

  # wild-only / elite-only = NO record at that position in the other file.
  # (Triallelic sites have a record in both, so they are counted separately.)
  wild_only   <- setdiff(setdiff(colnames(W$mat), colnames(E$mat)), tri_wild_keys)
  elite_only  <- setdiff(setdiff(colnames(E$mat), colnames(W$mat)), tri_elite_keys)

  # What the configured elite lines actually carry at each triallelic site:
  #   0/0              -> REF        (same base as the wild REF, so comparable)
  #   1/1 or 0/1       -> Triallelic (the elite-only third allele is present)
  tri_code <- matrix(NA_real_, nrow = nrow(E$mat), ncol = length(tri_wild_keys),
                     dimnames = list(rownames(E$mat), tri_wild_keys))
  for (k in seq_along(tri_wild_keys)) {
    x <- E$mat[, tri_elite_keys[k]]
    tri_code[, k] <- ifelse(is.na(x), NA_real_,
                            ifelse(x == CODE_REF, CODE_REF, CODE_TRIALLELIC))
    tri_tbl[[paste(gid, k)]] <- data.frame(
      gene_id = gid, short_name = short, site = tri_wild_keys[k],
      wild_alleles  = paste0(tri$wild_ref[k],  "/", tri$wild_alt[k]),
      elite_alleles = paste0(tri$elite_ref[k], "/", tri$elite_alt[k]),
      line_name = rownames(E$mat),
      genotype = ifelse(is.na(x), "missing",
                 ifelse(x == CODE_REF, paste0("REF (", tri$elite_ref[k], ")"),
                 ifelse(x == CODE_HET, paste0("HET with third allele (", tri$elite_alt[k], ")"),
                        paste0("third allele (", tri$elite_alt[k], ")")))),
      stringsAsFactors = FALSE)
  }
  if (length(tri_wild_keys))
    msg("  triallelic sites: ", length(tri_wild_keys),
        " | calls carrying the elite-only third allele among configured lines: ",
        sum(tri_code == CODE_TRIALLELIC, na.rm = TRUE))

  e_is_indel <- nchar(ef$REF) > 1 |
    vapply(strsplit(ef$ALT, ","), function(a) any(nchar(a) > 1), logical(1))

  # ---- monomorphic-site counts (added 2026-09-23, requested for the paper) --
  # A site is MONOMORPHIC here if, among the rows that actually enter this
  # comparison, it has at least one call and every call is the same genotype.
  # The reference set differs per file, deliberately:
  #   wild  -- only the accessions crosshap ASSIGNED to a haplotype group for this
  #            gene (`ind$Ind`), because unassigned accessions are dropped from the
  #            test and from the figure, so a site that varies only among them
  #            carries nothing. The set therefore differs from gene to gene.
  #   elite -- only the 5 CONFIGURED lines, not the 136-line pool, because those
  #            five are the only rows drawn.
  # Counted over all sites of each file, over the sites kept in `shared_sites`,
  # and over the sites removed from it, so the three always add up.
  mono_frac <- function(mat, keys) {
    if (!length(keys)) return(0L)
    sub <- mat[, intersect(keys, colnames(mat)), drop = FALSE]
    sum(apply(sub, 2, function(x) { x <- x[!is.na(x)]
                                    length(x) > 0L && length(unique(x)) == 1L }))
  }
  W_assigned <- W$mat[intersect(ind$Ind, rownames(W$mat)), , drop = FALSE]
  wild_removed  <- setdiff(colnames(W$mat), shared_keys)   # wild-only + triallelic
  elite_removed <- setdiff(colnames(E$mat), shared_keys)   # elite-only + triallelic
  msg("  monomorphic | wild (", nrow(W_assigned), " assigned accessions): ",
      mono_frac(W_assigned, colnames(W$mat)), " of ", ncol(W$mat),
      " | elite (", nrow(E$mat), " lines): ", mono_frac(E$mat, colnames(E$mat)),
      " of ", ncol(E$mat))

  overlap_tbl[[gid]] <- data.frame(
    gene_id = gid, short_name = short, trait = trait, chr = gwr$chr,
    win_start = gwr$win_start, win_end = gwr$win_end,
    n_wild_snps = ncol(W$mat), n_elite_records = ncol(E$mat),
    n_shared = length(shared_keys), n_wild_only = length(wild_only),
    n_elite_only = length(elite_only),
    n_elite_only_indel = sum(e_is_indel[match(elite_only, E$key)], na.rm = TRUE),
    n_shared_positions = nrow(ac), n_alleles_identical = n_ident,
    n_alleles_swapped = n_swap, n_alleles_triallelic = n_altd,
    n_alleles_ref_differs = n_refd,
    # monomorphic counts -- see mono_frac() above for the two reference sets
    n_wild_assigned_accessions       = nrow(W_assigned),
    n_wild_mono_all                  = mono_frac(W_assigned, colnames(W$mat)),
    n_wild_mono_kept                 = mono_frac(W_assigned, shared_keys),
    n_wild_mono_removed              = mono_frac(W_assigned, wild_removed),
    n_elite_lines_shown              = nrow(E$mat),
    n_elite_mono_all                 = mono_frac(E$mat, colnames(E$mat)),
    n_elite_mono_kept                = mono_frac(E$mat, shared_keys),
    n_elite_mono_removed             = mono_frac(E$mat, elite_removed),
    stringsAsFactors = FALSE)
  msg("  wild ", ncol(W$mat), " SNPs | elite ", ncol(E$mat), " records | shared ",
      length(shared_keys), " | wild-only ", length(wild_only), " | elite-only ",
      length(elite_only))

  # ---- wild group consensus -----------------------------------------------
  groups <- sort(unique(ind$hap))
  cons <- t(vapply(groups, function(g) {
    ids <- intersect(ind$Ind[ind$hap == g], rownames(W$mat))
    consensus_row(W$mat[ids, , drop = FALSE])
  }, numeric(ncol(W$mat))))
  rownames(cons) <- groups; colnames(cons) <- colnames(W$mat)

  # group n and mean phenotype (for the violin labels)
  gstat <- ind %>% group_by(hap) %>%
    summarise(n = dplyr::n(), mean_pheno = mean(Pheno, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()

  # closest real accession to each consensus -- verification only, not published
  for (g in groups) {
    ids <- intersect(ind$Ind[ind$hap == g], rownames(W$mat))
    d <- vapply(ids, function(s) frac_disagree(W$mat[s, ], cons[g, ]), numeric(1))
    best <- names(which.min(d))
    rep_tbl[[paste(gid, g)]] <- data.frame(
      gene_id = gid, short_name = short, hap = g, n_members = length(ids),
      closest_accession = best,
      pct_agreement_with_consensus = round(100 * (1 - min(d, na.rm = TRUE)), 2),
      n_consensus_na = sum(is.na(cons[g, ])),
      stringsAsFactors = FALSE)
  }

  # ---- elite genotypes, long form (for the paper) --------------------------
  el <- as.data.frame(as.table(E$mat[, intersect(colnames(E$mat), colnames(W$mat)), drop = FALSE]),
                      stringsAsFactors = FALSE)
  if (nrow(el)) {
    names(el) <- c("line_name", "site_key", "code")
    el$gene_id <- gid; el$short_name <- short
    el$genotype <- c("REF", "ALT", "HET")[el$code + 1]
    el$genotype[is.na(el$code)] <- "missing"
    elite_long[[gid]] <- el[, c("gene_id", "short_name", "line_name", "site_key", "genotype")]
  }

  # ---- assemble the three versions ----------------------------------------
  for (ver in P$VERSIONS) {
    if (ver == "shared_sites") {
      cols <- shared_keys
      cons_v <- cons[, cols, drop = FALSE]
      el_v   <- E$mat[, cols, drop = FALSE]
    } else {
      cols <- colnames(W$mat)                      # every wild SNP
      cons_v <- cons[, cols, drop = FALSE]
      el_v   <- matrix(NA_real_, nrow = nrow(E$mat), ncol = length(cols),
                       dimnames = list(rownames(E$mat), cols))
      hit <- intersect(cols, colnames(E$mat))
      el_v[, hit] <- E$mat[, hit, drop = FALSE]
      # Triallelic columns: the elite record exists, so the line's REAL genotype
      # is drawn -- never a fill -- in both filled versions (added 2026-09-14).
      if (length(tri_wild_keys))
        el_v[, tri_wild_keys] <- tri_code[, tri_wild_keys, drop = FALSE]
      absent <- wild_only                          # no elite record at the position
      if (length(absent))
        el_v[, absent] <- if (ver == "filled_marked") CODE_NORECORD else CODE_REF
    }
    # order columns by genomic position
    ord <- order(as.integer(sub("^[^:]+:([0-9]+):.*$", "\\1", cols)))
    cons_v <- cons_v[, ord, drop = FALSE]; el_v <- el_v[, ord, drop = FALSE]

    saveRDS(list(gene_id = gid, short_name = short, trait = trait, version = ver,
                 chr = gwr$chr, gene_start = gwr$gene_start, gene_end = gwr$gene_end,
                 win_start = gwr$win_start, win_end = gwr$win_end,
                 consensus = cons_v, elite = el_v, group_stats = gstat,
                 indfile = ind, n_shared = length(shared_keys),
                 n_wild_only = length(wild_only), n_elite_only = length(elite_only),
                 n_triallelic = length(tri_wild_keys)),
            file.path(P$DIR_MATRICES, paste0(gid, "__", ver, ".rds")))
    msg("  [", ver, "] ", nrow(cons_v), " group rows x ", ncol(cons_v), " sites")
  }
}

write_tsv <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

write_tsv(do.call(rbind, allele_tbl),  file.path(P$DIR_TABLES, "allele_concordance.tsv"))
write_tsv(do.call(rbind, overlap_tbl), file.path(P$DIR_TABLES, "site_overlap_summary.tsv"))
write_tsv(do.call(rbind, rep_tbl),     file.path(P$DIR_TABLES, "consensus_representatives.tsv"))
write_tsv(do.call(rbind, elite_long),  file.path(P$DIR_TABLES, "elite_genotypes_long.tsv"))
write_tsv(if (length(tri_tbl)) do.call(rbind, tri_tbl) else
            data.frame(gene_id = character(), site = character(), line_name = character(),
                       genotype = character()),
          file.path(P$DIR_TABLES, "triallelic_sites_elite_genotypes.tsv"))

msg("DONE -> ", P$DIR_MATRICES)
