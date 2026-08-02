# figure_utils.R
# Helpers for the publication haplotype figures (violin + per-group genotype heatmap).
# Wild-only (290 accessions); NO elite/cultivated lines anywhere in these figures.
# Ported and cleaned from the legacy Wild_And_Elite presentation tool; elite rows,
# elite-guided representative selection, and the "npg" palette have been removed.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

# ---- Okabe-Ito colorblind-safe palette (pale yellow demoted; black last) -----
OKABE_ITO <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00",
               "#56B4E9", "#F0E442", "#000000")

# Heatmap tile colors (kept from the legacy figures; good REF/ALT contrast).
HEAT_COLS <- c(Reference = "#d8cca0", Alternate = "#2F4F4F", Missing = "grey85")

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path)
}

# Order haplotype letters A,B,C... (numeric-aware fallback).
order_haps <- function(haps) {
  haps <- unique(as.character(haps))
  hnum <- suppressWarnings(as.numeric(haps))
  if (all(!is.na(hnum))) haps[order(hnum)] else sort(haps)
}

p_to_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p <= 1e-4) return("****")
  if (p <= 1e-3) return("***")
  if (p <= 1e-2) return("**")
  if (p <= 5e-2) return("*")
  "ns"
}

# ---------------------------------------------------------------------------
# Robust VCF reader: reproduces crosshap::read_vcf's transforms exactly so the
# reconstructed variant set/order matches the cached haplotyping IDs. Strips the
# "##" meta lines first (so fread never mis-detects a comma separator), collapses
# each genotype to its first field, "/"->"|", numeric POS, "-"->"." in colnames.
# ---------------------------------------------------------------------------
read_vcf_robust <- function(path) {
  cmd <- paste("zcat -f", shQuote(path), "| grep -v '^##'")
  vcf <- data.table::fread(cmd = cmd, header = TRUE, sep = "\t", nThread = 4) %>%
    as.data.frame(stringsAsFactors = FALSE)
  gt_cols <- setdiff(colnames(vcf),
                     c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))
  for (cc in colnames(vcf)) {
    vcf[[cc]] <- gsub(":.*", "", gsub("/", "|", as.character(vcf[[cc]])))
  }
  vcf$POS <- as.numeric(vcf$POS)
  colnames(vcf) <- gsub("-", ".", colnames(vcf))
  vcf
}

gt_to_bin01 <- function(g) {
  out <- rep(NA_real_, length(g))
  out[g %in% c("0|0", "0/0")] <- 0
  out[g %in% c("1|1", "1/1")] <- 1
  out  # het and missing -> NA (wild accessions are near-fully homozygous)
}

# Reconstruct the samples x variants 0/1 matrix using EXACTLY the variants and
# order crosshap haplotyped on, asserting identity against the cached IDs.
build_geno_matrix <- function(raw_path, common_ids, expected_ids, MGmin) {
  meta <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  vcf <- read_vcf_robust(raw_path)
  ord <- order(vcf[["#CHROM"]], vcf[["POS"]], vcf[["REF"]], vcf[["ALT"]])
  vcf <- vcf[ord, , drop = FALSE]
  id_base <- paste0(vcf[["#CHROM"]], ":", vcf[["POS"]])
  vcf <- vcf[id_base %in% common_ids, , drop = FALSE]
  id_base <- paste0(vcf[["#CHROM"]], ":", vcf[["POS"]])
  ids <- make.unique(id_base)
  if (!identical(ids, expected_ids)) {
    stop("Reconstructed variant IDs do not match cached vcf_raw_ids (raw VCF / cache mismatch).")
  }
  sample_cols <- setdiff(colnames(vcf), meta)
  gt <- vcf[, sample_cols, drop = FALSE]
  m <- vapply(gt, gt_to_bin01, numeric(nrow(gt)))      # variants x samples
  m <- matrix(as.numeric(m), nrow = nrow(gt), ncol = length(sample_cols),
              dimnames = list(ids, sample_cols))
  t(m)                                                  # samples x variants
}

# Representative wild accession for a haplotype group: fewest missing genotypes,
# then alphabetical (deterministic). No elite guidance (wild-only figures).
pick_representative <- function(m_group) {
  miss <- rowSums(is.na(m_group))
  zero <- names(miss)[miss == 0]
  if (length(zero) > 0) return(sort(zero)[1])
  names(sort(miss))[1]
}

# ---------------------------------------------------------------------------
# Violin panel: phenotype distribution per haplotype group (hap 0 excluded),
# boxplot inside, n labels, Kruskal-Wallis P subtitle, pairwise Wilcoxon (Holm)
# brackets vs the LARGEST group, drawn as significance stars. Returns the plot
# plus all numbers needed for the results-chapter CSV.
# ---------------------------------------------------------------------------
build_violin <- function(indfile, title, y_label, base_size = 20) {
  d <- as.data.frame(indfile, stringsAsFactors = FALSE)
  d$hap   <- as.character(d$hap)
  d$Pheno <- suppressWarnings(as.numeric(d$Pheno))
  d <- d[d$hap != "0" & !is.na(d$Pheno), , drop = FALSE]
  if (nrow(d) == 0) stop("No grouped, phenotyped rows for violin: ", title)

  hap_levels <- order_haps(d$hap)
  n_tab <- d %>% count(hap, name = "n") %>%
    mutate(hap = as.character(hap)) %>% arrange(match(hap, hap_levels))

  # Reference = largest haplotype group (ties -> first by letter).
  ref_hap <- n_tab$hap[order(-n_tab$n, match(n_tab$hap, hap_levels))][1]

  lab_map <- setNames(paste0(n_tab$hap, "\n(n=", n_tab$n, ")"), n_tab$hap)
  d$hap_label <- factor(lab_map[d$hap], levels = unname(lab_map[hap_levels]))

  kw <- kruskal.test(Pheno ~ hap, data = d)
  k  <- length(hap_levels); n_tot <- nrow(d)
  eps2 <- as.numeric((kw$statistic - k + 1) / (n_tot - k))   # epsilon-squared effect size

  # group descriptive stats
  gstats <- d %>% group_by(hap) %>%
    summarise(n = dplyr::n(), mean = mean(Pheno), sd = sd(Pheno), median = median(Pheno),
              q25 = quantile(Pheno, .25), q75 = quantile(Pheno, .75),
              min = min(Pheno), max = max(Pheno), .groups = "drop") %>%
    arrange(match(hap, hap_levels))

  # pairwise Wilcoxon vs reference, Holm-adjusted
  others <- hap_levels[hap_levels != ref_hap]
  pw <- data.frame()
  if (length(others) > 0) {
    pw <- do.call(rbind, lapply(others, function(h) {
      x <- d$Pheno[d$hap == ref_hap]; y <- d$Pheno[d$hap == h]
      praw <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
      data.frame(hap = h, p_raw = praw, stringsAsFactors = FALSE)
    }))
    pw$p_holm   <- p.adjust(pw$p_raw, method = "holm")
    pw$p.signif <- vapply(pw$p_holm, p_to_stars, character(1))
    pw$group1   <- lab_map[[ref_hap]]
    pw$group2   <- lab_map[pw$hap]
    yt <- max(d$Pheno, na.rm = TRUE); yr <- diff(range(d$Pheno, na.rm = TRUE))
    if (!is.finite(yr) || yr == 0) yr <- 1
    pw$y.position <- yt + seq_len(nrow(pw)) * (0.08 * yr)
  }

  pal <- setNames(OKABE_ITO[seq_along(hap_levels)], unname(lab_map[hap_levels]))

  p <- ggviolin(d, x = "hap_label", y = "Pheno", fill = "hap_label",
                add = "boxplot", add.params = list(fill = "white", width = 0.12)) +
    scale_fill_manual(values = pal) +
    labs(title = title,
         subtitle = paste0("Kruskal-Wallis P = ", formatC(kw$p.value, format = "e", digits = 2)),
         x = "Haplotype group", y = y_label) +
    theme_bw(base_size = base_size) +
    theme(legend.position = "none",
          plot.title = element_text(size = 27, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 19, face = "bold", hjust = 0),
          axis.title = element_text(size = 23, face = "bold"),
          axis.text  = element_text(size = 20, face = "bold"),
          panel.grid.major.x = element_blank())

  if (nrow(pw) > 0) {
    p <- p + stat_pvalue_manual(pw, label = "p.signif", xmin = "group1", xmax = "group2",
                                y.position = "y.position", tip.length = 0.01,
                                bracket.size = 0.5, size = 9, fontface = "bold")
  }

  list(plot = p, ref_hap = ref_hap, hap_levels = hap_levels,
       kw_stat = as.numeric(kw$statistic), kw_df = as.integer(kw$parameter),
       kw_p = as.numeric(kw$p.value), eps2 = eps2, n_total = n_tot, n_groups = k,
       gstats = gstats, pairwise = pw)
}

# ---------------------------------------------------------------------------
# Heatmap panel: one representative-genotype row per haplotype group (wild-only;
# no cultivar rows). REF=tan, ALT=slate, Missing=grey. SNP positions not labeled.
# ---------------------------------------------------------------------------
build_heatmap <- function(m01, hap_map, hap_levels, base_size = 18, tile_border = 0.9) {
  rep_info <- list(); rows <- list()
  for (h in hap_levels) {
    ids <- intersect(hap_map$Ind[hap_map$hap == h], rownames(m01))
    if (length(ids) == 0) next
    rep_id <- pick_representative(m01[ids, , drop = FALSE])
    vals <- as.numeric(m01[rep_id, ])
    rep_info[[length(rep_info) + 1]] <- data.frame(
      hap = h, representative_accession = rep_id, n_snps = length(vals),
      n_alt_in_rep = sum(vals == 1, na.rm = TRUE),
      n_ref_in_rep = sum(vals == 0, na.rm = TRUE),
      n_missing_in_rep = sum(is.na(vals)), stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      row_label = paste0("Haplotype group ", h),
      variant_index = seq_len(ncol(m01)), value = vals, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) stop("No haplotype-group representatives for heatmap.")
  rep_info <- do.call(rbind, rep_info)
  long <- do.call(rbind, rows)
  long$fill_key <- factor(ifelse(is.na(long$value), "Missing",
                           ifelse(long$value == 0, "Reference", "Alternate")),
                          levels = c("Reference", "Alternate", "Missing"))
  row_order <- paste0("Haplotype group ", hap_levels)
  row_order <- row_order[row_order %in% long$row_label]
  long$row_label <- factor(long$row_label, levels = rev(row_order))

  p <- ggplot(long, aes(variant_index, row_label, fill = fill_key)) +
    geom_tile(color = "white", linewidth = tile_border, height = 0.86) +
    scale_fill_manual(values = HEAT_COLS, breaks = c("Reference","Alternate","Missing"),
                      labels = c("REF","ALT","Missing"), drop = FALSE, name = NULL) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(axis.text.y = element_text(size = 22, face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 18, face = "bold"),
          legend.key.size = grid::unit(1.1, "lines"),
          legend.position = "right")

  list(plot = p, rep_info = rep_info)
}
