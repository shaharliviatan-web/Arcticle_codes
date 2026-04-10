suppressPackageStartupMessages({
  library(dplyr)
  library(vcfR)
  library(pheatmap)
  library(grid)
  library(RColorBrewer)
})

gt_to_bin01 <- function(gt_chr) {
  out <- rep(NA_real_, length(gt_chr))
  out[gt_chr %in% c("0/0", "0|0")] <- 0
  out[gt_chr %in% c("1/1", "1|1")] <- 1
  out[gt_chr %in% c("./.", ".|.", ".")] <- NA_real_
  out
}

pairwise_hamming_ignore_na <- function(m01) {
  stopifnot(is.matrix(m01))
  n <- nrow(m01)
  d <- matrix(0, n, n, dimnames = list(rownames(m01), rownames(m01)))

  for (i in seq_len(n - 1)) {
    xi <- m01[i, ]
    for (j in (i + 1):n) {
      xj <- m01[j, ]
      ok <- !is.na(xi) & !is.na(xj)
      n_ok <- sum(ok)
      d[i, j] <- d[j, i] <- if (n_ok == 0) 1 else sum(xi[ok] != xj[ok]) / n_ok
    }
  }

  as.dist(d)
}

mk_palette <- function(n) {
  if (n <= 9) {
    brewer.pal(max(3, n), "Set1")[seq_len(n)]
  } else {
    grDevices::colorRampPalette(brewer.pal(9, "Set1"))(n)
  }
}

write_heatmaps_for_eps <- function(
  HapObject,
  label,
  gene_file,
  title,
  trait,
  eps,
  MGmin,
  raw_path,
  common_ids,
  vcf_raw_ids,
  out_pdf
) {
  if (is.null(HapObject[[label]])) {
    return(invisible(NULL))
  }

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

  hap_map <- HapObject[[label]]$Indfile %>%
    select(Ind, hap) %>%
    distinct() %>%
    mutate(hap = as.character(hap))

  var_map <- HapObject[[label]]$Varfile %>%
    select(ID, MGs, meanr2) %>%
    mutate(MarkerGroup = ifelse(MGs == "0", "Unlinked", as.character(MGs)))

  v <- vcfR::read.vcfR(raw_path, verbose = FALSE)
  gt <- vcfR::extract.gt(v, element = "GT", as.numeric = FALSE)

  fix <- as.data.frame(vcfR::getFIX(v), stringsAsFactors = FALSE)
  var_base_all <- paste0(fix$CHROM, ":", fix$POS)

  keep_idx <- var_base_all %in% common_ids
  if (sum(keep_idx) < MGmin) {
    stop("Too few variants after common-ID filtering in vcfR matrix.")
  }

  gt <- gt[keep_idx, , drop = FALSE]
  fix <- fix[keep_idx, , drop = FALSE]
  var_base_filt <- var_base_all[keep_idx]
  var_ids <- make.unique(var_base_filt)

  if (length(var_ids) != length(vcf_raw_ids) || any(var_ids != vcf_raw_ids)) {
    stop(
      "Variant ID mismatch between heatmap matrix and CrossHap IDs.\n",
      "First heatmap IDs: ", paste(head(var_ids, 5), collapse = ", "),
      "\nFirst crosshap IDs: ", paste(head(vcf_raw_ids, 5), collapse = ", ")
    )
  }

  gt_bin <- apply(gt, 2, gt_to_bin01)
  gt_bin <- matrix(
    as.numeric(gt_bin),
    nrow = nrow(gt),
    ncol = ncol(gt),
    dimnames = list(var_ids, colnames(gt))
  )

  m01_all <- t(gt_bin)

  in_both <- intersect(rownames(m01_all), hap_map$Ind)
  if (length(in_both) == 0) {
    stop("No matching sample IDs between VCF and CrossHap Indfile.")
  }

  m01_all <- m01_all[in_both, , drop = FALSE]
  hap_map <- hap_map %>% filter(Ind %in% in_both)

  annot_col <- data.frame(row.names = colnames(m01_all))
  annot_col$MarkerGroup <- "Unlinked"

  matched <- intersect(colnames(m01_all), var_map$ID)
  if (length(matched) > 0) {
    annot_col$MarkerGroup[match(matched, rownames(annot_col))] <-
      var_map$MarkerGroup[match(matched, var_map$ID)]
  }

  mg_counts <- table(annot_col$MarkerGroup)

  r2_tbl <- var_map %>%
    filter(MarkerGroup != "Unlinked") %>%
    group_by(MarkerGroup) %>%
    summarize(avg_r2 = round(mean(meanr2, na.rm = TRUE), 2), .groups = "drop")

  r2_map <- setNames(r2_tbl$avg_r2, r2_tbl$MarkerGroup)

  annot_col$MarkerGroup <- vapply(annot_col$MarkerGroup, function(x) {
    if (x == "Unlinked") {
      n_unlinked <- if ("Unlinked" %in% names(mg_counts)) as.integer(mg_counts[["Unlinked"]]) else 0L
      paste0("Unlinked (n=", n_unlinked, ")")
    } else {
      paste0(
        x,
        " (n=", as.integer(mg_counts[[x]]),
        ", R2=", format(r2_map[[x]], nsmall = 2),
        ")"
      )
    }
  }, character(1))

  unique_mgs <- unique(annot_col$MarkerGroup)
  mg_colors <- setNames(mk_palette(length(unique_mgs)), unique_mgs)
  unlinked_name <- unique_mgs[grepl("^Unlinked", unique_mgs)]
  if (length(unlinked_name) > 0) {
    mg_colors[unlinked_name] <- "#E0E0E0"
  }

  ann_colors <- list(MarkerGroup = mg_colors)

  group_stats <- hap_map %>%
    group_by(hap) %>%
    summarize(n_total = dplyr::n(), .groups = "drop") %>%
    arrange(desc(n_total), hap)

  cell_h_pt <- 10
  extra_h_in <- 2.5
  max_group_n <- max(group_stats$n_total)
  pdf_h <- (max_group_n * cell_h_pt) / 72 + extra_h_in
  pdf_h <- min(80, max(8, pdf_h))

  gene_clean_name <- sub("\\.vcf\\.gz$", "", gene_file)

  hm_colors <- c("#FFFACD", "#2F4F4F")
  hm_breaks <- c(-0.5, 0.5, 1.5)

  grDevices::pdf(out_pdf, width = 24, height = pdf_h, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (k in seq_len(nrow(group_stats))) {
    g <- group_stats$hap[k]
    n_total <- group_stats$n_total[k]

    samples_g <- hap_map$Ind[hap_map$hap == g]
    samples_g <- intersect(samples_g, rownames(m01_all))
    if (length(samples_g) == 0) next

    m01_g <- m01_all[samples_g, , drop = FALSE]

    disp_names <- sub("_.*", "", rownames(m01_g))
    rownames(m01_g) <- make.unique(disp_names)

    if (nrow(m01_g) >= 2) {
      d <- pairwise_hamming_ignore_na(m01_g)
      hc <- hclust(d, method = "average")
      cluster_rows_arg <- hc
    } else {
      cluster_rows_arg <- FALSE
    }

    main_title <- paste0(
      gene_clean_name, " | ", title, "\n",
      "epsilon=", eps, " | Haplotype group: ", g, " | n_total=", n_total, "\n",
      "Genotypes: Light Yellow = Reference | Dark Teal = Alternate | Grey = Missing (NA)"
    )

    ph <- pheatmap::pheatmap(
      m01_g,
      color = hm_colors,
      breaks = hm_breaks,
      na_col = "grey70",
      cluster_rows = cluster_rows_arg,
      cluster_cols = FALSE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      annotation_col = annot_col,
      annotation_colors = ann_colors,
      fontsize_row = 8,
      border_color = "white",
      main = main_title,
      cellheight = cell_h_pt,
      legend = FALSE,
      silent = TRUE
    )

    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
  }

  invisible(out_pdf)
}