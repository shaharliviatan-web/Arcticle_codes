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

  if (n <= 1) {
    return(as.dist(d))
  }

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

cluster_block <- function(m01_block) {
  if (is.null(m01_block) || !is.matrix(m01_block) || nrow(m01_block) == 0) {
    return(list(mat = NULL, hc = NULL))
  }

  if (nrow(m01_block) == 1) {
    return(list(mat = m01_block, hc = NULL))
  }

  d <- pairwise_hamming_ignore_na(m01_block)
  hc <- hclust(d, method = "average")

  list(
    mat = m01_block[hc$order, , drop = FALSE],
    hc = hc
  )
}

make_display_names <- function(sample_ids) {
  disp_names <- sub("_.*", "", sample_ids)
  blank_idx <- grepl("^__separator__", sample_ids)

  if (any(blank_idx)) {
    disp_names[blank_idx] <- " "
  }

  non_blank_idx <- which(!blank_idx)
  if (length(non_blank_idx) > 0) {
    disp_names[non_blank_idx] <- make.unique(disp_names[non_blank_idx])
  }

  disp_names
}

shift_hclust_merge <- function(merge_mat, leaf_offset, node_offset) {
  out <- merge_mat

  if (length(out) == 0) {
    return(out)
  }

  out[out < 0] <- out[out < 0] - leaf_offset
  out[out > 0] <- out[out > 0] + node_offset
  out
}

build_combined_row_hclust <- function(labels, n_wild, n_separator_rows, n_elite, hc_wild, hc_elite) {
  n_total <- length(labels)
  if (n_total <= 1) {
    return(FALSE)
  }

  merge_list <- list()
  height_vec <- numeric(0)
  node_count <- 0L

  max_height_wild <- 0
  max_height_elite <- 0

  if (n_wild > 1 && !is.null(hc_wild)) {
    merge_list[[length(merge_list) + 1]] <- hc_wild$merge
    height_vec <- c(height_vec, hc_wild$height)
    node_count <- node_count + nrow(hc_wild$merge)
    max_height_wild <- max(hc_wild$height)
    root_wild <- node_count
  } else if (n_wild == 1) {
    root_wild <- -1L
  } else {
    root_wild <- NULL
  }

  elite_leaf_offset <- n_wild + n_separator_rows

  if (n_elite > 1 && !is.null(hc_elite)) {
    elite_merge_shifted <- shift_hclust_merge(
      hc_elite$merge,
      leaf_offset = elite_leaf_offset,
      node_offset = node_count
    )
    merge_list[[length(merge_list) + 1]] <- elite_merge_shifted
    height_vec <- c(height_vec, hc_elite$height)
    node_count <- node_count + nrow(hc_elite$merge)
    max_height_elite <- max(hc_elite$height)
    root_elite <- node_count
  } else if (n_elite == 1) {
    root_elite <- -(elite_leaf_offset + 1L)
  } else {
    root_elite <- NULL
  }

  next_height <- max(max_height_wild, max_height_elite, 0) + 0.5
  current_root <- NULL

  if (!is.null(root_wild)) {
    current_root <- root_wild
  }

  if (n_separator_rows > 0) {
    separator_indices <- seq.int(from = n_wild + 1L, length.out = n_separator_rows)

    for (separator_idx in separator_indices) {
      sep_ref <- -separator_idx
      if (is.null(current_root)) {
        current_root <- sep_ref
      } else {
        merge_list[[length(merge_list) + 1]] <- matrix(c(current_root, sep_ref), nrow = 1)
        height_vec <- c(height_vec, next_height)
        node_count <- node_count + 1L
        current_root <- node_count
        next_height <- next_height + 0.5
      }
    }
  }

  if (!is.null(root_elite)) {
    if (is.null(current_root)) {
      current_root <- root_elite
    } else {
      merge_list[[length(merge_list) + 1]] <- matrix(c(current_root, root_elite), nrow = 1)
      height_vec <- c(height_vec, next_height)
      node_count <- node_count + 1L
      current_root <- node_count
    }
  }

  merge_mat <- do.call(rbind, merge_list)
  if (is.null(merge_mat)) {
    return(FALSE)
  }

  structure(
    list(
      merge = merge_mat,
      height = height_vec,
      order = seq_len(n_total),
      labels = labels,
      method = "average",
      call = match.call(),
      dist.method = "custom"
    ),
    class = "hclust"
  )
}

build_group_heatmap_grob <- function(
  g,
  n_total,
  n_wild,
  n_elite,
  samples_g,
  m01_all,
  ann_row,
  annot_col,
  ann_colors,
  hm_colors,
  hm_breaks,
  cell_h_pt,
  target
) {
  m01_g <- m01_all[samples_g, , drop = FALSE]
  ann_g <- ann_row[samples_g, , drop = FALSE]

  wild_ids <- rownames(ann_g)[ann_g$Population == "Wild"]
  elite_ids <- rownames(ann_g)[ann_g$Population == "Elite"]

  wild_obj <- cluster_block(if (length(wild_ids) > 0) m01_g[wild_ids, , drop = FALSE] else NULL)
  elite_obj <- cluster_block(if (length(elite_ids) > 0) m01_g[elite_ids, , drop = FALSE] else NULL)

  m01_wild <- wild_obj$mat
  m01_elite <- elite_obj$mat
  hc_wild <- wild_obj$hc
  hc_elite <- elite_obj$hc

  row_blocks <- list()
  ann_blocks <- list()
  n_separator_rows <- 0L

  if (!is.null(m01_wild) && nrow(m01_wild) > 0) {
    row_blocks[[length(row_blocks) + 1]] <- m01_wild
    ann_blocks[[length(ann_blocks) + 1]] <- ann_row[rownames(m01_wild), , drop = FALSE]
  }

  if (!is.null(m01_wild) && nrow(m01_wild) > 0 &&
      !is.null(m01_elite) && nrow(m01_elite) > 0) {
    separator_rows <- matrix(
      2,
      nrow = 2,
      ncol = ncol(m01_g),
      dimnames = list(c("__separator1__", "__separator2__"), colnames(m01_g))
    )
    separator_ann <- data.frame(
      Population = c(" ", " "),
      row.names = c("__separator1__", "__separator2__"),
      stringsAsFactors = FALSE
    )

    row_blocks[[length(row_blocks) + 1]] <- separator_rows
    ann_blocks[[length(ann_blocks) + 1]] <- separator_ann
    n_separator_rows <- 2L
  }

  if (!is.null(m01_elite) && nrow(m01_elite) > 0) {
    row_blocks[[length(row_blocks) + 1]] <- m01_elite
    ann_blocks[[length(ann_blocks) + 1]] <- ann_row[rownames(m01_elite), , drop = FALSE]
  }

  m01_plot <- do.call(rbind, row_blocks)
  ann_plot <- do.call(rbind, ann_blocks)

  display_names <- make_display_names(rownames(m01_plot))
  rownames(m01_plot) <- display_names
  rownames(ann_plot) <- display_names

  row_hc <- build_combined_row_hclust(
    labels = rownames(m01_plot),
    n_wild = if (is.null(m01_wild)) 0L else nrow(m01_wild),
    n_separator_rows = n_separator_rows,
    n_elite = if (is.null(m01_elite)) 0L else nrow(m01_elite),
    hc_wild = hc_wild,
    hc_elite = hc_elite
  )

  main_title <- paste0(
    make_gene_name(target$gene_file), " | ", target$title, "\n",
    "epsilon=", target$epsilon, " | Haplotype group: ", g,
    " | n_total=", n_total, " (wild=", n_wild, ", elite=", n_elite, ")\n",
    "Top: clustered wild rows | Bottom: clustered elite rows | White rows = separator"
  )

  ph <- pheatmap::pheatmap(
    m01_plot,
    color = hm_colors,
    breaks = hm_breaks,
    na_col = "grey85",
    cluster_rows = row_hc,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    border_color = "white",
    annotation_row = ann_plot,
    annotation_col = annot_col,
    annotation_colors = ann_colors,
    main = main_title,
    cellheight = cell_h_pt,
    treeheight_row = 50,
    legend = FALSE,
    silent = TRUE
  )

  list(
    gtable = ph$gtable,
    height_weight = max(6, nrow(m01_plot) + 5)
  )
}

write_heatmaps_for_target <- function(HapObject, target, merged_path, common_pos, vcf_merged_ids, out_pdf) {
  target <- normalize_target(target)
  label <- paste0("Haplotypes_MGmin", target$MGmin, "_E", target$epsilon)
  if (is.null(HapObject[[label]])) return(invisible(NULL))

  ensure_dir(dirname(out_pdf))

  hap_map <- HapObject[[label]]$Indfile %>%
    select(Ind, hap) %>%
    distinct() %>%
    mutate(hap = as.character(hap))

  var_map <- HapObject[[label]]$Varfile %>%
    select(ID, MGs, meanr2) %>%
    mutate(MarkerGroup = ifelse(MGs == "0", "Unlinked", as.character(MGs)))

  v <- vcfR::read.vcfR(merged_path, verbose = FALSE)
  gt <- vcfR::extract.gt(v, element = "GT", as.numeric = FALSE)
  fix <- as.data.frame(vcfR::getFIX(v), stringsAsFactors = FALSE)
  var_base_all <- paste0(fix$CHROM, ":", fix$POS)

  keep_idx <- var_base_all %in% common_pos
  gt <- gt[keep_idx, , drop = FALSE]
  var_base_filt <- var_base_all[keep_idx]
  var_ids <- make.unique(var_base_filt)

  if (length(var_ids) != length(vcf_merged_ids) || any(var_ids != vcf_merged_ids)) {
    stop("Variant ID mismatch between heatmap matrix and CrossHap IDs")
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
    stop("No matching sample IDs between VCF and CrossHap Indfile")
  }

  m01_all <- m01_all[in_both, , drop = FALSE]
  hap_map <- hap_map %>% filter(Ind %in% in_both)

  ann_row <- data.frame(
    Population = ifelse(grepl("^HS", rownames(m01_all)), "Wild", "Elite"),
    row.names = rownames(m01_all),
    stringsAsFactors = FALSE
  )

  annot_col <- data.frame(row.names = colnames(m01_all), stringsAsFactors = FALSE)
  annot_col$MarkerGroup <- "Unlinked"
  matched <- intersect(colnames(m01_all), var_map$ID)
  if (length(matched) > 0) {
    annot_col$MarkerGroup[match(matched, rownames(annot_col))] <- var_map$MarkerGroup[match(matched, var_map$ID)]
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
        ", R2=", format(r2_map[[x]], nsmall = 2), ")"
      )
    }
  }, character(1))

  unique_mgs <- unique(annot_col$MarkerGroup)
  mg_colors <- setNames(mk_palette(length(unique_mgs)), unique_mgs)
  unlinked_name <- unique_mgs[grepl("^Unlinked", unique_mgs)]
  if (length(unlinked_name) > 0) {
    mg_colors[unlinked_name] <- "#E0E0E0"
  }

  ann_colors <- list(
    Population = c(Wild = "#1b9e77", Elite = "#d95f02", " " = "#FFFFFF"),
    MarkerGroup = mg_colors
  )

  group_stats <- hap_map %>%
    mutate(
      Population = ifelse(grepl("^HS", Ind), "Wild", "Elite"),
      hap_num = suppressWarnings(as.numeric(hap)),
      hap_is_zero = !is.na(hap_num) & hap_num == 0
    ) %>%
    group_by(hap, hap_num, hap_is_zero) %>%
    summarize(
      n_total = n(),
      n_wild = sum(Population == "Wild"),
      n_elite = sum(Population == "Elite"),
      .groups = "drop"
    ) %>%
    arrange(hap_is_zero, desc(n_total), hap_num, hap)

  hm_colors <- c("#FFFACD", "#2F4F4F", "#FFFFFF")
  hm_breaks <- c(-0.5, 0.5, 1.5, 2.5)
  cell_h_pt <- 10

  group_grobs <- list()
  group_heights <- c()

  for (k in seq_len(nrow(group_stats))) {
    g <- group_stats$hap[k]
    n_total <- group_stats$n_total[k]
    n_wild <- group_stats$n_wild[k]
    n_elite <- group_stats$n_elite[k]

    samples_g <- hap_map$Ind[hap_map$hap == g]
    samples_g <- intersect(samples_g, rownames(m01_all))
    if (length(samples_g) == 0) {
      next
    }

    group_obj <- build_group_heatmap_grob(
      g = g,
      n_total = n_total,
      n_wild = n_wild,
      n_elite = n_elite,
      samples_g = samples_g,
      m01_all = m01_all,
      ann_row = ann_row,
      annot_col = annot_col,
      ann_colors = ann_colors,
      hm_colors = hm_colors,
      hm_breaks = hm_breaks,
      cell_h_pt = cell_h_pt,
      target = target
    )

    group_grobs[[length(group_grobs) + 1]] <- group_obj$gtable
    group_heights <- c(group_heights, group_obj$height_weight)
  }

  if (length(group_grobs) == 0) {
    return(invisible(NULL))
  }

  between_group_space_in <- 0.45

  total_height_units <- sum(group_heights)
  pdf_h <- (total_height_units * cell_h_pt) / 72 + 2 + (length(group_grobs) - 1) * between_group_space_in
  pdf_h <- min(120, max(12, pdf_h))

  grDevices::pdf(out_pdf, width = 24, height = pdf_h, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  layout_heights <- vector("list", 2 * length(group_grobs) - 1)
  panel_rows <- integer(length(group_grobs))

  for (i in seq_along(group_grobs)) {
    layout_heights[[2 * i - 1]] <- grid::unit(group_heights[i], "null")
    panel_rows[i] <- 2 * i - 1

    if (i < length(group_grobs)) {
      layout_heights[[2 * i]] <- grid::unit(between_group_space_in, "in")
    }
  }

  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = length(layout_heights),
        ncol = 1,
        heights = do.call(grid::unit.c, layout_heights)
      )
    )
  )

  for (i in seq_along(group_grobs)) {
    grid::pushViewport(grid::viewport(layout.pos.row = panel_rows[i], layout.pos.col = 1))
    grid::grid.draw(group_grobs[[i]])
    grid::upViewport()
  }

  grid::upViewport()

  invisible(out_pdf)
}