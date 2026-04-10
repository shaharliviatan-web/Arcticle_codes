suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(vcfR)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# ==============================================================================
# CONFIG
# ==============================================================================

RUN_DIR <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_Haplotype_Analysis/Wild_And_Elite/Target_Genes_analysis_Crosshap_Heatmap_Violin_Graphs/Runs/8_Genes_Analysis_Upgraded_Heatmaps_V2"

TRAIT <- "starch"
GENE_FILE <- "HORVU.MOREX.r3.3HG0301770.vcf.gz"
GENE_TITLE <- "Alpha-glucan phosphorylase alpha 1-4"

EPSILON <- 0.8
MGMIN <- 3

REFERENCE_HAP_ID <- "C"

# Violin stays like the current tool: hap 0 excluded.
EXCLUDE_HAP0_FROM_VIOLIN <- TRUE

# Heatmap presentation excludes haplotype 0.
INCLUDE_HAP0_IN_HEATMAP <- FALSE

OUTPUT_DIR <- file.path(RUN_DIR, "Presentation_PNG")
PNG_WIDTH_IN <- 18
PNG_HEIGHT_IN <- 16
PNG_DPI <- 300

VIOLIN_BASE_SIZE <- 18
HEATMAP_BASE_SIZE <- 16
HEATMAP_TILE_BORDER <- 0.9

# ==============================================================================
# HELPERS
# ==============================================================================

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

make_gene_name <- function(gene_file) {
  sub("\\.vcf\\.gz$", "", gene_file)
}

format_p_plain <- function(p, digits = 4, min_print = 1e-4) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
  if (p < min_print) return(paste0("<", formatC(min_print, format = "f", digits = digits)))
  formatC(p, format = "f", digits = digits)
}

p_to_signif_symbol <- function(p) {
  if (is.na(p)) return("ns")
  if (p <= 0.0001) return("****")
  if (p <= 0.001) return("***")
  if (p <= 0.01) return("**")
  if (p <= 0.05) return("*")
  "ns"
}

load_pheno_wild_elite <- function(pheno_file) {
  raw_pheno <- fread(pheno_file, header = FALSE, fill = TRUE)
  colnames(raw_pheno) <- c("FID", "IID", "trait")
  raw_pheno$Ind <- ifelse(grepl("^HS", raw_pheno$FID), paste(raw_pheno$FID, raw_pheno$IID, sep = "_"), raw_pheno$FID)
  raw_pheno$Pheno <- suppressWarnings(as.numeric(raw_pheno$trait))
  unique(raw_pheno[, .(Ind, Pheno)])
}

gt_to_bin01 <- function(gt_chr) {
  out <- rep(NA_real_, length(gt_chr))
  out[gt_chr %in% c("0/0", "0|0")] <- 0
  out[gt_chr %in% c("1/1", "1|1")] <- 1
  out[gt_chr %in% c("./.", ".|.", ".")] <- NA_real_
  out
}

order_haps <- function(haps) {
  haps <- unique(as.character(haps))
  hnum <- suppressWarnings(as.numeric(haps))
  if (all(!is.na(hnum))) {
    haps[order(hnum)]
  } else {
    sort(haps)
  }
}

short_sample_name <- function(x) {
  sub("_.*", "", x)
}

pick_representative <- function(m01_group) {
  miss_n <- rowSums(is.na(m01_group))
  zero_miss <- names(miss_n)[miss_n == 0]

  if (length(zero_miss) > 0) {
    return(sort(zero_miss)[1])
  }

  names(sort(miss_n, na.last = TRUE))[1]
}

mean_distance_to_elites <- function(wild_vec, elite_mat) {
  if (is.null(elite_mat) || !is.matrix(elite_mat) || nrow(elite_mat) == 0) {
    return(Inf)
  }

  dists <- apply(elite_mat, 1, function(elite_vec) {
    ok <- !is.na(wild_vec) & !is.na(elite_vec)
    n_ok <- sum(ok)
    if (n_ok == 0) {
      return(1)
    }
    sum(wild_vec[ok] != elite_vec[ok]) / n_ok
  })

  mean(dists)
}

pick_representative_guided <- function(wild_mat, elite_mat = NULL) {
  if (is.null(wild_mat) || !is.matrix(wild_mat) || nrow(wild_mat) == 0) {
    return(NULL)
  }

  if (is.null(elite_mat) || !is.matrix(elite_mat) || nrow(elite_mat) == 0) {
    return(pick_representative(wild_mat))
  }

  candidate_tbl <- data.frame(
    sample_id = rownames(wild_mat),
    mean_elite_distance = apply(wild_mat, 1, function(x) mean_distance_to_elites(x, elite_mat)),
    missing_n = rowSums(is.na(wild_mat)),
    stringsAsFactors = FALSE
  )

  zero_missing_tbl <- candidate_tbl[candidate_tbl$missing_n == 0, , drop = FALSE]
  if (nrow(zero_missing_tbl) > 0) {
    candidate_tbl <- zero_missing_tbl
  }

  candidate_tbl <- candidate_tbl[order(candidate_tbl$mean_elite_distance, candidate_tbl$missing_n, candidate_tbl$sample_id), , drop = FALSE]
  candidate_tbl$sample_id[1]
}

get_group_representative_info <- function(group_id, hap_map, m01_all) {
  hap_map_group <- hap_map %>% filter(hap == group_id)

  wild_ids <- hap_map_group %>%
    filter(Population == "Wild") %>%
    pull(Ind)

  elite_ids <- hap_map_group %>%
    filter(Population == "Elite") %>%
    pull(Ind)

  if (length(wild_ids) == 0) {
    return(NULL)
  }

  wild_mat <- m01_all[wild_ids, , drop = FALSE]
  elite_mat <- NULL
  if (length(elite_ids) > 0) {
    elite_ids <- sort(intersect(elite_ids, rownames(m01_all)))
    if (length(elite_ids) > 0) {
      elite_mat <- m01_all[elite_ids, , drop = FALSE]
    }
  }

  rep_id <- pick_representative_guided(wild_mat, elite_mat)
  if (is.null(rep_id) || !nzchar(rep_id)) {
    return(NULL)
  }

  list(
    group_id = group_id,
    rep_id = rep_id,
    elite_ids = elite_ids
  )
}

build_violin_plot <- function(HapObject, target, pheno, reference_hap_id, exclude_hap0 = TRUE) {
  label <- paste0("Haplotypes_MGmin", target$MGmin, "_E", target$epsilon)
  indfile <- HapObject[[label]]$Indfile

  if (is.null(indfile) || nrow(indfile) == 0) {
    stop("No Indfile data for ", label)
  }

  plot_df <- as.data.frame(indfile, stringsAsFactors = FALSE)
  if (!"Pheno" %in% names(plot_df)) {
    plot_df <- left_join(plot_df, pheno, by = "Ind")
  }

  plot_df$hap <- as.character(plot_df$hap)
  plot_df$Pheno <- suppressWarnings(as.numeric(plot_df$Pheno))

  if (exclude_hap0) {
    plot_df <- plot_df %>% filter(hap != "0")
  }

  wild_data <- plot_df %>% filter(!is.na(Pheno))
  elite_data <- plot_df %>% filter(is.na(Pheno))

  if (nrow(wild_data) == 0) {
    stop("No wild rows available for violin plot")
  }

  hap_levels <- order_haps(wild_data$hap)

  n_table <- wild_data %>%
    count(hap, name = "n_wild") %>%
    mutate(hap = as.character(hap)) %>%
    arrange(match(hap, hap_levels))

  hap_label_map <- setNames(paste0(n_table$hap, "\n(n=", n_table$n_wild, ")"), n_table$hap)

  wild_data$hap_label <- hap_label_map[wild_data$hap]
  elite_data$hap_label <- hap_label_map[elite_data$hap]

  label_levels <- unname(hap_label_map[hap_levels])
  wild_data$hap_label <- factor(wild_data$hap_label, levels = label_levels)
  elite_data$hap_label <- factor(elite_data$hap_label, levels = label_levels)

  if (!(reference_hap_id %in% unique(wild_data$hap))) {
    stop("REFERENCE_HAP_ID not found among wild haplotypes in violin data: ", reference_hap_id)
  }

  ref_group_label <- hap_label_map[[reference_hap_id]]

  kw_fit <- kruskal.test(Pheno ~ hap, data = wild_data)

  other_haps <- hap_levels[hap_levels != reference_hap_id]
  pairwise_tbl <- data.frame()

  if (length(other_haps) > 0) {
    pairwise_tbl <- do.call(rbind, lapply(other_haps, function(hap_id) {
      x <- wild_data$Pheno[wild_data$hap == reference_hap_id]
      y <- wild_data$Pheno[wild_data$hap == hap_id]
      p_raw <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)

      data.frame(
        group1 = hap_label_map[[reference_hap_id]],
        group2 = hap_label_map[[hap_id]],
        p = p_raw,
        stringsAsFactors = FALSE
      )
    }))

    pairwise_tbl$p.adj <- p.adjust(pairwise_tbl$p, method = "holm")
    pairwise_tbl$p.signif <- vapply(pairwise_tbl$p.adj, p_to_signif_symbol, character(1))

    y_top <- max(wild_data$Pheno, na.rm = TRUE)
    y_bottom <- min(wild_data$Pheno, na.rm = TRUE)
    y_range <- y_top - y_bottom
    if (!is.finite(y_range) || y_range == 0) y_range <- 1

    pairwise_tbl$y.position <- y_top + seq_len(nrow(pairwise_tbl)) * (0.08 * y_range)
  }

  y_range <- max(wild_data$Pheno, na.rm = TRUE) - min(wild_data$Pheno, na.rm = TRUE)
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  dummy_y <- min(wild_data$Pheno, na.rm = TRUE) - (y_range * 0.15)
  elite_data$DummyY <- dummy_y

  subtitle_text <- paste0(
    "Kruskal-Wallis P = ", formatC(kw_fit$p.value, format = "e", digits = 2)
  )

  p_violin <- ggviolin(
    wild_data,
    x = "hap_label",
    y = "Pheno",
    fill = "hap_label",
    palette = "npg",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
    geom_jitter(
      data = elite_data,
      aes(x = hap_label, y = DummyY),
      color = "black",
      size = 2,
      width = 0.1
    ) +
    geom_text_repel(
      data = elite_data,
      aes(x = hap_label, y = DummyY, label = Ind),
      size = 5.0,
      color = "#2ca25f",
      fontface = "bold",
      max.overlaps = Inf,
      direction = "y",
      nudge_y = -(y_range * 0.05)
    ) +
    labs(
      title = target$title,
      subtitle = subtitle_text,
      x = "Haplotype Group",
      y = "Phenotype Value"
    ) +
    theme_bw(base_size = VIOLIN_BASE_SIZE) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 17, face = "bold"),
      axis.title = element_text(size = 19, face = "bold"),
      axis.text.x = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      panel.grid.major.x = element_blank()
    )

  if (nrow(pairwise_tbl) > 0) {
    p_violin <- p_violin + stat_pvalue_manual(
      pairwise_tbl,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.45,
      size = 7,
      fontface = "bold"
    )
  }

  list(
    plot = p_violin,
    kw_p = kw_fit$p.value
  )
}

extract_genotype_matrix <- function(merged_path, common_pos, vcf_merged_ids, hap_map) {
  v <- vcfR::read.vcfR(merged_path, verbose = FALSE)
  gt <- vcfR::extract.gt(v, element = "GT", as.numeric = FALSE)
  fix <- as.data.frame(vcfR::getFIX(v), stringsAsFactors = FALSE)

  var_base_all <- paste0(fix$CHROM, ":", fix$POS)
  keep_idx <- var_base_all %in% common_pos

  gt <- gt[keep_idx, , drop = FALSE]
  var_base_filt <- var_base_all[keep_idx]
  var_ids <- make.unique(var_base_filt)

  if (length(var_ids) != length(vcf_merged_ids) || any(var_ids != vcf_merged_ids)) {
    stop("Variant ID mismatch between genotype matrix and cached CrossHap IDs")
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

  m01_all[in_both, , drop = FALSE]
}

make_group_heatmap_plot <- function(group_id, rep_id, elite_ids, m01_all, keep_cols, target_n_rows = NULL) {
  if (sum(keep_cols) == 0) {
    return(NULL)
  }

  rep_mat <- m01_all[rep_id, keep_cols, drop = FALSE]
  elite_mat <- NULL
  if (!is.null(elite_mat) && nrow(elite_mat) > 0) {
    elite_mat <- elite_mat[, keep_cols, drop = FALSE]
  }

  if (length(elite_ids) > 0) {
    elite_ids <- sort(intersect(elite_ids, rownames(m01_all)))
    if (length(elite_ids) > 0) {
      elite_mat <- m01_all[elite_ids, keep_cols, drop = FALSE]
    }
  }

  row_blocks <- list(rep_mat)
  row_names <- c(paste0("Haplotype group ", group_id))

  if (!is.null(elite_mat) && nrow(elite_mat) > 0) {
    sep_mat <- matrix(2, nrow = 1, ncol = ncol(rep_mat))
    colnames(sep_mat) <- colnames(rep_mat)
    rownames(sep_mat) <- "__separator__"
    row_blocks[[length(row_blocks) + 1]] <- sep_mat
    row_blocks[[length(row_blocks) + 1]] <- elite_mat
    row_names <- c(row_names, "__separator__", elite_ids)
  }

  plot_mat <- do.call(rbind, row_blocks)
  rownames(plot_mat) <- row_names

  if (!is.null(target_n_rows) && nrow(plot_mat) < target_n_rows) {
    n_pad <- target_n_rows - nrow(plot_mat)
    pad_mat <- matrix(2, nrow = n_pad, ncol = ncol(plot_mat))
    colnames(pad_mat) <- colnames(plot_mat)
    rownames(pad_mat) <- paste0("__pad__", seq_len(n_pad))
    plot_mat <- rbind(plot_mat, pad_mat)
  }

  long_df <- do.call(rbind, lapply(seq_len(nrow(plot_mat)), function(i) {
    vals <- plot_mat[i, ]
    fill_key <- ifelse(
      is.na(vals), "Missing",
      ifelse(vals == 0, "Reference",
      ifelse(vals == 1, "Alternate", "Separator"))
    )

    data.frame(
      row_label = rownames(plot_mat)[i],
      variant_index = seq_along(vals),
      fill_key = fill_key,
      stringsAsFactors = FALSE
    )
  }))

  long_df$fill_key <- factor(
    long_df$fill_key,
    levels = c("Reference", "Alternate", "Separator", "Missing")
  )

  long_df$row_label <- factor(long_df$row_label, levels = rev(rownames(plot_mat)))

  p <- ggplot(long_df, aes(x = variant_index, y = row_label, fill = fill_key)) +
    geom_tile(color = "white", linewidth = HEATMAP_TILE_BORDER, height = 0.88) +
    geom_point(
      data = data.frame(
        dummy_x = c(-1, -1, -1),
        dummy_y = c(NA, NA, NA),
        fill_key = factor(
          c("Reference", "Alternate", "Missing"),
          levels = c("Reference", "Alternate", "Separator", "Missing")
        )
      ),
      aes(x = dummy_x, y = dummy_y, fill = fill_key),
      inherit.aes = FALSE,
      alpha = 0,
      shape = 22,
      size = 0,
      stroke = 0,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = c(
        Reference = "#d8cca0",
        Alternate = "#2F4F4F",
        Separator = "#FFFFFF",
        Missing = "grey85"
      ),
      breaks = c("Reference", "Alternate", "Missing"),
      labels = c("REF", "ALT", "Missing"),
      drop = FALSE
    ) +
    guides(fill = guide_legend(override.aes = list(
      fill = c("#d8cca0", "#2F4F4F", "grey85"),
      color = c("grey40", "grey40", "grey40"),
      alpha = 1,
      shape = c(22, 22, 22),
      size = c(5, 5, 5),
      stroke = c(0.6, 0.6, 0.6)
    ))) +
    labs(
      x = NULL,
      y = NULL
    ) +
    scale_y_discrete(labels = function(x) ifelse(grepl("^__separator__$|^__pad__", x) | is.na(x), "", x)) +
    theme_minimal(base_size = HEATMAP_BASE_SIZE) +
    theme(
      axis.text.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 14, face = "bold")
    )

  list(
    plot = p,
    n_rows = nrow(plot_mat)
  )
}

# ==============================================================================
# MAIN
# ==============================================================================

gene_name <- make_gene_name(GENE_FILE)

cache_dir <- file.path(RUN_DIR, "Cache", TRAIT, gene_name)
hap_rds <- file.path(cache_dir, "HapObject.rds")
ctx_rds <- file.path(cache_dir, "RunContext.rds")

if (!file.exists(hap_rds)) stop("Missing HapObject.rds: ", hap_rds)
if (!file.exists(ctx_rds)) stop("Missing RunContext.rds: ", ctx_rds)

HapObject <- readRDS(hap_rds)
ctx <- readRDS(ctx_rds)

target <- list(
  trait = TRAIT,
  gene_file = GENE_FILE,
  merged_vcf_file = basename(ctx$merged_path),
  title = GENE_TITLE,
  epsilon = EPSILON,
  MGmin = MGMIN
)

pheno <- ctx$pheno
if (is.null(pheno) || nrow(pheno) == 0) {
  pheno <- load_pheno_wild_elite(ctx$pheno_file)
}

label <- paste0("Haplotypes_MGmin", target$MGmin, "_E", target$epsilon)
if (is.null(HapObject[[label]])) {
  stop("HapObject does not contain label: ", label)
}

hap_map <- HapObject[[label]]$Indfile %>%
  select(Ind, hap) %>%
  distinct() %>%
  mutate(
    hap = as.character(hap),
    Population = ifelse(grepl("^HS", Ind), "Wild", "Elite")
  )

m01_all <- extract_genotype_matrix(
  merged_path = ctx$merged_path,
  common_pos = ctx$common_pos,
  vcf_merged_ids = ctx$vcf_merged_ids,
  hap_map = hap_map
)

violin_obj <- build_violin_plot(
  HapObject = HapObject,
  target = target,
  pheno = pheno,
  reference_hap_id = as.character(REFERENCE_HAP_ID),
  exclude_hap0 = EXCLUDE_HAP0_FROM_VIOLIN
)

heatmap_haps <- order_haps(hap_map$hap)
if (!INCLUDE_HAP0_IN_HEATMAP) {
  heatmap_haps <- setdiff(heatmap_haps, "0")
}

group_rep_info <- lapply(heatmap_haps, function(hap_id) {
  get_group_representative_info(
    group_id = hap_id,
    hap_map = hap_map,
    m01_all = m01_all
  )
})

group_rep_info <- Filter(Negate(is.null), group_rep_info)

if (length(group_rep_info) == 0) {
  stop("No group representatives could be selected")
}

# Show all SNPs for the gene. Representatives are still chosen by elite similarity
# plus low missingness, but missing tiles are allowed in the final presentation.
all_gene_keep_cols <- rep(TRUE, ncol(m01_all))

group_plot_info <- lapply(group_rep_info, function(info) {
  make_group_heatmap_plot(
    group_id = info$group_id,
    rep_id = info$rep_id,
    elite_ids = info$elite_ids,
    m01_all = m01_all,
    keep_cols = all_gene_keep_cols,
    target_n_rows = NULL
  )
})

group_plot_info <- Filter(Negate(is.null), group_plot_info)

if (length(group_plot_info) == 0) {
  stop("No group heatmaps could be created")
}

group_plots <- lapply(seq_along(group_plot_info), function(i) {
  p <- group_plot_info[[i]]$plot
  if (i == length(group_plot_info)) {
    p + theme(legend.position = "right")
  } else {
    p + theme(legend.position = "none")
  }
})

group_heights <- vapply(group_plot_info, function(x) x$n_rows + 0.35, numeric(1))

heatmap_patch <- wrap_plots(group_plots, ncol = 1, heights = group_heights)

final_plot <- (
  violin_obj$plot /
  heatmap_patch
) +
  plot_layout(heights = c(3.0, 2.0))

ensure_dir(OUTPUT_DIR)

out_png <- file.path(
  OUTPUT_DIR,
  paste0(
    gene_name,
    "__MGmin", MGMIN,
    "__E", EPSILON,
    "__ref_", REFERENCE_HAP_ID,
    "__presentation.png"
  )
)

ggsave(
  filename = out_png,
  plot = final_plot,
  width = PNG_WIDTH_IN,
  height = PNG_HEIGHT_IN,
  dpi = PNG_DPI,
  bg = "white"
)

message("Saved: ", out_png)