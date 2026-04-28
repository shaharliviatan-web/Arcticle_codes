suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
})

build_violin_components <- function(HapObject, target, pheno, pairwise_method = "holm") {
  label <- paste0("Haplotypes_MGmin", target$MGmin, "_E", target$epsilon)
  indfile <- HapObject[[label]]$Indfile

  if (is.null(indfile) || nrow(indfile) == 0) {
    stop("No Indfile data for ", label)
  }

  stats_df <- as.data.frame(indfile, stringsAsFactors = FALSE)
  if (!"Pheno" %in% names(stats_df)) {
    stats_df <- dplyr::left_join(stats_df, pheno, by = "Ind")
  }

  stats_df$hap <- as.character(stats_df$hap)
  stats_df$Pheno <- suppressWarnings(as.numeric(stats_df$Pheno))

  plot_data <- stats_df %>% filter(hap != "0")
  wild_data <- plot_data %>% filter(!is.na(Pheno))
  elite_data <- plot_data %>% filter(is.na(Pheno))

  if (nrow(wild_data) == 0) {
    stop("No wild rows after filtering hap!=0")
  }

  n_table <- wild_data %>% count(hap, name = "n_wild")
  plot_data <- plot_data %>%
    left_join(n_table, by = "hap") %>%
    mutate(hap_label = paste0(hap, "\n(n=", n_wild, ")"))

  wild_data <- wild_data %>% left_join(plot_data %>% select(Ind, hap_label), by = "Ind")
  elite_data <- elite_data %>% left_join(plot_data %>% select(Ind, hap_label), by = "Ind")

  wild_groups <- unique(wild_data$hap_label)
  if (length(wild_groups) <= 1) {
    stop("Not enough groups after filtering hap!=0")
  }

  kw_fit <- kruskal.test(Pheno ~ hap, data = wild_data)

  ref_group_label <- NULL
  if (!is.null(target$force_ref)) {
    match_idx <- grep(paste0("^", target$force_ref), wild_groups)
    if (length(match_idx) > 0) ref_group_label <- wild_groups[match_idx[1]]
  }

  if (is.null(ref_group_label)) {
    match_idx <- grep("^A\n", wild_groups)
    if (length(match_idx) > 0) {
      ref_group_label <- wild_groups[match_idx[1]]
    } else {
      ref_group_label <- wild_groups[1]
    }
  }

  other_groups <- wild_groups[wild_groups != ref_group_label]
  pairwise_tbl <- data.frame()
  if (length(other_groups) > 0) {
    pairwise_tbl <- do.call(rbind, lapply(other_groups, function(group_label) {
      x <- wild_data$Pheno[wild_data$hap_label == ref_group_label]
      y <- wild_data$Pheno[wild_data$hap_label == group_label]
      p_raw <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
      data.frame(group1 = ref_group_label, group2 = group_label, p = p_raw, stringsAsFactors = FALSE)
    }))
    pairwise_tbl$p.adj <- p.adjust(pairwise_tbl$p, method = pairwise_method)
    pairwise_tbl$p.signif <- vapply(pairwise_tbl$p.adj, p_to_signif_symbol, character(1))

    y_top <- max(wild_data$Pheno, na.rm = TRUE)
    y_range <- y_top - min(wild_data$Pheno, na.rm = TRUE)
    if (!is.finite(y_range) || y_range == 0) y_range <- 1
    pairwise_tbl$y.position <- y_top + seq_len(nrow(pairwise_tbl)) * (0.08 * y_range)
  }

  y_range <- max(wild_data$Pheno, na.rm = TRUE) - min(wild_data$Pheno, na.rm = TRUE)
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  dummy_y <- min(wild_data$Pheno, na.rm = TRUE) - (y_range * 0.15)
  elite_data$DummyY <- dummy_y

  subtitle_text <- paste0(
    "Kruskal-Wallis P (wild, nominal) = ", formatC(kw_fit$p.value, format = "e", digits = 2),
    " | Pairwise symbols = Holm-corrected"
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
      size = 3.5,
      color = "#2ca25f",
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
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0, face = "italic")
    )

  if (nrow(pairwise_tbl) > 0) {
    p_violin <- p_violin + ggpubr::stat_pvalue_manual(
      pairwise_tbl,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.3,
      size = 5
    )
  }

  elite_assignments <- elite_data %>%
    select(Elite_Line = Ind, Assigned_Haplotype = hap) %>%
    distinct() %>%
    arrange(Assigned_Haplotype, Elite_Line)

  list(
    plot = p_violin,
    kw_p = kw_fit$p.value,
    n_ind = nrow(wild_data),
    n_groups = length(unique(wild_data$hap)),
    group_sizes = paste(sort(as.integer(table(wild_data$hap))), collapse = "|"),
    elite_assignments = elite_assignments,
    ref_group_label = ref_group_label,
    pairwise_table = pairwise_tbl
  )
}

write_combined_pdf <- function(HapObject, target, pheno, out_pdf, pairwise_method = "holm") {
  ensure_dir(dirname(out_pdf))
  label <- paste0("Haplotypes_MGmin", target$MGmin, "_E", target$epsilon)

  components <- build_violin_components(HapObject, target, pheno, pairwise_method = pairwise_method)

  grDevices::pdf(out_pdf, width = 16, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot.new()
  title(main = paste0(target$title, " | Summary"))
  txt <- c(
    paste0("Trait: ", target$trait),
    paste0("Gene: ", target$gene_file),
    paste0("Merged VCF: ", target$merged_vcf_file),
    paste0("MGmin: ", target$MGmin),
    paste0("Epsilon: ", target$epsilon),
    paste0("Kruskal-Wallis P (wild, nominal): ", format_p_plain(components$kw_p)),
    paste0("Wild individuals tested: ", components$n_ind),
    paste0("Wild haplotype groups: ", components$n_groups),
    paste0("Wild group sizes: ", components$group_sizes),
    paste0("Reference group for pairwise tests: ", components$ref_group_label),
    paste0("Pairwise correction: ", pairwise_method)
  )
  text(0.02, 0.96, paste(txt, collapse = "\n"), adj = c(0, 1), family = "mono")

  tryCatch({
    viz <- crosshap_viz(HapObject, epsilon = target$epsilon)
    print(viz + plot_annotation(title = paste0(target$title, " | CrossHap"), subtitle = paste0("MGmin=", target$MGmin, " | epsilon=", target$epsilon)))
  }, error = function(e) {
    plot.new()
    title(main = paste0(target$title, " | CrossHap plot failed"))
    text(0.02, 0.96, conditionMessage(e), adj = c(0, 1))
  })

  print(components$plot)

  list(
    kw_p = components$kw_p,
    n_ind = components$n_ind,
    n_groups = components$n_groups,
    group_sizes = components$group_sizes,
    elite_assignments = components$elite_assignments,
    pairwise_table = components$pairwise_table,
    out_pdf = out_pdf,
    label = label
  )
}
