suppressPackageStartupMessages({
  library(crosshap)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

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

extract_gene_name_from_out_pdf <- function(out_pdf) {
  b <- basename(out_pdf)
  b <- sub("\\.pdf$", "", b)
  m <- regexec("^CrosshapTree\\+Violin__(.+)__MGmin[0-9]+__EpsVec[0-9]+$", b)
  g <- regmatches(b, m)[[1]]
  if (length(g) >= 2) return(g[2])
  NA_character_
}

compute_kw_stats_for_label <- function(HapObject, label) {
  ind <- HapObject[[label]]$Indfile
  if (is.null(ind) || nrow(ind) == 0) {
    return(list(p = NA_real_, n_ind = 0L, n_groups = 0L, group_sizes = NA_character_))
  }

  ind2 <- ind %>%
    mutate(
      hap = as.character(hap),
      Pheno = suppressWarnings(as.numeric(Pheno))
    ) %>%
    filter(hap != "0") %>%
    filter(!is.na(Pheno))

  n_ind <- nrow(ind2)
  group_sizes <- as.integer(sort(table(ind2$hap)))
  n_groups <- length(group_sizes)

  if (n_ind == 0 || n_groups <= 1) {
    return(list(
      p = NA_real_,
      n_ind = as.integer(n_ind),
      n_groups = as.integer(n_groups),
      group_sizes = if (length(group_sizes) == 0) NA_character_ else paste(group_sizes, collapse = "|")
    ))
  }

  p <- tryCatch(
    kruskal.test(Pheno ~ hap, data = ind2)$p.value,
    error = function(e) NA_real_
  )

  list(
    p = p,
    n_ind = as.integer(n_ind),
    n_groups = as.integer(n_groups),
    group_sizes = paste(group_sizes, collapse = "|")
  )
}

add_header_to_crosshap_plot <- function(viz, header_title, header_subtitle) {
  if (inherits(viz, "ggplot")) {
    return(
      viz +
        labs(title = header_title, subtitle = header_subtitle) +
        theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0),
          plot.subtitle = element_text(size = 10, hjust = 0)
        )
    )
  }

  if (inherits(viz, "patchwork") && requireNamespace("patchwork", quietly = TRUE)) {
    return(
      viz +
        patchwork::plot_annotation(
          title = header_title,
          subtitle = header_subtitle,
          theme = theme(
            plot.title = element_text(face = "bold", size = 12, hjust = 0),
            plot.subtitle = element_text(size = 10, hjust = 0)
          )
        )
    )
  }

  viz
}

write_crosshap_tree_page <- function(HapObject, gene_name, MGmin, epsilon_vector) {
  if (!requireNamespace("clustree", quietly = TRUE)) {
    plot.new()
    title(main = paste0(gene_name, " | CrossHap tree (missing clustree package)"))
    text(
      0.02,
      0.90,
      paste0(
        "Install clustree to enable this page:\n",
        "install.packages('clustree')\n\n",
        "MGmin=", MGmin, " | eps=", paste(epsilon_vector, collapse = ",")
      ),
      adj = c(0, 1)
    )
    return(invisible(NULL))
  }

  tryCatch({
    p <- clustree::clustree_viz(HapObject, type = "hap") +
      ggtitle(
        paste0(gene_name, " | CrossHap tree"),
        subtitle = paste0("MGmin=", MGmin, " | eps=", paste(epsilon_vector, collapse = ","))
      )
    print(p)
  }, error = function(e) {
    plot.new()
    title(main = paste0(gene_name, " | CrossHap tree failed"))
    text(0.02, 0.90, conditionMessage(e), adj = c(0, 1))
  })

  invisible(NULL)
}

write_kw_summary_page <- function(title, trait, gene_file, MGmin, epsilon_vector, stats_by_eps) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  plot.new()
  title(main = paste0(title, " | Raw omnibus summary"))

  par(family = "mono")
  text(0.02, 0.96, paste0("Trait=", trait, " | Gene=", gene_file, " | MGmin=", MGmin), adj = c(0, 1), cex = 0.95)
  text(0.02, 0.91, paste0("eps vector=", paste(epsilon_vector, collapse = ",")), adj = c(0, 1), cex = 0.92)

  header <- sprintf("%-8s %-10s %-6s %-6s %-18s", "eps", "raw_p", "n", "groups", "group_sizes")
  text(0.02, 0.84, header, adj = c(0, 1), cex = 0.92)

  y0 <- 0.80
  dy <- 0.043

  for (i in seq_along(epsilon_vector)) {
    eps <- epsilon_vector[[i]]
    st <- stats_by_eps[[as.character(eps)]]
    line <- sprintf(
      "%-8s %-10s %-6d %-6d %-18s",
      as.character(eps),
      format_p_plain(st$p),
      st$n_ind,
      st$n_groups,
      st$group_sizes %||% "NA"
    )
    text(0.02, y0 - (i - 1) * dy, line, adj = c(0, 1), cex = 0.90)
  }
}

normalize_mgmin_test_stats <- function(mgmin_test_stats) {
  if (is.null(mgmin_test_stats)) return(NULL)

  df <- as.data.frame(mgmin_test_stats, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(df)

  df$epsilon <- as.character(df$epsilon)

  if ("MGmin" %in% names(df)) {
    df$MGmin <- as.integer(df$MGmin)
  }

  if ("trait" %in% names(df)) {
    df$trait <- as.character(df$trait)
  }

  if ("gene_file" %in% names(df)) {
    df$gene_file <- as.character(df$gene_file)
  }

  df
}

write_corrected_test_summary_page <- function(title, trait, gene_file, MGmin, epsilon_vector, stats_by_eps, mgmin_test_stats = NULL) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  plot.new()
  title(main = paste0(title, " | Corrected test summary"))

  mgmin_test_stats <- normalize_mgmin_test_stats(mgmin_test_stats)

  par(family = "mono")
  text(0.02, 0.96, paste0("Trait=", trait, " | Gene=", gene_file, " | MGmin=", MGmin), adj = c(0, 1), cex = 0.95)
  text(0.02, 0.91, "Within-gene Holm values are shown only for unique retained tests in the final statistics workflow.", adj = c(0, 1), cex = 0.88)

  header <- sprintf("%-8s %-10s %-10s %-8s %-18s", "eps", "raw_p", "holm_p", "groups", "group_sizes")
  text(0.02, 0.84, header, adj = c(0, 1), cex = 0.90)

  y0 <- 0.80
  dy <- 0.043

  for (i in seq_along(epsilon_vector)) {
    eps_chr <- as.character(epsilon_vector[[i]])
    st <- stats_by_eps[[eps_chr]]
    holm_p <- NA_real_

    if (!is.null(mgmin_test_stats) && nrow(mgmin_test_stats) > 0) {
      match_row <- mgmin_test_stats[
        mgmin_test_stats$trait == trait &
          mgmin_test_stats$gene_file == gene_file &
          as.integer(mgmin_test_stats$MGmin) == as.integer(MGmin) &
          mgmin_test_stats$epsilon == eps_chr,
        ,
        drop = FALSE
      ]

      if (nrow(match_row) > 1) {
        stop(
          paste0(
            "Corrected summary matched multiple rows for trait=",
            trait,
            ", gene_file=",
            gene_file,
            ", MGmin=",
            MGmin,
            ", epsilon=",
            eps_chr
          )
        )
      }

      if (nrow(match_row) == 1) {
        holm_p <- as.numeric(match_row$kw_p_holm_within_gene[1])
      }
    }

    line <- sprintf(
      "%-8s %-10s %-10s %-8s %-18s",
      eps_chr,
      format_p_plain(st$p),
      format_p_plain(holm_p),
      as.character(st$n_groups),
      st$group_sizes %||% "NA"
    )
    text(0.02, y0 - (i - 1) * dy, line, adj = c(0, 1), cex = 0.88)
  }
}

write_gene_level_summary_page <- function(title, trait, gene_file, gene_summary_row = NULL) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  plot.new()
  title(main = paste0(title, " | Gene-level corrected summary"))

  if (is.null(gene_summary_row) || nrow(gene_summary_row) == 0) {
    text(0.02, 0.95, "No final gene-level corrected summary is available for this gene.", adj = c(0, 1))
    return(invisible(NULL))
  }

  row <- as.data.frame(gene_summary_row, stringsAsFactors = FALSE)

  row <- row[
    row$trait == trait & row$gene_file == gene_file,
    ,
    drop = FALSE
  ]

  if (nrow(row) != 1) {
    stop(
      paste0(
        "Gene-level corrected summary matched ",
        nrow(row),
        " rows for trait=", trait,
        ", gene_file=", gene_file
      )
    )
  }

  row <- row[1, , drop = FALSE]
  lines <- c(
    paste0("Trait: ", trait),
    paste0("Gene: ", gene_file),
    paste0("Selected best MGmin: ", row$best_MGmin),
    paste0("Selected best epsilon: ", row$best_epsilon),
    paste0("Best raw p: ", format_p_plain(row$best_kw_p_raw)),
    paste0("Best within-gene Holm p: ", format_p_plain(row$best_kw_p_holm_within_gene)),
    paste0("Trait-level FDR p: ", format_p_plain(row$trait_fdr_p)),
    paste0("Trait-level Bonferroni p: ", format_p_plain(row$trait_bonferroni_p)),
    paste0("significant_fdr: ", as.character(row$significant_fdr)),
    paste0("significant_bonferroni: ", as.character(row$significant_bonferroni))
  )

  y <- 0.95
  for (line in lines) {
    text(0.02, y, line, adj = c(0, 1), cex = 0.95)
    y <- y - 0.075
  }

  invisible(NULL)
}

compute_pairwise_reference_tests <- function(plot_data) {
  groups <- sort(unique(plot_data$hap_label))
  if (length(groups) <= 1) return(NULL)

  global_mean <- mean(plot_data$Pheno, na.rm = TRUE)

  group_means <- plot_data %>%
    group_by(hap_label) %>%
    summarize(m = mean(Pheno, na.rm = TRUE), .groups = "drop") %>%
    mutate(diff = abs(m - global_mean))

  ref_group_label <- group_means$hap_label[which.max(group_means$diff)]
  other_groups <- groups[groups != ref_group_label]

  if (length(other_groups) == 0) return(NULL)

  raw_p <- vapply(other_groups, function(g) {
    tryCatch(
      wilcox.test(
        plot_data$Pheno[plot_data$hap_label == ref_group_label],
        plot_data$Pheno[plot_data$hap_label == g],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  }, numeric(1))

  p_adj <- p.adjust(raw_p, method = "holm")
  sig <- vapply(p_adj, p_to_signif_symbol, character(1))

  y_max <- max(plot_data$Pheno, na.rm = TRUE)
  y_min <- min(plot_data$Pheno, na.rm = TRUE)
  y_span <- max(1e-8, y_max - y_min)
  y_positions <- y_max + y_span * (0.10 + seq_along(other_groups) * 0.08)

  data.frame(
    group1 = ref_group_label,
    group2 = other_groups,
    p = raw_p,
    p.adj = p_adj,
    p.adj.signif = sig,
    y.position = y_positions,
    stringsAsFactors = FALSE
  )
}

plot_violin_with_holm <- function(HapObject, label, title_text) {
  plot_data <- HapObject[[label]]$Indfile
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    return(ggplot() + theme_void() + ggtitle(title_text, subtitle = "No Indfile data"))
  }

  plot_data <- plot_data %>%
    mutate(
      hap = as.character(hap),
      Pheno = suppressWarnings(as.numeric(Pheno))
    ) %>%
    filter(hap != "0") %>%
    filter(!is.na(Pheno))

  groups <- sort(unique(plot_data$hap))
  if (length(groups) <= 1) {
    return(ggplot() + theme_void() + ggtitle(title_text, subtitle = "Not enough groups after filtering hap!=0"))
  }

  fit <- kruskal.test(Pheno ~ hap, data = plot_data)
  subtitle_text <- paste0(
    "Kruskal-Wallis P (Nominal) = ",
    format_p_plain(fit$p.value),
    " | Pairwise: Wilcoxon reference-vs-others with Holm correction"
  )

  n_table <- table(plot_data$hap)
  plot_data$hap_label <- paste0(plot_data$hap, "\n(n=", n_table[plot_data$hap], ")")
  plot_data$hap_label <- factor(plot_data$hap_label, levels = unique(plot_data$hap_label))

  pairwise_df <- compute_pairwise_reference_tests(plot_data)

  p <- ggpubr::ggviolin(
    plot_data,
    x = "hap_label",
    y = "Pheno",
    fill = "hap_label",
    palette = "npg",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Haplotype Group",
      y = "Phenotype Value"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0, face = "italic")
    )

  if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      pairwise_df,
      label = "p.adj.signif",
      tip.length = 0.01,
      hide.ns = FALSE
    )
  }

  p
}

write_combined_pdf <- function(
  HapObject,
  out_pdf,
  title,
  trait,
  gene_file,
  gene_name = NULL,
  MGmin,
  epsilon_vector,
  mgmin_test_stats = NULL,
  gene_summary_row = NULL
) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(out_pdf, width = 10, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (is.null(gene_name) || !nzchar(gene_name)) {
    gene_name <- extract_gene_name_from_out_pdf(out_pdf)
  }
  if (is.na(gene_name) || !nzchar(gene_name)) {
    gene_name <- title
  }

  stats_by_eps <- setNames(vector("list", length(epsilon_vector)), as.character(epsilon_vector))
  for (eps in epsilon_vector) {
    label <- paste0("Haplotypes_MGmin", MGmin, "_E", eps)
    if (is.null(HapObject[[label]])) {
      stats_by_eps[[as.character(eps)]] <- list(
        p = NA_real_,
        n_ind = 0L,
        n_groups = 0L,
        group_sizes = NA_character_
      )
    } else {
      stats_by_eps[[as.character(eps)]] <- compute_kw_stats_for_label(HapObject, label)
    }
  }

  write_crosshap_tree_page(
    HapObject = HapObject,
    gene_name = gene_name,
    MGmin = MGmin,
    epsilon_vector = epsilon_vector
  )

  write_kw_summary_page(
    title = title,
    trait = trait,
    gene_file = gene_file,
    MGmin = MGmin,
    epsilon_vector = epsilon_vector,
    stats_by_eps = stats_by_eps
  )

  write_corrected_test_summary_page(
    title = title,
    trait = trait,
    gene_file = gene_file,
    MGmin = MGmin,
    epsilon_vector = epsilon_vector,
    stats_by_eps = stats_by_eps,
    mgmin_test_stats = mgmin_test_stats
  )

  write_gene_level_summary_page(
    title = title,
    trait = trait,
    gene_file = gene_file,
    gene_summary_row = gene_summary_row
  )

  for (eps in epsilon_vector) {
    label <- paste0("Haplotypes_MGmin", MGmin, "_E", eps)
    if (is.null(HapObject[[label]])) next

    st <- stats_by_eps[[as.character(eps)]]

    tryCatch({
      viz <- suppressWarnings(suppressMessages(crosshap::crosshap_viz(HapObject, epsilon = eps)))

      header_title <- paste0(gene_name, " | eps=", eps, " | MGmin=", MGmin)
      header_subtitle <- paste0(
        title,
        " | Kruskal-Wallis P (Nominal) = ", format_p_plain(st$p),
        " | n=", st$n_ind,
        " | groups=", st$n_groups,
        " | sizes=", st$group_sizes %||% "NA"
      )

      viz2 <- add_header_to_crosshap_plot(viz, header_title, header_subtitle)
      print(viz2)
    }, error = function(e) {
      plot.new()
      title(main = paste0("crosshap_viz failed | eps=", eps))
      text(0.02, 0.95, conditionMessage(e), adj = c(0, 1))
    })

    p_vio <- plot_violin_with_holm(
      HapObject = HapObject,
      label = label,
      title_text = paste0(title, " | eps=", eps)
    )
    print(p_vio)
  }

  invisible(out_pdf)
}