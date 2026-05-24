suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_path_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_path_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[[1]])))
} else {
  getwd()
}

results_dir <- if (length(args) >= 1) args[[1]] else "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_data/results_maf0.000001_geno1_thin0_pc5"
bim_file <- if (length(args) >= 2) args[[2]] else file.path(results_dir, "snp_map.bim")
output_base_dir <- if (length(args) >= 3) args[[3]] else script_dir
file_pattern <- if (length(args) >= 4) args[[4]] else "^morexV3_.*_gwas\\.ps$"

diagnostics_dir <- file.path(output_base_dir, "FDR_diagnostics")
top_snps_dir <- file.path(diagnostics_dir, "top_snps")
plots_dir <- file.path(diagnostics_dir, "plots")
ranked_top_dir <- file.path(plots_dir, "ranked_top")
ranked_logx_dir <- file.path(plots_dir, "ranked_logx")
pvalue_zooms_dir <- file.path(plots_dir, "pvalue_zooms")
qvalue_zooms_dir <- file.path(plots_dir, "qvalue_zooms")
qqplots_dir <- file.path(plots_dir, "qqplots")
summaries_dir <- file.path(diagnostics_dir, "summaries")

dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(top_snps_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ranked_top_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ranked_logx_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pvalue_zooms_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qvalue_zooms_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qqplots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

if (!file.exists(bim_file)) {
  stop("SNP map file not found: ", bim_file)
}

message("Loading SNP map: ", bim_file)

snp_map <- tryCatch(
  read.table(
    bim_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("CHR", "SNP_ID", "CM", "BP", "A1", "A2")
  ),
  error = function(error) {
    stop("Failed to read SNP map: ", conditionMessage(error))
  }
)

gwas_files <- sort(list.files(
  path = results_dir,
  pattern = file_pattern,
  full.names = TRUE
))

if (!length(gwas_files)) {
  stop("No GWAS result files matched pattern '", file_pattern, "' in ", results_dir)
}

message("Matched ", length(gwas_files), " GWAS result files.")

top_rank_limits <- c(1000, 10000, 100000)
pvalue_zoom_limits <- c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
qvalue_zoom_limits <- c(0.20, 0.10, 0.05)

extract_trait_name <- function(file_path) {
  trait_name <- basename(file_path)
  trait_name <- sub("_gwas\\.ps$", "", trait_name)
  trait_name <- sub("^morexV3_", "", trait_name)
  trait_name
}

build_trait_interpretation <- function(n_q_le_0_05, n_q_le_0_10, n_bonferroni_0_05, n_bonferroni_0_10) {
  if (n_q_le_0_05 > 0 && n_bonferroni_0_05 > 0) {
    return("This trait has SNPs passing BH-FDR q<=0.05 and Bonferroni 0.05.")
  }

  if (n_q_le_0_05 > 0) {
    return("This trait has SNPs passing BH-FDR q<=0.05 but not Bonferroni 0.05.")
  }

  if (n_q_le_0_10 > 0) {
    return("This trait has SNPs passing BH-FDR q<=0.10 but not Bonferroni.")
  }

  if (n_bonferroni_0_10 > 0) {
    return("This trait has SNPs passing Bonferroni 0.10 but not BH-FDR.")
  }

  "No SNPs pass BH-FDR or Bonferroni."
}

calc_lambda_gc <- function(p_values) {
  chisq_values <- qchisq(1 - p_values, df = 1)
  median(chisq_values, na.rm = TRUE) / qchisq(0.5, df = 1)
}

build_trait_stats <- function(valid_df, total_snps) {
  valid_snps <- nrow(valid_df)
  bonferroni_0_05_cutoff_raw_p <- 0.05 / valid_snps
  bonferroni_0_10_cutoff_raw_p <- 0.10 / valid_snps
  n_q_le_0_05 <- sum(valid_df$BH_q <= 0.05, na.rm = TRUE)
  n_q_le_0_10 <- sum(valid_df$BH_q <= 0.10, na.rm = TRUE)
  n_bonferroni_0_05 <- sum(valid_df$P <= bonferroni_0_05_cutoff_raw_p, na.rm = TRUE)
  n_bonferroni_0_10 <- sum(valid_df$P <= bonferroni_0_10_cutoff_raw_p, na.rm = TRUE)
  bh_q_0_05_cutoff_raw_p <- if (n_q_le_0_05 > 0) max(valid_df$P[valid_df$BH_q <= 0.05], na.rm = TRUE) else NA_real_
  bh_q_0_10_cutoff_raw_p <- if (n_q_le_0_10 > 0) max(valid_df$P[valid_df$BH_q <= 0.10], na.rm = TRUE) else NA_real_

  list(
    total_snps = total_snps,
    valid_snps = valid_snps,
    min_raw_p_value = min(valid_df$P, na.rm = TRUE),
    max_minus_log10_p = max(valid_df$log10P, na.rm = TRUE),
    min_bh_q_value = min(valid_df$BH_q, na.rm = TRUE),
    n_q_le_0_05 = n_q_le_0_05,
    n_q_le_0_10 = n_q_le_0_10,
    n_bonferroni_0_05 = n_bonferroni_0_05,
    n_bonferroni_0_10 = n_bonferroni_0_10,
    bonferroni_0_05_cutoff_raw_p = bonferroni_0_05_cutoff_raw_p,
    bonferroni_0_10_cutoff_raw_p = bonferroni_0_10_cutoff_raw_p,
    bonferroni_0_05_cutoff_minus_log10_p = -log10(bonferroni_0_05_cutoff_raw_p),
    bonferroni_0_10_cutoff_minus_log10_p = -log10(bonferroni_0_10_cutoff_raw_p),
    bh_q_0_05_cutoff_raw_p = bh_q_0_05_cutoff_raw_p,
    bh_q_0_10_cutoff_raw_p = bh_q_0_10_cutoff_raw_p,
    bh_q_0_05_cutoff_minus_log10_p = if (is.na(bh_q_0_05_cutoff_raw_p)) NA_real_ else -log10(bh_q_0_05_cutoff_raw_p),
    bh_q_0_10_cutoff_minus_log10_p = if (is.na(bh_q_0_10_cutoff_raw_p)) NA_real_ else -log10(bh_q_0_10_cutoff_raw_p),
    lambda_gc = calc_lambda_gc(valid_df$P)
  )
}

build_ranked_annotation_text <- function(stats, rank_limit = NULL) {
  rank_label <- if (is.null(rank_limit)) {
    paste("Valid SNPs =", format(stats$valid_snps, big.mark = ","))
  } else {
    paste("Top ranks shown =", format(rank_limit, big.mark = ","), "of", format(stats$valid_snps, big.mark = ","))
  }

  paste(
    rank_label,
    paste("q<=0.05 hits =", stats$n_q_le_0_05),
    paste("q<=0.10 hits =", stats$n_q_le_0_10),
    paste("Bonferroni 0.05 hits =", stats$n_bonferroni_0_05),
    paste("Min p =", signif(stats$min_raw_p_value, 4)),
    paste("Min BH q =", signif(stats$min_bh_q_value, 4)),
    sep = "\n"
  )
}

build_zoom_annotation_text <- function(lines) {
  paste(lines, collapse = "\n")
}

make_histogram_plot <- function(plot_data, x_var, title_text, x_label, output_file, binwidth = NULL, limits = NULL, vlines = NULL) {
  histogram_plot <- ggplot(plot_data, aes(x = .data[[x_var]])) +
    geom_histogram(
      bins = if (is.null(binwidth)) 50 else NULL,
      binwidth = binwidth,
      fill = "steelblue4",
      color = "white",
      linewidth = 0.2
    ) +
    labs(title = title_text, x = x_label, y = "SNP count") +
    theme_bw(base_size = 12)

  if (!is.null(limits)) {
    histogram_plot <- histogram_plot + coord_cartesian(xlim = limits)
  }

  if (!is.null(vlines)) {
    for (vline in vlines) {
      histogram_plot <- histogram_plot +
        geom_vline(xintercept = vline$value, color = vline$color, linetype = vline$linetype, linewidth = 0.7)
    }
  }

  ggsave(output_file, histogram_plot, width = 8, height = 5, dpi = 300)
}

make_zoom_histogram_plot <- function(plot_data, x_var, title_text, x_label, output_file, range_limit, annotation_text, vlines = NULL) {
  binwidth <- if (range_limit > 0) range_limit / 80 else NULL

  histogram_plot <- ggplot(plot_data, aes(x = .data[[x_var]])) +
    geom_histogram(
      bins = if (is.null(binwidth)) 80 else NULL,
      binwidth = binwidth,
      fill = "steelblue4",
      color = "white",
      linewidth = 0.2
    ) +
    labs(title = title_text, x = x_label, y = "SNP count") +
    coord_cartesian(xlim = c(0, range_limit)) +
    theme_bw(base_size = 12) +
    annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.02, vjust = 1.1, size = 3.7)

  if (!is.null(vlines)) {
    for (vline in vlines) {
      histogram_plot <- histogram_plot +
        geom_vline(
          xintercept = vline$value,
          color = vline$color,
          linetype = vline$linetype,
          linewidth = 0.7
        )
    }
  }

  ggsave(output_file, histogram_plot, width = 8, height = 5, dpi = 300)
}

make_ranked_bh_plot_core <- function(ranked_df, stats, trait_name, output_file, title_suffix = NULL, use_log_rank = FALSE, rank_limit = NULL) {
  max_y <- max(
    ranked_df$observed_log10p,
    ranked_df$bh_005,
    ranked_df$bh_010,
    stats$bonferroni_0_05_cutoff_minus_log10_p,
    stats$bonferroni_0_10_cutoff_minus_log10_p,
    na.rm = TRUE
  )

  x_values <- if (use_log_rank) log10(ranked_df$rank) else ranked_df$rank
  x_label <- if (use_log_rank) expression(log[10]("SNP rank")) else "SNP rank"
  title_text <- paste0(
    "Ranked BH-FDR Diagnostic: ",
    trait_name,
    if (!is.null(title_suffix)) paste0(" (", title_suffix, ")") else ""
  )

  png(output_file, width = 8, height = 5, units = "in", res = 300)
  par(mar = c(5, 5, 4, 7) + 0.1)
  plot(
    x_values,
    ranked_df$observed_log10p,
    type = "l",
    col = "black",
    lwd = 1,
    xlab = x_label,
    ylab = expression(-log[10](italic(p))),
    main = title_text,
    ylim = c(0, ceiling(max_y) + 1)
  )
  lines(x_values, ranked_df$bh_005, col = "firebrick", lwd = 1.2, lty = 2)
  lines(x_values, ranked_df$bh_010, col = "royalblue3", lwd = 1.2, lty = 2)
  abline(h = stats$bonferroni_0_05_cutoff_minus_log10_p, col = "darkgreen", lwd = 1.1, lty = 3)
  abline(h = stats$bonferroni_0_10_cutoff_minus_log10_p, col = "darkolivegreen4", lwd = 1.1, lty = 3)
  legend(
    "topright",
    inset = c(-0.3, 0),
    xpd = TRUE,
    legend = c(
      "Observed -log10(p)",
      "BH q=0.05 ranked threshold",
      "BH q=0.10 ranked threshold",
      "Bonferroni 0.05 cutoff",
      "Bonferroni 0.10 cutoff"
    ),
    col = c("black", "firebrick", "royalblue3", "darkgreen", "darkolivegreen4"),
    lty = c(1, 2, 2, 3, 3),
    lwd = c(1, 1.2, 1.2, 1.1, 1.1),
    bty = "n",
    cex = 0.85
  )
  mtext(build_ranked_annotation_text(stats, rank_limit = rank_limit), side = 4, line = 0.5, adj = 0, cex = 0.8)
  dev.off()
}

make_ranked_bh_plot <- function(valid_df, trait_name, output_file) {
  ranked_df <- valid_df[order(valid_df$P, decreasing = FALSE), c("P"), drop = FALSE]
  ranked_df$rank <- seq_len(nrow(ranked_df))
  ranked_df$observed_log10p <- -log10(ranked_df$P)
  ranked_df$bh_005 <- -log10((ranked_df$rank / nrow(ranked_df)) * 0.05)
  ranked_df$bh_010 <- -log10((ranked_df$rank / nrow(ranked_df)) * 0.10)

  make_ranked_bh_plot_core(
    ranked_df = ranked_df,
    stats = attr(valid_df, "trait_stats"),
    trait_name = trait_name,
    output_file = output_file
  )
}

make_ranked_top_plots <- function(valid_df, trait_name, output_dir, rank_limits) {
  ranked_df <- valid_df[order(valid_df$P, decreasing = FALSE), c("P"), drop = FALSE]
  ranked_df$rank <- seq_len(nrow(ranked_df))
  ranked_df$observed_log10p <- -log10(ranked_df$P)
  ranked_df$bh_005 <- -log10((ranked_df$rank / nrow(ranked_df)) * 0.05)
  ranked_df$bh_010 <- -log10((ranked_df$rank / nrow(ranked_df)) * 0.10)

  for (rank_limit in rank_limits) {
    current_limit <- min(rank_limit, nrow(ranked_df))
    subset_df <- ranked_df[seq_len(current_limit), , drop = FALSE]
    output_file <- file.path(output_dir, paste0(trait_name, "_ranked_top_", format(current_limit, scientific = FALSE, trim = TRUE), ".png"))
    make_ranked_bh_plot_core(
      ranked_df = subset_df,
      stats = attr(valid_df, "trait_stats"),
      trait_name = trait_name,
      output_file = output_file,
      title_suffix = paste0("top ", format(current_limit, big.mark = ","), " ranks"),
      rank_limit = current_limit
    )
  }
}

make_ranked_logx_plot <- function(valid_df, trait_name, output_file) {
  ranked_df <- valid_df[order(valid_df$P, decreasing = FALSE), c("P"), drop = FALSE]
  ranked_df$rank <- seq_len(nrow(ranked_df))
  ranked_df$observed_log10p <- -log10(ranked_df$P)
  ranked_df$bh_005 <- -log10((ranked_df$rank / nrow(ranked_df)) * 0.05)
  ranked_df$bh_010 <- -log10((ranked_df$rank / nrow(ranked_df)) * 0.10)

  make_ranked_bh_plot_core(
    ranked_df = ranked_df,
    stats = attr(valid_df, "trait_stats"),
    trait_name = trait_name,
    output_file = output_file,
    title_suffix = "log10 rank axis",
    use_log_rank = TRUE
  )
}

make_pvalue_zoom_plots <- function(valid_df, trait_name, output_dir, zoom_limits) {
  stats <- attr(valid_df, "trait_stats")
  vlines <- list(
    list(value = stats$bonferroni_0_05_cutoff_raw_p, color = "darkgreen", linetype = "dashed"),
    list(value = stats$bonferroni_0_10_cutoff_raw_p, color = "darkolivegreen4", linetype = "dashed")
  )

  if (!is.na(stats$bh_q_0_05_cutoff_raw_p)) {
    vlines <- c(vlines, list(list(value = stats$bh_q_0_05_cutoff_raw_p, color = "firebrick", linetype = "dotted")))
  }

  if (!is.na(stats$bh_q_0_10_cutoff_raw_p)) {
    vlines <- c(vlines, list(list(value = stats$bh_q_0_10_cutoff_raw_p, color = "royalblue3", linetype = "dotted")))
  }

  for (zoom_limit in zoom_limits) {
    zoom_df <- valid_df[valid_df$P <= zoom_limit, c("P"), drop = FALSE]
    annotation_text <- build_zoom_annotation_text(c(
      paste("SNPs in range =", format(nrow(zoom_df), big.mark = ",")),
      paste("Bonferroni 0.05 =", signif(stats$bonferroni_0_05_cutoff_raw_p, 4)),
      paste("Bonferroni 0.10 =", signif(stats$bonferroni_0_10_cutoff_raw_p, 4)),
      if (!is.na(stats$bh_q_0_05_cutoff_raw_p)) paste("Observed BH-FDR q<=0.05 cutoff =", signif(stats$bh_q_0_05_cutoff_raw_p, 4)) else "Observed BH-FDR q<=0.05 cutoff = none",
      if (!is.na(stats$bh_q_0_10_cutoff_raw_p)) paste("Observed BH-FDR q<=0.10 cutoff =", signif(stats$bh_q_0_10_cutoff_raw_p, 4)) else "Observed BH-FDR q<=0.10 cutoff = none"
    ))

    output_file <- file.path(output_dir, paste0(trait_name, "_pvalue_zoom_le_", format(zoom_limit, scientific = TRUE), ".png"))

    make_zoom_histogram_plot(
      plot_data = zoom_df,
      x_var = "P",
      title_text = paste("Raw P-value Zoom:", trait_name, "(p <=", format(zoom_limit, scientific = TRUE), ")"),
      x_label = "Raw p-value",
      output_file = output_file,
      range_limit = zoom_limit,
      annotation_text = annotation_text,
      vlines = vlines
    )
  }
}

make_qvalue_zoom_plots <- function(valid_df, trait_name, output_dir, zoom_limits) {
  for (zoom_limit in zoom_limits) {
    zoom_df <- valid_df[valid_df$BH_q <= zoom_limit, c("BH_q"), drop = FALSE]
    annotation_text <- build_zoom_annotation_text(c(
      paste("SNPs in range =", format(nrow(zoom_df), big.mark = ",")),
      paste("Valid SNPs =", format(nrow(valid_df), big.mark = ","))
    ))

    output_file <- file.path(output_dir, paste0(trait_name, "_qvalue_zoom_le_", gsub("\\.", "_", sprintf("%.2f", zoom_limit)), ".png"))

    make_zoom_histogram_plot(
      plot_data = zoom_df,
      x_var = "BH_q",
      title_text = paste("BH Q-value Zoom:", trait_name, "(q <=", zoom_limit, ")"),
      x_label = "BH q-value",
      output_file = output_file,
      range_limit = zoom_limit,
      annotation_text = annotation_text,
      vlines = list(
        list(value = 0.05, color = "firebrick", linetype = "dashed"),
        list(value = 0.10, color = "royalblue3", linetype = "dashed")
      )
    )
  }
}

make_qq_plot <- function(valid_df, trait_name, output_file) {
  sorted_p <- sort(valid_df$P)
  expected <- -log10(ppoints(length(sorted_p)))
  observed <- -log10(sorted_p)
  max_axis <- ceiling(max(expected, observed, na.rm = TRUE) * 1.05)

  png(output_file, width = 8, height = 8, units = "in", res = 300)
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(
    expected,
    observed,
    pch = 20,
    cex = 0.25,
    col = "navy",
    xlab = "Expected -log10(p)",
    ylab = "Observed -log10(p)",
    main = paste("QQ Plot:", trait_name),
    xlim = c(0, max_axis),
    ylim = c(0, max_axis)
  )
  abline(0, 1, col = "red3", lty = 2, lwd = 1.5)
  legend(
    "topleft",
    legend = c(
      paste("Valid SNPs:", format(nrow(valid_df), big.mark = ",")),
      paste("Lambda GC:", format(round(attr(valid_df, "trait_stats")$lambda_gc, 3), nsmall = 3))
    ),
    bty = "n"
  )
  dev.off()
}

write_trait_interpretation <- function(trait_name, stats, output_file) {
  interpretation_df <- data.frame(
    trait_name = trait_name,
    min_raw_p = stats$min_raw_p_value,
    min_bh_q = stats$min_bh_q_value,
    n_q_le_0_05 = stats$n_q_le_0_05,
    n_q_le_0_10 = stats$n_q_le_0_10,
    n_bonferroni_0_05 = stats$n_bonferroni_0_05,
    n_bonferroni_0_10 = stats$n_bonferroni_0_10,
    observed_bh_q_0_05_cutoff_raw_p = stats$bh_q_0_05_cutoff_raw_p,
    observed_bh_q_0_05_cutoff_minus_log10_p = stats$bh_q_0_05_cutoff_minus_log10_p,
    observed_bh_q_0_10_cutoff_raw_p = stats$bh_q_0_10_cutoff_raw_p,
    observed_bh_q_0_10_cutoff_minus_log10_p = stats$bh_q_0_10_cutoff_minus_log10_p,
    automatic_interpretation = build_trait_interpretation(
      stats$n_q_le_0_05,
      stats$n_q_le_0_10,
      stats$n_bonferroni_0_05,
      stats$n_bonferroni_0_10
    ),
    stringsAsFactors = FALSE
  )

  write.csv(interpretation_df, output_file, row.names = FALSE)
}

process_gwas_file <- function(gwas_file, snp_map_df, top_snps_output_dir, plots_output_dir) {
  trait_name <- extract_trait_name(gwas_file)
  message("\n--- Processing trait: ", trait_name, " ---")

  gwas_results <- tryCatch(
    read.table(
      gwas_file,
      header = FALSE,
      stringsAsFactors = FALSE,
      col.names = c("SNP_ID_DUMMY", "BETA", "SE", "P")
    ),
    error = function(error) {
      warning("Failed to read ", gwas_file, ": ", conditionMessage(error))
      NULL
    }
  )

  if (is.null(gwas_results)) {
    return(NULL)
  }

  total_snps <- nrow(gwas_results)

  if (!total_snps) {
    warning("Skipping empty GWAS file: ", gwas_file)
    return(NULL)
  }

  if (total_snps != nrow(snp_map_df)) {
    warning(
      "Skipping trait ", trait_name,
      " because row count does not match SNP map. GWAS rows = ", total_snps,
      ", SNP map rows = ", nrow(snp_map_df)
    )
    return(NULL)
  }

  merged_df <- data.frame(
    SNP_ID = snp_map_df$SNP_ID,
    CHR = snp_map_df$CHR,
    BP = snp_map_df$BP,
    P = gwas_results$P,
    stringsAsFactors = FALSE
  )

  valid_mask <- is.finite(merged_df$P) & !is.na(merged_df$P) & merged_df$P > 0 & merged_df$P <= 1
  valid_df <- merged_df[valid_mask, , drop = FALSE]
  valid_snps <- nrow(valid_df)

  if (!valid_snps) {
    warning("Skipping trait ", trait_name, " because it has no valid p-values.")
    return(NULL)
  }

  valid_df$BH_q <- p.adjust(valid_df$P, method = "BH")
  valid_df$log10P <- -log10(valid_df$P)

  trait_stats <- build_trait_stats(valid_df, total_snps)
  attr(valid_df, "trait_stats") <- trait_stats

  top_n <- min(100, valid_snps)
  top_snps_df <- valid_df[order(valid_df$P, decreasing = FALSE), c("SNP_ID", "CHR", "BP", "P", "BH_q", "log10P")]
  top_snps_df <- head(top_snps_df, top_n)
  colnames(top_snps_df) <- c("SNP_ID", "CHR", "BP", "raw_P", "BH_q_value", "minus_log10_P")

  top_snps_file <- file.path(top_snps_output_dir, paste0("top100_", trait_name, ".csv"))
  write.csv(top_snps_df, top_snps_file, row.names = FALSE)

  plot_data <- data.frame(P = valid_df$P, BH_q = valid_df$BH_q)
  zoom_plot_data <- plot_data[plot_data$P <= 0.01, , drop = FALSE]

  make_histogram_plot(
    plot_data = plot_data,
    x_var = "P",
    title_text = paste("Raw P-value Histogram:", trait_name),
    x_label = "Raw p-value",
    output_file = file.path(plots_output_dir, paste0(trait_name, "_raw_p_histogram.png")),
    limits = c(0, 1)
  )

  if (nrow(zoom_plot_data)) {
    make_histogram_plot(
      plot_data = zoom_plot_data,
      x_var = "P",
      title_text = paste("Zoomed Raw P-value Histogram (<= 0.01):", trait_name),
      x_label = "Raw p-value",
      output_file = file.path(plots_output_dir, paste0(trait_name, "_raw_p_histogram_zoom_0.01.png")),
      binwidth = 0.00025,
      limits = c(0, 0.01)
    )
  } else {
    warning("Trait ", trait_name, " has no p-values <= 0.01; saving empty zoom plot frame.")
    empty_zoom_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No valid p-values <= 0.01", size = 5) +
      labs(
        title = paste("Zoomed Raw P-value Histogram (<= 0.01):", trait_name),
        x = "Raw p-value",
        y = "SNP count"
      ) +
      theme_void()
    ggsave(
      file.path(plots_output_dir, paste0(trait_name, "_raw_p_histogram_zoom_0.01.png")),
      empty_zoom_plot,
      width = 8,
      height = 5,
      dpi = 300
    )
  }

  make_histogram_plot(
    plot_data = plot_data,
    x_var = "BH_q",
    title_text = paste("BH Q-value Histogram:", trait_name),
    x_label = "BH q-value",
    output_file = file.path(plots_output_dir, paste0(trait_name, "_bh_q_histogram.png")),
    limits = c(0, 1),
    vlines = list(
      list(value = 0.05, color = "firebrick", linetype = "dashed"),
      list(value = 0.10, color = "royalblue3", linetype = "dashed")
    )
  )

  make_ranked_bh_plot(
    valid_df = valid_df,
    trait_name = trait_name,
    output_file = file.path(plots_output_dir, paste0(trait_name, "_ranked_bh_diagnostic.png"))
  )

  make_ranked_top_plots(
    valid_df = valid_df,
    trait_name = trait_name,
    output_dir = ranked_top_dir,
    rank_limits = top_rank_limits
  )

  make_ranked_logx_plot(
    valid_df = valid_df,
    trait_name = trait_name,
    output_file = file.path(ranked_logx_dir, paste0(trait_name, "_ranked_log10x.png"))
  )

  make_pvalue_zoom_plots(
    valid_df = valid_df,
    trait_name = trait_name,
    output_dir = pvalue_zooms_dir,
    zoom_limits = pvalue_zoom_limits
  )

  make_qvalue_zoom_plots(
    valid_df = valid_df,
    trait_name = trait_name,
    output_dir = qvalue_zooms_dir,
    zoom_limits = qvalue_zoom_limits
  )

  make_qq_plot(
    valid_df = valid_df,
    trait_name = trait_name,
    output_file = file.path(qqplots_dir, paste0(trait_name, "_qqplot.png"))
  )

  write_trait_interpretation(
    trait_name = trait_name,
    stats = trait_stats,
    output_file = file.path(summaries_dir, paste0(trait_name, "_interpretation.csv"))
  )

  message(
    "Trait summary: min p = ", signif(trait_stats$min_raw_p_value, 4),
    "; min q = ", signif(trait_stats$min_bh_q_value, 4),
    "; q<=0.05 hits = ", trait_stats$n_q_le_0_05,
    "; q<=0.10 hits = ", trait_stats$n_q_le_0_10,
    "; Bonferroni 0.05 hits = ", trait_stats$n_bonferroni_0_05
  )

  data.frame(
    trait_name = trait_name,
    total_snps = total_snps,
    valid_snps = valid_snps,
    min_raw_p_value = trait_stats$min_raw_p_value,
    max_minus_log10_p = trait_stats$max_minus_log10_p,
    min_bh_q_value = trait_stats$min_bh_q_value,
    n_q_le_0_05 = trait_stats$n_q_le_0_05,
    n_q_le_0_10 = trait_stats$n_q_le_0_10,
    n_bonferroni = trait_stats$n_bonferroni_0_05,
    n_bonferroni_0_10 = trait_stats$n_bonferroni_0_10,
    bonferroni_0_05_cutoff_raw_p = trait_stats$bonferroni_0_05_cutoff_raw_p,
    bonferroni_0_10_cutoff_raw_p = trait_stats$bonferroni_0_10_cutoff_raw_p,
    bh_q_0_05_cutoff_raw_p = trait_stats$bh_q_0_05_cutoff_raw_p,
    bh_q_0_05_cutoff_minus_log10_p = trait_stats$bh_q_0_05_cutoff_minus_log10_p,
    bh_q_0_10_cutoff_raw_p = trait_stats$bh_q_0_10_cutoff_raw_p,
    bh_q_0_10_cutoff_minus_log10_p = trait_stats$bh_q_0_10_cutoff_minus_log10_p,
    lambda_gc = trait_stats$lambda_gc,
    stringsAsFactors = FALSE
  )
}

summary_rows <- lapply(
  gwas_files,
  process_gwas_file,
  snp_map_df = snp_map,
  top_snps_output_dir = top_snps_dir,
  plots_output_dir = plots_dir
)

summary_rows <- summary_rows[!vapply(summary_rows, is.null, logical(1))]

if (!length(summary_rows)) {
  stop("No traits were processed successfully. Check warnings above.")
}

summary_df <- do.call(rbind, summary_rows)
summary_file <- file.path(diagnostics_dir, "combined_fdr_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

message("\nSaved combined summary: ", summary_file)
message("Saved top SNP tables in: ", top_snps_dir)
message("Saved diagnostic plots in: ", plots_dir)
message("Saved per-trait summaries in: ", summaries_dir)
message("\nFDR diagnostics finished successfully.")
