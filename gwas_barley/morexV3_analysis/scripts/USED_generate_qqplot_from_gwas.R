suppressPackageStartupMessages({
  library(stringr)
  library(data.table)
  library(parallel)
})

results_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_data/results_maf0.000001_geno1_thin0_pc5"
output_dir  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_data/QQplots_maf0.000001_geno1_thin0_pc5"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ps_files <- c(
  file.path(results_dir, "morexV3_protein_corrected_V3_gwas.ps"),
  file.path(results_dir, "morexV3_starch_corrected_V3_gwas.ps"),
  file.path(results_dir, "morexV3_fiber_corrected_V3_gwas.ps"),
  file.path(results_dir, "morexV3_betaglucan_corrected_V3_gwas.ps")
)

clean_trait_name <- function(ps_file) {
  basename(ps_file) |>
    str_remove("_gwas\\.ps$") |>
    str_remove("^morexV3_") |>
    str_remove("_corrected_V3$") |>
    str_replace_all("_", " ") |>
    str_to_title()
}

calc_lambda_gc <- function(pvals) {
  chisq <- qchisq(1 - pvals, df = 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
}

make_qq_plot <- function(ps_file) {
  if (!file.exists(ps_file)) {
    warning("Missing file: ", ps_file)
    return(NULL)
  }

  trait_label <- clean_trait_name(ps_file)
  trait_file_stub <- basename(ps_file) |>
    str_remove("_gwas\\.ps$")

  message("Processing: ", trait_label)

  dt <- fread(
    ps_file,
    header = FALSE,
    select = 4,
    col.names = "P",
    nThread = 12,
    showProgress = TRUE
  )

  pvals <- dt[["P"]]
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]

  if (!length(pvals)) {
    warning("No valid p-values found in: ", ps_file)
    return(NULL)
  }

  pvals <- sort(pvals)
  observed <- -log10(pvals)
  expected <- -log10(ppoints(length(pvals)))
  lambda_gc <- calc_lambda_gc(pvals)

  max_axis <- ceiling(max(c(observed, expected), na.rm = TRUE) * 1.05)
  output_file <- file.path(output_dir, paste0("qqplot_", trait_file_stub, ".png"))

  png(output_file, width = 8, height = 8, units = "in", res = 300)
  par(mar = c(5, 5, 4, 2) + 0.1)

  plot(
    expected,
    observed,
    pch = 20,
    cex = 0.35,
    col = "navy",
    xlab = "Expected -log10(P)",
    ylab = "Observed -log10(P)",
    main = paste("QQ Plot for", trait_label),
    xlim = c(0, max_axis),
    ylim = c(0, max_axis)
  )

  abline(0, 1, col = "red3", lty = 2, lwd = 2)

  legend(
    "topleft",
    legend = c(
      paste("Valid SNPs:", format(length(pvals), big.mark = ",")),
      paste("Lambda GC:", format(round(lambda_gc, 3), nsmall = 3))
    ),
    bty = "n",
    cex = 0.9
  )

  dev.off()
  message("Saved: ", output_file)
}

mclapply(ps_files, make_qq_plot, mc.cores = 4)