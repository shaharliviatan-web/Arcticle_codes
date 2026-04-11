# ===================================================================
# Manhattan Plot Script (V5 - Clean Titles, Scientific FDR)
# ===================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(qqman)
  library(stringr)
})

# --- Part 1: Setup ---
results_base_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_data"
plots_base_dir   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_Plots"

if (!dir.exists(plots_base_dir)) dir.create(plots_base_dir, recursive = TRUE)

target_results_folder <- "results_maf0.000001_geno1_thin0_pc5"
current_results_dir <- file.path(results_base_dir, target_results_folder)

if (!dir.exists(current_results_dir)) stop("Target directory not found.")

plot_sub_dir <- target_results_folder %>%
                sub("results_", "", .) %>%
                sub("maf", "MAF", .) %>%
                sub("geno", "GENO", .) %>%
                sub("thin", "THIN", .) %>%
                sub("pc", "PC", .)

current_plots_dir <- file.path(plots_base_dir, plot_sub_dir)
if (!dir.exists(current_plots_dir)) dir.create(current_plots_dir, recursive = TRUE)

# --- Part 2: Load Map ---
bim_file <- file.path(current_results_dir, "snp_map.bim")
if (!file.exists(bim_file)) stop("snp_map.bim not found.")

snp_map <- read.table(bim_file, header = FALSE, col.names = c("CHR", "SNP_ID", "CM", "BP", "A1", "A2"))
snp_count <- nrow(snp_map)
bonferroni_log10 <- -log10(0.05 / snp_count)

ps_files <- list.files(path = current_results_dir, pattern = "\\.ps$", full.names = TRUE)
fdr_levels <- c(0.05, 0.10)

# --- Part 3: Main Loop ---
for (ps_file in ps_files) {
  
  phenotype_name_full <- basename(ps_file) %>% str_remove("_gwas\\.ps$") %>% str_remove("^morexV3_")
  cat("\n--- Processing:", phenotype_name_full, "---\n")
  
  gwas_results <- read.table(ps_file, header = FALSE, col.names = c("SNP_ID_DUMMY", "BETA", "SE", "P"))
  
  # Prepare Plot Data
  gwas_data_for_plotting <- data.frame(snp_map, P = gwas_results$P)
  gwas_data_for_plotting_cleaned <- gwas_data_for_plotting %>%
    mutate(CHR_NUM = as.numeric(gsub("[^0-9]", "", as.character(CHR))))
  
  pheno_display_name <- str_replace_all(phenotype_name_full, "_corrected_V3", "") %>%
                        str_replace_all("_", " ") %>% str_to_title()
  
  max_observed <- max(na.omit(-log10(gwas_data_for_plotting_cleaned$P)))

  # --- INNER LOOP (FDR Levels) ---
  for (q_val in fdr_levels) {
      
      # Filename keeps the FDR info so you don't overwrite files
      output_filename <- file.path(current_plots_dir, paste0("manhattan_", phenotype_name_full, "_pc5_FDR", q_val, ".png"))
      
      # ---------------------------------------------------------
      # SCIENTIFIC FDR CALCULATION
      # ---------------------------------------------------------
      p_values <- gwas_results$P
      q_values <- p.adjust(p_values, method = "BH")
      
      sig_indices <- which(q_values <= q_val)
      
      is_theoretical <- FALSE
      
      if (length(sig_indices) > 0) {
        # CASE A: Real hits
        fdr_threshold_p <- max(p_values[sig_indices], na.rm = TRUE)
        fdr_log10 <- -log10(fdr_threshold_p)
        label_suffix <- ""
      } else {
        # CASE B: No hits (Theoretical Line)
        m <- length(na.omit(p_values))
        fdr_threshold_p <- q_val / m
        fdr_log10 <- -log10(fdr_threshold_p)
        is_theoretical <- TRUE
        label_suffix <- " (Theoretical)"
      }
      
      # Determine Y-limit
      plot_y_limit <- ceiling(max(max_observed, bonferroni_log10, 6, fdr_log10, na.rm = TRUE)) + 1
      
      # TITLE RESTORED TO ORIGINAL FORMAT (No FDR mentioned)
      plot_title <- paste("Manhattan Plot for", pheno_display_name)
      
      png(output_filename, width = 12, height = 7, units = 'in', res = 300)
      
      manhattan(
        gwas_data_for_plotting_cleaned,
        main = plot_title,
        chr = "CHR_NUM", bp = "BP", p = "P", snp = "SNP_ID",
        col = c("blue4", "orange3"), 
        suggestiveline = FALSE,
        genomewideline = bonferroni_log10, 
        ylim = c(0, plot_y_limit)
      )
      
      # Red fixed threshold
      abline(h = 6, col = "red", lty = 2, lwd = 1.5)
      
      # Blue FDR Line (Real or Theoretical)
      line_type <- if(is_theoretical) 3 else 2
      abline(h = fdr_log10, col = "blue2", lty = line_type, lwd = 1.5)
      
      annotation_text <- paste(
        paste("SNPs =", format(snp_count, big.mark = ",")),
        paste("Bonferroni =", round(bonferroni_log10, 2)),
        paste0("FDR ", q_val, label_suffix, " = ", round(fdr_log10, 2)),
        paste("Threshold = 6.00"),
        sep = "\n"
      )
      mtext(text = annotation_text, side = 4, line = 0.5, cex = 0.7, adj = 0)
      
      dev.off()
  }
}

cat("\nSUCCESS: V5 Script Finished (Clean Titles).\n")