# ===================================================================
# Manhattan Plot Script for V3 (v1.1 - Corrected Naming for V3)
# - This version is specifically adapted for the V3 naming convention.
# - It KEEPS the '_corrected_V3' suffix in the phenotype name.
# - This ensures that newly generated plots match the existing ones
#   and that the skip logic works correctly.
# ===================================================================

# --- Part 1: Setup ---
library(dplyr)
library(qqman)
library(stringr)

# --- Define separate base directories for input and output ---
results_base_dir <- "../data"
plots_base_dir   <- "../GWAS_Plots"

cat("Reading GWAS results from:", normalizePath(results_base_dir), "\n")
cat("Saving output plots to:", normalizePath(plots_base_dir), "\n\n")


# --- Part 2: Automated Plotting Loop ---

# Find all the 'results_...' sub-directories inside the 'data' folder
results_directories <- list.dirs(path = results_base_dir, full.names = TRUE, recursive = FALSE)
results_directories <- results_directories[grepl("results_maf", basename(results_directories))]

if (length(results_directories) == 0) {
  stop("ERROR: No 'results_maf...' directories found in ../data.")
}

cat("Found", length(results_directories), "run directories to process.\n\n")

# --- OUTER LOOP: Iterate over each GWAS RESULTS directory ---
for (current_results_dir in results_directories) {
  
  cat("====================================================\n")
  cat("Processing source directory:", basename(current_results_dir), "\n")
  
  # Define and create the corresponding OUTPUT directory for plots
  plot_dir_name <- sub("results_", "", basename(current_results_dir)) %>%
                   sub("maf", "MAF", .) %>%
                   sub("geno", "GENO", .) %>%
                   sub("thin", "THIN", .) %>%
                   sub("pc", "PC", .)
  
  current_plots_dir <- file.path(plots_base_dir, plot_dir_name)
  
  if (!dir.exists(current_plots_dir)) {
    cat("--> Creating output directory:", basename(current_plots_dir), "\n")
    dir.create(current_plots_dir, recursive = TRUE)
  }

  bim_file <- file.path(current_results_dir, "snp_map.bim")
  if (!file.exists(bim_file)) {
    cat("WARNING: snp_map.bim not found in", basename(current_results_dir), ". Skipping.\n\n")
    next
  }
  
  snp_map <- read.table(bim_file, header = FALSE, col.names = c("CHR", "SNP_ID", "CM", "BP", "A1", "A2"))
  snp_count <- nrow(snp_map)
  bonferroni_p_value <- 0.05 / snp_count
  bonferroni_log10 <- -log10(bonferroni_p_value)
  
  ps_files <- list.files(path = current_results_dir, pattern = "\\.ps$", full.names = TRUE)
  
  # --- INNER LOOP: Iterate over each phenotype file ---
  for (ps_file in ps_files) {
    
    # ===================================================================
    # --- THIS IS THE CRITICAL FIX for V3 NAMING ---
    # Parse the name to keep the '_corrected_V3' part.
    phenotype_name_full <- basename(ps_file) %>% 
                           str_remove("_gwas\\.ps$") %>% 
                           str_remove("^morexV3_")
    # ===================================================================

    cat("\n--- Processing phenotype:", phenotype_name_full, "---\n")
    
    # Extract PC number for the filename
    pc_string <- str_extract(basename(current_plots_dir), "PC[0-9]+") %>% str_to_lower()
    if (is.na(pc_string)) { pc_string <- "pc_unknown" }
    
    # Build the correct, full output filename
    output_filename <- file.path(current_plots_dir, paste0("manhattan_", phenotype_name_full, "_", pc_string, ".png"))
    
    # Robust, file-by-file skip logic will now work correctly
    if (file.exists(output_filename)) {
        cat("Plot", basename(output_filename), "already exists. Skipping.\n")
        next
    }
    
    gwas_results <- read.table(
        ps_file, 
        header = FALSE, 
        col.names = c("SNP_ID_DUMMY", "BETA", "SE", "P"),
        stringsAsFactors = FALSE
    )

    if (nrow(snp_map) != nrow(gwas_results)) {
        cat("WARNING: Row count mismatch for", phenotype_name_full, ". Skipping.\n")
        next
    }
    
    gwas_data_for_plotting <- data.frame(snp_map, P = gwas_results$P)

    gwas_data_for_plotting_cleaned <- gwas_data_for_plotting %>%
      mutate(CHR = as.character(CHR)) %>%
      mutate(CHR_NUM = as.numeric(gsub("[^0-9]", "", CHR)))
      
    run_params_string <- str_replace_all(basename(current_plots_dir), "_", " ")
    pheno_title_name <- str_to_title(str_replace_all(phenotype_name_full, "_", " "))
    plot_title <- paste("Manhattan Plot for", pheno_title_name, "\n(", run_params_string, ")")
    
    max_y_value <- max(na.omit(-log10(gwas_data_for_plotting_cleaned$P)))
    plot_y_limit <- ceiling(max(max_y_value, bonferroni_log10)) + 1
    
    cat("Saving plot to:", output_filename, "\n")
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
    
    annotation_text <- paste(
      paste("SNPs =", format(snp_count, big.mark = ",")),
      paste("Bonf. threshold (-log10p) =", round(bonferroni_log10, 2)),
      paste("Bonf. p-value =", format(bonferroni_p_value, scientific = TRUE, digits = 3)),
      sep = "\n"
    )
    mtext(text = annotation_text, side = 4, line = 0.5, cex = 0.7, adj = 0)
    
    dev.off()
  }
   cat("\nFinished processing all phenotypes in this directory.\n\n")
}

cat("====================================================\n")
cat("ALL V3 PLOTTING JOBS FINISHED SUCCESSFULLY!\n")
cat("====================================================\n")

