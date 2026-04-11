# --- 1. Load Libraries ---
library(data.table)
library(ggplot2)

# --- 2. Configuration ---
# Read the SMALL, pre-calculated summary file (this is very fast)
SUMMARY_FILE <- "genome_ld_decay_profile.txt"

# Define the new R^2 threshold you want to plot
NEW_CRITICAL_R2 <- 0.2

# Define the new output filename
NEW_PNG_FILE <- "ld_decay_profile_R0.2_line.png"

# --- 3. Read Summary Data ---
cat("Reading existing summary file:", SUMMARY_FILE, "\n")
ld_decay_profile <- fread(SUMMARY_FILE)

# --- 4. Plot the LD Decay Profile (with new line) ---
cat("Generating new plot with R^2 line at", NEW_CRITICAL_R2, "\n")

# Re-create the plot using the loaded data
ld_plot <- ggplot(ld_decay_profile, aes(x = Distance_kb, y = Avg_R2)) +
    # 1. Add the smoothed LOESS curve
    geom_smooth(method = "loess", se = FALSE, color = "darkred", size = 1.2) +
    # 2. Add points for the binned averages
    geom_point(color = "gray30", alpha = 0.5) +
    
    # --- This is the updated line ---
    # 3. Add the new critical R2 threshold line at 0.2
    geom_hline(yintercept = NEW_CRITICAL_R2, linetype = "dashed", color = "blue", size = 0.8) +
    
    labs(
        title = "LD Decay Profile (Wild Barley, R^2 = 0.2 threshold)",
        x = "Physical Distance (kb)",
        y = expression(Average ~ R^2)
    ) +
    # Use max distance from the actual data for the x-axis limit
    xlim(0, max(ld_decay_profile$Distance_kb, na.rm = TRUE)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# --- 5. Save the new plot ---
ggsave(NEW_PNG_FILE, plot = ld_plot, width = 10, height = 6, units = "in", dpi = 300)

cat("--- Success! New plot saved to:", NEW_PNG_FILE, "---\n")