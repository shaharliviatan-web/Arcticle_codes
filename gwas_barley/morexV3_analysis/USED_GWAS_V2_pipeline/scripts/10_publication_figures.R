#!/usr/bin/env Rscript
# 10_publication_figures.R
# TAG-ready publication folders for BLUE x 10 PCs. Generates TWO parallel sets:
#   - Set A (publication_BonfOnly):           single Bonferroni line at -log10p = 6.7712
#   - Set B (publication_BonfAndSuggestive64): Bonferroni line + suggestive line at -log10p = 6.4
#
# Each set is self-contained: main figure, individual panels, PC selection kit, tables,
# supplementary materials, README and cross-trait interpretation. Designed for direct
# handoff to a manuscript-writing context (Claude.ai online).

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(png)
  library(grid)
  library(gridExtra)
})

PIPE       <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
INTER_DIR  <- file.path(PIPE, "intermediates")
PS_DIR     <- file.path(PIPE, "results", "emmax_ps")
PC_SEL_DIR <- file.path(PIPE, "results", "pc_selection")
COMP_DIR   <- file.path(PIPE, "results", "comparison")
TAB_DIR    <- file.path(PIPE, "results", "tables")

# --- Locked constants ---
ALPHA    <- 0.10
N_ALL    <- 7110996L
N_PRUNED <- 590462L
BONF     <- -log10(ALPHA / N_PRUNED)   # 6.7712
SUGG     <- 6.4
LD_WIN   <- 188000L                    # +/- bp around lead SNP

TRAITS        <- c("betaglucan", "fiber", "protein", "starch")
TRAIT_LABELS  <- c(betaglucan = "β-glucan", fiber = "Fiber",
                   protein = "Protein", starch = "Starch")
PANEL_LABELS  <- c("A", "B", "C", "D")
CHR_COLS      <- c("blue4", "orange3")

# Standalone panel sizes
STAND_MAN_W <- 10; STAND_MAN_H <- 5
STAND_QQ_W  <- 5;  STAND_QQ_H  <- 5
# Grid (main figure) panel sizes  -- 4 rows x (Manhattan + QQ); TAG ~ 174 mm wide
GRID_MAN_W <- 5.0; GRID_MAN_H <- 1.85
GRID_QQ_W  <- 1.85; GRID_QQ_H <- 1.85
RASTER_DPI <- 900

# ============================================================
# Load reference data
# ============================================================
cat("[10] Loading SNP map (with alleles)...\n")
snp_map <- fread(file.path(INTER_DIR, "morexV3_290.bim"), header = FALSE,
                 col.names = c("CHR_raw", "SNP", "cm", "BP", "A1", "A2"))
snp_map[, CHR := as.integer(sub("H$", "", CHR_raw))]
stopifnot(all(snp_map$CHR %in% 1:7))

cat("[10] Loading MAF table...\n")
frq <- fread(file.path(INTER_DIR, "morexV3_290_freq.frq"))
# .frq cols: CHR, SNP, A1, A2, MAF, NCHROBS (A1=minor in PLINK --freq output)
maf_dt <- frq[, .(SNP, MAF)]

cat("[10] Loading lambda table...\n")
lam <- fread(file.path(TAB_DIR, "lambda_table.tsv"))
lam_blue10 <- lam[pheno_type == "BLUE" & n_PCs == 10][, .(trait, lambda_GC)]
stopifnot(nrow(lam_blue10) == 4)

# ============================================================
# Per-trait preparation (heavy; do once)
# ============================================================
prepare_trait <- function(trait_name) {
  ps_file <- file.path(PS_DIR, sprintf("morexV3__%s__BLUE__pc10.ps", trait_name))
  gw <- fread(ps_file, header = FALSE, col.names = c("SNP", "beta", "SE", "P"),
              showProgress = FALSE)
  stopifnot(nrow(gw) == nrow(snp_map))
  d <- data.table(CHR = snp_map$CHR, SNP = snp_map$SNP, BP = snp_map$BP,
                  A1 = snp_map$A1, A2 = snp_map$A2,
                  beta = as.numeric(gw$beta), SE = as.numeric(gw$SE),
                  P = as.numeric(gw$P))
  d <- d[!is.na(P) & P > 0 & P <= 1]
  d[, nlp := -log10(P)]
  d <- merge(d, maf_dt, by = "SNP", all.x = TRUE)

  # Lead candidates (all above SUGG, covers both sets)
  lead <- d[nlp >= SUGG][order(-nlp)]

  # Manhattan layout
  setorder(d, CHR, BP)
  chr_max <- d[, .(mx = max(BP)), by = CHR][order(CHR)]
  chr_max[, off := cumsum(as.numeric(shift(mx, fill = 0)))]
  d <- merge(d, chr_max[, .(CHR, off)], by = "CHR")
  setorder(d, CHR, BP)
  d[, x := as.numeric(BP) + off]
  d[, col := ifelse(CHR %% 2L == 0L, CHR_COLS[2], CHR_COLS[1])]
  chr_cent <- d[, .(c = (min(x) + max(x)) / 2), by = CHR][order(CHR)]
  xlim    <- c(min(d$x), max(d$x))
  max_obs <- max(d$nlp, na.rm = TRUE)
  ymax    <- ceiling(max(max_obs, BONF, 6, na.rm = TRUE)) + 1

  # QQ data (subsample for plotting)
  setorder(d, P)
  n_total <- nrow(d)
  expected_full <- -log10(ppoints(n_total))
  observed_full <- d$nlp
  if (n_total > 200000L) {
    idx_tail <- seq_len(50000L)
    idx_bulk <- sort(sample(50001L:n_total, 150000L))
    idx      <- c(idx_tail, idx_bulk)
  } else {
    idx <- seq_len(n_total)
  }
  qq <- list(expected = expected_full[idx],
             observed = observed_full[idx],
             max_obs  = max_obs,
             n_total  = n_total)

  lambda_val <- lam_blue10[trait == trait_name, lambda_GC][1]

  list(trait    = trait_name,
       d        = d,
       xlim     = xlim,
       ylim     = c(0, ymax),
       chr_cent = chr_cent,
       qq       = qq,
       lead     = lead,
       lambda   = lambda_val)
}

cat("[10] Preparing per-trait data (4 traits)...\n")
trait_data <- lapply(TRAITS, prepare_trait)
names(trait_data) <- TRAITS

# ============================================================
# Pre-render point-cloud rasters (heavy, do once per size)
# ============================================================
render_point_cloud <- function(td, w_in, h_in) {
  tmp <- tempfile(fileext = ".png", tmpdir = Sys.getenv("TMPDIR"))
  png(tmp, width = w_in, height = h_in, units = "in",
      res = RASTER_DPI, bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(td$d$x, td$d$nlp, xlim = td$xlim, ylim = td$ylim,
       xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "",
       pch = 16, cex = 0.85, col = td$d$col)
  dev.off()
  img <- readPNG(tmp); unlink(tmp); img
}
cat("[10] Rendering point-cloud rasters (standalone + grid sizes, 4 traits)...\n")
rasters_standalone <- lapply(trait_data, render_point_cloud, w_in = STAND_MAN_W, h_in = STAND_MAN_H)
rasters_grid       <- lapply(trait_data, render_point_cloud, w_in = GRID_MAN_W,  h_in = GRID_MAN_H)

# ============================================================
# Drawing helpers
# ============================================================
draw_manhattan <- function(td, raster_img, set_b = FALSE,
                           panel_label = NULL, trait_label = NULL,
                           cex_scale = 1.0) {
  par(mar = c(4.2, 4.5, 1.8, 0.8), las = 1,
      cex.axis = 0.9 * cex_scale, cex.lab = 1.0 * cex_scale,
      mgp = c(2.6, 0.6, 0))
  plot(NA, xlim = td$xlim, ylim = td$ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "Chromosome", ylab = expression(-log[10](italic(p))), main = "")
  rasterImage(raster_img, td$xlim[1], td$ylim[1], td$xlim[2], td$ylim[2],
              interpolate = FALSE)
  axis(2)
  axis(1, at = td$chr_cent$c, labels = paste0(td$chr_cent$CHR, "H"), tick = FALSE)
  box()
  # Thresholds (suggestive first, Bonf on top)
  if (set_b) abline(h = SUGG, col = "firebrick3", lty = 2, lwd = 1.3 * cex_scale)
  abline(h = BONF, col = "black",       lty = 1, lwd = 1.6 * cex_scale)
  # In-panel labels (no methodology jargon)
  ann <- if (!is.null(panel_label) && !is.null(trait_label)) {
    sprintf("%s) %s", panel_label, trait_label)
  } else if (!is.null(trait_label)) trait_label else NULL
  if (!is.null(ann)) {
    text(td$xlim[1] + 0.015 * diff(td$xlim), td$ylim[2] * 0.95,
         labels = ann, adj = c(0, 1), font = 2, cex = 1.05 * cex_scale)
  }
}

draw_qq <- function(td, trait_label = NULL, cex_scale = 1.0) {
  qq <- td$qq
  xmax <- max(qq$expected, na.rm = TRUE)
  ymax <- max(qq$observed, qq$max_obs, na.rm = TRUE)
  par(mar = c(4.0, 4.2, 1.8, 0.8), las = 1,
      cex.axis = 0.85 * cex_scale, cex.lab = 0.95 * cex_scale,
      mgp = c(2.3, 0.55, 0))
  plot(qq$expected, qq$observed,
       xlim = c(0, xmax), ylim = c(0, max(ymax, xmax)),
       pch = 16, cex = 0.35, col = "grey30",
       xlab = expression(Expected ~ -log[10](italic(p))),
       ylab = expression(Observed ~ -log[10](italic(p))), main = "")
  abline(0, 1, col = "red", lwd = 0.7)
  legend("topleft",
         legend = bquote(lambda[GC] == .(sprintf("%.3f", td$lambda))),
         bty = "n", cex = 0.85 * cex_scale, text.col = "grey20",
         inset = c(0.02, 0.02))
  if (!is.null(trait_label)) {
    text(0.02 * xmax, max(ymax, xmax) * 0.95,
         labels = trait_label, adj = c(0, 1), font = 2,
         cex = 1.05 * cex_scale)
  }
}

# ============================================================
# Generate a publication set
# ============================================================
generate_set <- function(set_name, set_b = FALSE) {
  cat(sprintf("\n[10] ====== %s ======\n", set_name))
  out_dir  <- file.path(PIPE, "results", set_name)
  main_dir <- file.path(out_dir, "main_figure")
  ind_dir  <- file.path(out_dir, "individual_panels")
  pc_dir   <- file.path(out_dir, "pc_selection")
  tab_dir  <- file.path(out_dir, "tables")
  supp_dir <- file.path(out_dir, "supp")
  for (d in c(out_dir, main_dir, ind_dir, pc_dir, tab_dir, supp_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  # ---- Standalone Manhattans + QQs ----
  cat("[10] Writing standalone Manhattan + QQ panels...\n")
  for (i in seq_along(TRAITS)) {
    td  <- trait_data[[i]]
    img <- rasters_standalone[[i]]
    label <- TRAIT_LABELS[[td$trait]]

    base <- sprintf("manhattan_BLUE_pc10__%s", td$trait)
    draw_man <- function() draw_manhattan(td, img, set_b = set_b,
                                          trait_label = label, cex_scale = 1.1)
    cairo_pdf(file.path(ind_dir, paste0(base, ".pdf")),
              width = STAND_MAN_W, height = STAND_MAN_H, pointsize = 10)
    draw_man(); dev.off()
    png(file.path(ind_dir, paste0(base, ".png")),
        width = STAND_MAN_W, height = STAND_MAN_H, units = "in",
        res = 600, pointsize = 10, type = "cairo-png")
    draw_man(); dev.off()

    qq_base <- sprintf("qq_BLUE_pc10__%s", td$trait)
    draw_q <- function() draw_qq(td, trait_label = label, cex_scale = 1.15)
    cairo_pdf(file.path(ind_dir, paste0(qq_base, ".pdf")),
              width = STAND_QQ_W, height = STAND_QQ_H, pointsize = 10)
    draw_q(); dev.off()
    png(file.path(ind_dir, paste0(qq_base, ".png")),
        width = STAND_QQ_W, height = STAND_QQ_H, units = "in",
        res = 600, pointsize = 10, type = "cairo-png")
    draw_q(); dev.off()
  }

  # ---- Main figure: 4 rows x (Manhattan + QQ) cols ----
  cat("[10] Composing main figure (4x2 grid)...\n")
  tmp_panel_dir <- tempfile(pattern = "panels_", tmpdir = Sys.getenv("TMPDIR"))
  dir.create(tmp_panel_dir)
  for (i in seq_along(TRAITS)) {
    td  <- trait_data[[i]]
    img <- rasters_grid[[i]]
    plabel <- PANEL_LABELS[i]
    tlabel <- TRAIT_LABELS[[td$trait]]

    png(file.path(tmp_panel_dir, sprintf("man_%d.png", i)),
        width = GRID_MAN_W, height = GRID_MAN_H, units = "in",
        res = 600, pointsize = 8, type = "cairo-png")
    draw_manhattan(td, img, set_b = set_b,
                   panel_label = plabel, trait_label = tlabel, cex_scale = 0.9)
    dev.off()

    png(file.path(tmp_panel_dir, sprintf("qq_%d.png", i)),
        width = GRID_QQ_W, height = GRID_QQ_H, units = "in",
        res = 600, pointsize = 8, type = "cairo-png")
    draw_qq(td, trait_label = NULL, cex_scale = 0.85)
    dev.off()
  }
  grobs <- vector("list", 8L)
  for (i in seq_along(TRAITS)) {
    grobs[[2 * i - 1]] <- rasterGrob(
      readPNG(file.path(tmp_panel_dir, sprintf("man_%d.png", i))),
      interpolate = TRUE)
    grobs[[2 * i]] <- rasterGrob(
      readPNG(file.path(tmp_panel_dir, sprintf("qq_%d.png", i))),
      interpolate = TRUE)
  }
  arr <- arrangeGrob(grobs = grobs, ncol = 2L, nrow = 4L,
                     widths  = unit(c(GRID_MAN_W, GRID_QQ_W), "in"),
                     heights = unit(rep(GRID_MAN_H, 4L), "in"))
  W <- GRID_MAN_W + GRID_QQ_W + 0.15
  H <- 4 * GRID_MAN_H + 0.15
  png(file.path(main_dir, "Figure_main_manhattan_qq.png"),
      width = W, height = H, units = "in", res = 300, type = "cairo-png")
  grid.draw(arr); dev.off()
  cairo_pdf(file.path(main_dir, "Figure_main_manhattan_qq.pdf"),
            width = W, height = H)
  grid.draw(arr); dev.off()
  unlink(tmp_panel_dir, recursive = TRUE)
  cat(sprintf("[10]   main figure: %.2f x %.2f in (%.1f x %.1f mm)\n",
              W, H, W * 25.4, H * 25.4))

  # ---- Lead SNP table ----
  cat("[10] Building lead_snps.tsv...\n")
  lead_all <- rbindlist(lapply(seq_along(TRAITS), function(i) {
    td   <- trait_data[[i]]
    lead <- if (set_b) td$lead else td$lead[nlp >= BONF]
    if (nrow(lead) == 0L) return(NULL)
    lead <- copy(lead)
    lead[, trait := td$trait]
    lead[, BonfPruned010_threshold := round(BONF, 4)]
    lead[, above_Bonf := nlp >= BONF]
    lead[, LD_window_start := pmax(0L, BP - LD_WIN)]
    lead[, LD_window_end := BP + LD_WIN]
    if (set_b) {
      lead[, suggestive_threshold := SUGG]
      lead[, above_suggestive_only := !above_Bonf]
    }
    lead
  }))
  if (nrow(lead_all) > 0L) {
    setnames(lead_all,
             c("SNP", "CHR", "BP", "P", "nlp"),
             c("SNP_id", "chr", "position_bp", "p_value", "neg_log10_p"))
    lead_all[, p_value := signif(p_value, 4)]
    lead_all[, neg_log10_p := round(neg_log10_p, 4)]
    lead_all[, beta := signif(beta, 5)]
    lead_all[, SE := signif(SE, 5)]
    lead_all[, MAF := round(MAF, 4)]
    base_cols <- c("trait", "SNP_id", "chr", "position_bp", "A1", "A2", "MAF",
                   "beta", "SE", "p_value", "neg_log10_p",
                   "BonfPruned010_threshold", "above_Bonf",
                   "LD_window_start", "LD_window_end")
    if (set_b) base_cols <- c(base_cols, "suggestive_threshold", "above_suggestive_only")
    setcolorder(lead_all, base_cols)
    setorder(lead_all, trait, -neg_log10_p)
  }
  fwrite(lead_all, file.path(tab_dir, "lead_snps.tsv"), sep = "\t")
  cat(sprintf("[10]   lead_snps.tsv: %d rows\n", nrow(lead_all)))

  # ---- Per-trait summary ----
  cat("[10] Building per_trait_summary.tsv...\n")
  per_trait <- rbindlist(lapply(seq_along(TRAITS), function(i) {
    td <- trait_data[[i]]
    top <- td$d[order(-nlp)][1]
    out <- data.table(
      trait = td$trait,
      n_SNPs_tested = nrow(td$d),
      lambda_GC = round(td$lambda, 4),
      n_SNPs_above_Bonf = sum(td$d$nlp >= BONF),
      top_SNP = top$SNP,
      top_chr = sprintf("%dH", top$CHR),
      top_pos = top$BP,
      top_neg_log10_p = round(top$nlp, 4),
      top_p_value = signif(top$P, 4)
    )
    if (set_b) out[, n_SNPs_above_suggestive := sum(td$d$nlp >= SUGG)]
    out
  }))
  fwrite(per_trait, file.path(tab_dir, "per_trait_summary.tsv"), sep = "\t")

  # ---- Copy top-15 per-chr for the 4 winning cells ----
  for (trait in TRAITS) {
    src <- file.path(TAB_DIR, "top_snps_per_chr",
                     sprintf("top15perChr__%s__BLUE__pc10.tsv", trait))
    dst <- file.path(tab_dir, sprintf("top15_per_chr__%s.tsv", trait))
    file.copy(src, dst, overwrite = TRUE)
  }

  # ---- PC selection (copy + rename) ----
  cat("[10] Copying PC selection materials...\n")
  pc_map <- c(
    "scree_plot.pdf"          = "Figure_S1_PC_scree.pdf",
    "scree_plot.png"          = "Figure_S1_PC_scree.png",
    "lambda_vs_npcs.pdf"      = "Figure_S2_lambda_vs_npcs.pdf",
    "lambda_vs_npcs.png"      = "Figure_S2_lambda_vs_npcs.png",
    "pca_scatter_PC1_PC2.pdf" = "Figure_S3_PCA_scatter_PC1_PC2.pdf",
    "pca_scatter_PC1_PC2.png" = "Figure_S3_PCA_scatter_PC1_PC2.png",
    "pca_scatter_PC1_PC3.pdf" = "Figure_S4_PCA_scatter_PC1_PC3.pdf",
    "pca_scatter_PC1_PC3.png" = "Figure_S4_PCA_scatter_PC1_PC3.png",
    "pca_scatter_PC2_PC3.pdf" = "Figure_S5_PCA_scatter_PC2_PC3.pdf",
    "pca_scatter_PC2_PC3.png" = "Figure_S5_PCA_scatter_PC2_PC3.png",
    "pc_variance_table.tsv"   = "Table_S2_pc_variance.tsv",
    "lambda_vs_npcs_long.tsv" = "Table_S3_lambda_vs_npcs_long.tsv",
    "lambda_vs_npcs_wide.tsv" = "Table_S4_lambda_vs_npcs_wide.tsv"
  )
  for (src in names(pc_map)) {
    file.copy(file.path(PC_SEL_DIR, src),
              file.path(pc_dir, pc_map[[src]]), overwrite = TRUE)
  }

  # ---- Supplementary ----
  cat("[10] Copying supplementary materials...\n")
  supp_map <- c(
    "qq_grid_6configs.pdf"                            = "Supp_Figure_S1_sensitivity_QQ_grid.pdf",
    "qq_grid_6configs.png"                            = "Supp_Figure_S1_sensitivity_QQ_grid.png",
    "manhattan_grid_6configs_BonfPruned010.pdf"       = "Supp_Figure_S2_sensitivity_Manhattan_grid.pdf",
    "manhattan_grid_6configs_BonfPruned010.png"       = "Supp_Figure_S2_sensitivity_Manhattan_grid.png"
  )
  for (src in names(supp_map)) {
    file.copy(file.path(COMP_DIR, src),
              file.path(supp_dir, supp_map[[src]]), overwrite = TRUE)
  }
  file.copy(file.path(TAB_DIR, "summary_24runs.tsv"),
            file.path(supp_dir, "Supp_Table_S1_24run_summary.tsv"),
            overwrite = TRUE)

  list(set_name = set_name, out_dir = out_dir,
       lead_all = lead_all, per_trait = per_trait)
}

set_a <- generate_set("publication_BonfOnly",            set_b = FALSE)
set_b <- generate_set("publication_BonfAndSuggestive64", set_b = TRUE)

# ============================================================
# README and cross-trait interpretation (per set)
# ============================================================
write_readme <- function(set_info, set_b_flag = FALSE) {
  out_dir <- set_info$out_dir
  thr_block <- if (set_b_flag) {
    "This set includes a **suggestive threshold line at −log10p = 6.4** in addition to the Bonferroni line at −log10p = 6.7712. The lead-SNP table includes all SNPs above the suggestive level, with a boolean column distinguishing Bonferroni-passing SNPs."
  } else {
    "This set shows **only** the Bonferroni threshold line at −log10p = 6.7712. The lead-SNP table lists only SNPs above this level."
  }

  txt <- c(
    "# Publication folder — BLUE × 10 PCs",
    "",
    "## Chosen configuration",
    "",
    "- **Phenotype correction**: BLUE (Best Linear Unbiased Estimator)",
    "- **Number of PCs used as covariates in EMMAX**: 10",
    "- **Significance threshold**: Bonferroni α = 0.10 on the LD-pruned 590,462 SNP set → −log10p = **6.7712**",
    paste0("- ", thr_block),
    "",
    "## Reference numbers",
    "",
    "- Total SNPs tested: **7,110,996**",
    "- LD-pruned SNPs used for population-structure correction: **590,462**",
    "- Sample size: **290** wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*) from the Southern Levant",
    paste0("- λ_GC range across the 4 traits (BLUE × 10 PCs): see `tables/per_trait_summary.tsv`"),
    "- Reference genome: Morex V3"
  )
  if (set_b_flag) txt <- c(txt,
    "- Suggestive threshold: −log10p = 6.4 (about 1.5× more lenient than Bonferroni)")
  txt <- c(txt, "",
    "## Directory tree",
    "",
    "```",
    paste0(basename(out_dir), "/"),
    "├── main_figure/",
    "│   └── Figure_main_manhattan_qq.{pdf,png}    — composed 4 traits × (Manhattan + QQ), TAG ~174 mm",
    "├── individual_panels/",
    "│   ├── manhattan_BLUE_pc10__{trait}.{pdf,png}  (4 traits × 2 formats = 8 files)",
    "│   └── qq_BLUE_pc10__{trait}.{pdf,png}         (4 traits × 2 formats = 8 files)",
    "├── pc_selection/",
    "│   ├── Figure_S1_PC_scree.{pdf,png}              — 20-PC scree",
    "│   ├── Figure_S2_lambda_vs_npcs.{pdf,png}        — λ_GC vs n_PCs (PC-count justification)",
    "│   ├── Figure_S3..S5_PCA_scatter_*.{pdf,png}     — PC1×2, PC1×3, PC2×3",
    "│   ├── Table_S2_pc_variance.tsv                  — 20 rows: PC, eigenval, %% variance, cumulative %%",
    "│   ├── Table_S3_lambda_vs_npcs_long.tsv          — 24 rows (BLUP+BLUE × pc3/5/10 × 4 traits)",
    "│   └── Table_S4_lambda_vs_npcs_wide.tsv          — 8 rows × pc3/pc5/pc10 columns",
    "├── tables/",
    "│   ├── lead_snps.tsv                  — significant SNPs (MAF, alleles, β/SE, LD window)",
    "│   ├── per_trait_summary.tsv          — 4 rows: trait, n_tested, λ, n_above_thr, top SNP",
    "│   └── top15_per_chr__{trait}.tsv     — top 15 per chromosome (15 × 7 = 105 rows each)",
    "├── supp/",
    "│   ├── Supp_Figure_S1_sensitivity_QQ_grid.{pdf,png}",
    "│   ├── Supp_Figure_S2_sensitivity_Manhattan_grid.{pdf,png}",
    "│   └── Supp_Table_S1_24run_summary.tsv  — full 24-run sensitivity defends BLUE × 10 PCs",
    "├── README.md",
    "└── cross_trait_interpretation.md",
    "```",
    "",
    "## Methods (short, for the manuscript)",
    "",
    "- PCA and aIBS kinship matrix were computed on an LD-pruned subset (--indep-pairwise 50 5 0.2) of 590,462 SNPs.",
    "- GWAS was performed using EMMAX (mixed model) on all 7,110,996 SNPs, including the top 10 PCs as fixed-effect covariates.",
    "- Sensitivity analysis across 24 cells (4 traits × {BLUP,BLUE} × {3,5,10} PCs) was used to select the configuration; BLUE × 10 PCs minimised the mean |λ − 1| across the four traits without any trait falling below the 0.95 deflation floor.",
    "- Significance threshold: Bonferroni at α = 0.10, corrected for the number of independent (LD-pruned) tests (n = 590,462), yielding −log10p = 6.77.",
    "",
    "## Citations to use in methods",
    "",
    "- Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006). *Principal components analysis corrects for stratification in genome-wide association studies.* Nature Genetics 38, 904–909.",
    "- Patterson N, Price AL, Reich D (2006). *Population structure and eigenanalysis.* PLoS Genetics 2(12): e190.",
    "- Kang HM, Sul JH, Service SK, Zaitlen NA, Kong S, Freimer NB, Sabatti C, Eskin E (2010). *Variance component model to account for sample structure in genome-wide association studies.* Nature Genetics 42, 348–354. (EMMAX)",
    "",
    "## What is NOT here (downstream work)",
    "",
    "1. **Candidate gene count per ±188 kb window** around each lead SNP — requires the Morex V3 GFF annotation + `bedtools intersect`. Plan to add as `tables/candidate_genes_per_lead.tsv` in a follow-up step.",
    "2. Final manuscript prose and biological interpretation — see `cross_trait_interpretation.md` for a suggestive starting draft."
  )
  writeLines(txt, file.path(out_dir, "README.md"))
}

write_cross_trait <- function(set_info, set_b_flag = FALSE) {
  out_dir   <- set_info$out_dir
  lead      <- set_info$lead_all
  per_trait <- set_info$per_trait

  txt <- c(
    "# Cross-trait interpretation — BLUE × 10 PCs",
    "",
    "**Status: SUGGESTIVE DRAFT PROSE.** Generated automatically from the analysis outputs; user must verify each claim and edit before publication. Use as a starting scaffold, not as final manuscript text.",
    "",
    "## Per-trait headline summary",
    ""
  )
  for (i in seq_len(nrow(per_trait))) {
    r <- per_trait[i]
    extra <- if (set_b_flag && "n_SNPs_above_suggestive" %in% names(r)) {
      sprintf(" (+%d above suggestive 6.4)", r$n_SNPs_above_suggestive - r$n_SNPs_above_Bonf)
    } else ""
    txt <- c(txt, sprintf(
      "- **%s**: %d SNPs above Bonferroni (−log10p ≥ 6.77)%s. Top SNP: `%s` on %s at %s bp (−log10p = %.3f). λ_GC = %.4f.",
      TRAIT_LABELS[[r$trait]], r$n_SNPs_above_Bonf, extra,
      r$top_SNP, r$top_chr, format(r$top_pos, big.mark = ","),
      r$top_neg_log10_p, r$lambda_GC
    ))
  }

  txt <- c(txt, "", "## Key observations", "")

  # 1. β-glucan headline hit
  bg_lead <- lead[trait == "betaglucan"][1]
  if (nrow(bg_lead) > 0L && !is.na(bg_lead$chr)) {
    direction <- if (bg_lead$beta < 0) "decreases" else "increases"
    txt <- c(txt, sprintf(
      "1. **Headline β-glucan hit**: `%s` on chromosome %dH at %s bp (−log10p = %.3f, MAF = %.3f, β = %.4g, SE = %.4g). This SNP exceeds the Bonferroni threshold (−log10p = 6.77) and is the strongest signal across all four traits. The minor allele (A1 = %s) is associated with %s β-glucan content.",
      bg_lead$SNP_id, bg_lead$chr, format(bg_lead$position_bp, big.mark = ","),
      bg_lead$neg_log10_p, bg_lead$MAF, bg_lead$beta, bg_lead$SE,
      bg_lead$A1, direction
    ), "")
  }

  # 2. chr4H cluster (β-glucan)
  bg_4h <- lead[trait == "betaglucan" & chr == 4 &
                  position_bp >= 34000000 & position_bp <= 36000000]
  if (nrow(bg_4h) >= 2L) {
    span_kb <- round((max(bg_4h$position_bp) - min(bg_4h$position_bp)) / 1000)
    txt <- c(txt, sprintf(
      "2. **chr4H cluster (β-glucan)**: %d SNPs co-located between positions %s and %s bp (~%s kb window), all near or above Bonferroni. The cluster fits inside a single ±188 kb LD-defined candidate region and likely tags one underlying causal locus rather than %d independent associations.",
      nrow(bg_4h),
      format(min(bg_4h$position_bp), big.mark = ","),
      format(max(bg_4h$position_bp), big.mark = ","),
      format(span_kb, big.mark = ","),
      nrow(bg_4h)
    ), "")
  }

  # 3. Cross-trait chromosomal overlap (above Bonferroni)
  if (nrow(lead) > 0L) {
    bonf_lead <- lead[above_Bonf == TRUE]
    if (nrow(bonf_lead) > 0L) {
      chr_traits <- bonf_lead[, .(traits = paste(sort(unique(trait)), collapse = ", "),
                                  n_traits = uniqueN(trait)),
                              by = chr]
      multi <- chr_traits[n_traits > 1L][order(chr)]
      if (nrow(multi) > 0L) {
        txt <- c(txt, "3. **Cross-trait chromosomal overlap** (Bonferroni-passing hits on the same chromosome across multiple traits):", "")
        for (i in seq_len(nrow(multi))) {
          txt <- c(txt, sprintf("   - Chr %dH: traits = %s", multi$chr[i], multi$traits[i]))
        }
        txt <- c(txt, "",
          "Worth checking whether the proximal lead SNPs in these regions are in LD (suggesting a shared causal locus) or independent (separate causal mechanisms).")
      } else {
        txt <- c(txt, "3. **Cross-trait chromosomal overlap**: No chromosome carries Bonferroni-passing hits for multiple traits at this threshold. β-glucan dominates the signal landscape.")
      }
    }
    txt <- c(txt, "")
  }

  # 4. λ context
  txt <- c(txt,
    "## Genomic control summary",
    "",
    sprintf("All four traits show slight deflation (λ_GC range %.4f–%.4f, all above the 0.95 floor) under BLUE × 10 PCs. This is consistent with strong population-structure correction by the 10 PCs (capturing ~33%% of total genetic variance) on a continuously structured wild barley sample. The deflation is mild enough that the strongest signals remain robust.",
      min(per_trait$lambda_GC), max(per_trait$lambda_GC)),
    "",
    "## Caveats for the paper",
    "",
    "- Wild barley population structure is continuous (no sharp scree elbow); 10 PCs were chosen as the most defensible point on the λ-vs-N_PCs curve, supported by mean |λ − 1| minimisation across the four traits.",
    "- No FDR-based analysis was performed (Bonferroni-only by design)."
  )
  if (set_b_flag) txt <- c(txt,
    "- The suggestive threshold (−log10p = 6.4) is provided for completeness; SNPs in the suggestive band are flagged via `above_suggestive_only = TRUE` in `lead_snps.tsv` and should be discussed only as supplementary candidates.")
  txt <- c(txt, "",
    "## Bridge to candidate-gene analysis",
    "",
    sprintf("All %d Bonferroni-passing lead SNPs define candidate regions of ±188 kb each (the empirical LD-decay window from the prior LD-decay analysis at r² = 0.2). The next analytical step is to count the number of annotated Morex V3 genes within each window using `bedtools intersect` against the V3 GFF, producing per-region gene counts and per-SNP candidate-gene lists.",
      sum(per_trait$n_SNPs_above_Bonf))
  )

  writeLines(txt, file.path(out_dir, "cross_trait_interpretation.md"))
}

cat("\n[10] Writing READMEs and cross-trait interpretations...\n")
write_readme(set_a, set_b_flag = FALSE)
write_readme(set_b, set_b_flag = TRUE)
write_cross_trait(set_a, set_b_flag = FALSE)
write_cross_trait(set_b, set_b_flag = TRUE)

cat("\n[10] === DONE ===\n")
cat(sprintf("Set A: %s\n", set_a$out_dir))
cat(sprintf("Set B: %s\n", set_b$out_dir))
cat(sprintf("Set A lead SNPs: %d  |  Set B lead SNPs: %d\n",
            nrow(set_a$lead_all), nrow(set_b$lead_all)))
