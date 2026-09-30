# 00_config.R — shared paths, constants and helpers for every R script of this step.
# Sourced first by scripts 02–10. Read-only on the rest of the project: every input
# below is read, and every output goes inside this folder.
#
# Reorganised 2026-09-30 (user request): the 2026-09-25 check was split into ordered
# scripts with one config; the Desert test became two-sided; the haplotype groups of
# the four presented genes were added; figures redrawn to the TAG figure spec.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

.args <- commandArgs(trailingOnly = FALSE)
HERE  <- normalizePath(file.path(dirname(sub("^--file=", "", .args[grep("^--file=", .args)])), ".."))
PROJ  <- normalizePath(file.path(HERE, ".."))
TAB   <- file.path(HERE, "results", "tables")
FIG   <- file.path(HERE, "results", "figures")
INTER <- file.path(HERE, "intermediates")
LOGS  <- file.path(HERE, "logs")
for (d in c(TAB, FIG, INTER, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Inputs (read-only) -----------------------------------------------------------------
PLINK        <- "/usr/local/bin/plink"          # bare `plink` resolves to a broken v0.76
GWAS_DIR     <- file.path(PROJ, "01_USED_GWAS_V2_pipeline")
BFILE        <- file.path(GWAS_DIR, "intermediates/morexV3_290")
FRQ_FILE     <- file.path(GWAS_DIR, "intermediates/morexV3_290_freq.frq")
KIN_FILE     <- file.path(GWAS_DIR, "intermediates/morexV3_kinship.aIBS.kinf")   # rows/cols = .fam order
LOCI_FILE    <- file.path(GWAS_DIR, "results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv")
REGION_FILE  <- file.path(PROJ, "00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables/A4_site_boxplot_data.csv")
PHENO_DIR    <- "/mnt/data/shahar/gwas_barley/data/inputs"                   # <trait>_corrected_V3.pheno = the EMMAX BLUPs
STEP04_CACHE <- file.path(PROJ, "04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Cache")
GDSL_ASSIGN  <- file.path(PROJ, "03_01_7H_branch_Starch_Fiber_shared_signal_explore/results/tables/mgmin3_haplotype_assignment.tsv")

# ---- Constants ----------------------------------------------------------------------------
SEED   <- 20260925
N_PERM <- 10000                                 # site-permutation replicates
REGION_ORDER <- c("North", "Coast", "Desert",
                  "HZ1 (North-Coast)", "HZ2 (North-Desert)", "HZ3 (Coast-Desert)")
REGION_SHORT <- c("North", "Coast", "Desert", "HZ1", "HZ2", "HZ3")
TRAIT_ORDER  <- c("betaglucan", "fiber", "starch", "protein")
TRAIT_LAB    <- c(betaglucan = "β-glucan", fiber = "Fiber", starch = "Starch", protein = "Protein")

# The four genes presented in Results Ch. 3–4, with the crosshap run that defined their
# haplotype groups (Ch. 3: step 04, MGmin 2 / eps 0.6; Ch. 4: 03_01 branch, MGmin 3 / eps 0.9)
GENES <- tibble::tribble(
  ~gene,    ~gene_id,                     ~trait,           ~lead_SNP,       ~source,
  "GPAT6",  "HORVU.MOREX.r3.3HG0301300",  "fiber",          "3H:544084309",  "step04",
  "GH17",   "HORVU.MOREX.r3.5HG0487060",  "fiber",          "5H:462564771",  "step04",
  "PHT4;3", "HORVU.MOREX.r3.3HG0301710",  "starch",         "3H:546433616",  "step04",
  "GDSL",   "HORVU.MOREX.r3.7HG0729030",  "fiber; starch",  "7H:573606306; 7H:573606460", "03_01")
GENE_TAG <- c("3H:544084309" = "GPAT6", "5H:462564771" = "GH17", "3H:546433616" = "PHT4;3",
              "7H:573606306" = "GDSL", "7H:573606460" = "GDSL")

# ---- Figure style (as 08_USED_creating_figures/make_figure_4.R and _5.R) -------------------
FONT <- "Liberation Sans"                       # Arial-metric; Arial is not installed
PT   <- 9                                       # lettering size on the page (TAG: 8–12 pt); 9 pt as Figs. 4–5 (2026-09-30, was 8)
DPI  <- 600; W_MM <- 174; H_MAX <- 234
PAL_REGION <- c("#2166AC", "#1B7837", "#B2182B", "#67A9CF", "#762A83", "#E08214")  # Ch. 1 palette
REGION_TXT <- c("white", "white", "white", "#1A1A1A", "white", "#1A1A1A")
INK <- "#1F2430"; INK2 <- "#5B6170"; GREY <- "#C3C8D0"; GRID <- "#E3E6EA"
# Barcode cells (Fig. S1). Changed 2026-09-30 (user: major and no-call too similar, no-call too close to the
# white background), then swapped the same day (user): major = light tan, no call = gray. OKLab lightness
# minor 0.26 / no call 0.71 / major 0.89 vs white 1.00, so the three states separate in grayscale and for
# color-vision deficiency too, and major also differs in hue. Delta E (OKLab x100): major vs no call 20.5
# (was 12.3), major vs white 13.7, no call vs white 29.4 (was 4.6).
CELL_MINOR <- "#1F2430"; CELL_MAJOR <- "#F0D89E"; CELL_NOCALL <- "#9AA1AB"

# ---- Shared loaders -----------------------------------------------------------------------
load_accessions <- function() {
  a <- read_csv(REGION_FILE, show_col_types = FALSE) %>% distinct(region, site, location, short_Tag)
  stopifnot(nrow(a) == 290, !anyDuplicated(a$short_Tag))
  a %>% mutate(IID = paste0("HS", sub("_", "", short_Tag)),
               region = factor(region, levels = REGION_ORDER))
}
load_pheno <- function() {
  bind_rows(lapply(TRAIT_ORDER, function(t) {
    p <- read.table(file.path(PHENO_DIR, paste0(t, "_corrected_V3.pheno")), col.names = c("FID", "IID", "y"))
    data.frame(trait = t, IID = p$IID, y = p$y)
  }))
}

# ---- Site-permutation null ------------------------------------------------------------------
# Accessions stay in their site; the six region labels are shuffled among the 29 sites
# (North keeps 4 sites, every other region 5). One fixed matrix, shared by every test.
site_perm_matrix <- function(acc) {
  sr <- acc %>% distinct(site, region) %>% arrange(site)
  stopifnot(nrow(sr) == 29, !anyDuplicated(sr$site))
  set.seed(SEED)
  list(sites = sr$site, labels = replicate(N_PERM, sample(as.character(sr$region))))
}
chisq_stat <- function(g, r) {                  # g: group factor, r: region factor
  tab <- table(g, r); e <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  sum((tab - e)^2 / e, na.rm = TRUE)
}
# Two tests for one unit (a lead SNP's allele split, or a haplotype group vs the rest):
#   Desert  — statistic: Desert share of the focal set minus Desert share of the others;
#             two-sided p = 2 x the smaller tail (capped at 1), +1 correction.
#   regions — statistic: chi-square of focal/other x six regions; upper-tail p.
# `grp` is the grouping to test in the six-region test (two levels for a lead; all
# haplotype groups of a gene for the gene-level test).
perm_tests <- function(site, region, focal, grp, SP) {
  idx <- match(site, SP$sites); stopifnot(!anyNA(idx))
  dstat <- function(r) mean(r[focal] == "Desert") - mean(r[!focal] == "Desert")
  obs_d <- dstat(as.character(region))
  obs_c <- chisq_stat(grp, factor(region, levels = REGION_ORDER))
  pd <- numeric(N_PERM); pc <- numeric(N_PERM)
  for (b in seq_len(N_PERM)) {
    lab <- SP$labels[idx, b]
    pd[b] <- dstat(lab)
    pc[b] <- chisq_stat(grp, factor(lab, levels = REGION_ORDER))
  }
  up <- (1 + sum(pd >= obs_d - 1e-12)) / (N_PERM + 1)
  lo <- (1 + sum(pd <= obs_d + 1e-12)) / (N_PERM + 1)
  list(desert_diff = obs_d, p_desert = min(1, 2 * min(up, lo)),
       p_regions = (1 + sum(pc >= obs_c - 1e-12)) / (N_PERM + 1))
}
desert_or <- function(focal, region) {           # accession-level odds ratio (effect size only)
  unname(fisher.test(table(factor(focal, levels = c(TRUE, FALSE)),
                           factor(region == "Desert", levels = c(TRUE, FALSE))))$estimate)
}
region_counts <- function(region) {
  as.integer(table(factor(region, levels = REGION_ORDER)))
}

# ---- Table writer ------------------------------------------------------------------------------
# write_tsv prints some rounded doubles with binary noise (0.6749 -> 0.67490000000000006);
# as.character() gives the shortest form. Used for every results/tables/ output.
write_tab <- function(df, path) {
  df %>% mutate(across(where(is.double), ~ ifelse(is.na(.x), NA_character_, as.character(.x)))) %>%
    write_tsv(path, na = "NA")
}
