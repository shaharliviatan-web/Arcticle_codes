#!/usr/bin/env Rscript
# ============================================================================
# 07_fiber_starch_tradeoff_direction / 01_tradeoff_direction.R
#
# Question tested (per shared 7H candidate gene):
#   Q1. Are the crosshap haplotype groups IDENTICAL in the fiber run and the
#       starch run? (crosshap clusters on genotype only, so they must be -- this
#       verifies it empirically, individual by individual.)
#   Q2. Is the highest-mean-FIBER haplotype group the lowest-mean-STARCH group?
#       i.e. does the SAME haplotype move fiber up and starch down (tradeoff
#       direction)?
#
# Scope: the 3 genes that passed the FDR candidate-gene filter for BOTH fiber
# and starch at the 7H ~573.6 Mb locus (step 08 annotation master):
#   HORVU.MOREX.r3.7HG0729020, ...7HG0729090, ...7HG0729100
#
# Inputs (read-only):
#   ../04_runs/candidate_genes_1000bp_V1/Cache/{fiber,starch}/<gene>/MGmin_2/HapObject.rds
#       crosshap objects; $HapObject[[<Haplotypes_MGmin#_E#>]]$Indfile has, per
#       accession: Ind, hap (0 = unassigned; A/B/.. = haplotype group), Pheno
#       (the BLUP of THAT run's trait).
#   /mnt/data/shahar/gwas_barley/data/inputs/{fiber,starch}_corrected_V3.pheno
#       BLUP phenotypes, 290 accessions (overall correlation anchor only).
#
# Per-gene haplotype partition used = the per-gene BEST (most significant)
# crosshap setting from ../04_runs/.../Stats/per_test_stats.csv
# (is_best_within_gene == TRUE). For all 3 genes fiber and starch agree on the
# same best setting (they must -- the partition is shared), so one setting/gene.
#
# Outputs (results/tables/):
#   shared_genes.tsv               the 3 genes + lead SNPs + class per trait
#   haplotype_group_means.tsv      per gene x hap group: n, mean/sd/median fiber & starch
#   tradeoff_summary.tsv           per gene: grouping-identity, KW p, corr, direction call
#   overall_phenotypic_correlation.txt
#
# Conventions: TMPDIR pinned; absolute paths; ASCII only; no /tmp or $HOME writes.
# ============================================================================

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
options(stringsAsFactors = FALSE)

step_dir  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/07_fiber_starch_tradeoff_direction"
cache_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/04_runs/candidate_genes_1000bp_V1/Cache"
out_dir   <- file.path(step_dir, "results", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- shared genes, best per-gene crosshap setting, annotation + GWAS context ----
genes <- data.frame(
  gene_id      = c("HORVU.MOREX.r3.7HG0729020",
                   "HORVU.MOREX.r3.7HG0729090",
                   "HORVU.MOREX.r3.7HG0729100"),
  best_key     = c("Haplotypes_MGmin2_E0.8",
                   "Haplotypes_MGmin2_E0.85",
                   "Haplotypes_MGmin2_E0.85"),
  MGmin        = c(2L, 2L, 2L),
  epsilon      = c(0.8, 0.85, 0.85),
  annotation   = c("Probable high-affinity nitrate transporter 2.4",
                   "Pentatricopeptide repeat-containing protein At3g61360",
                   "NifU-like protein 3, chloroplastic"),
  evidence_tier= c("HIGH", "MEDIUM", "HIGH"),
  fiber_leadSNP  = "7H:573606306",
  fiber_class    = "significant",
  starch_leadSNP = "7H:573606460",
  starch_class   = "marginal"
)

read_ind <- function(trait, gene_id, key) {
  rds <- file.path(cache_dir, trait, gene_id, "MGmin_2", "HapObject.rds")
  x   <- readRDS(rds)
  ind <- x$HapObject[[key]]$Indfile
  ind[, c("Ind", "hap", "Pheno")]
}

group_means <- list()
summ        <- list()

for (i in seq_len(nrow(genes))) {
  g   <- genes$gene_id[i]
  key <- genes$best_key[i]

  fi <- read_ind("fiber",  g, key)
  si <- read_ind("starch", g, key)
  m  <- merge(fi, si, by = "Ind", suffixes = c("_fib", "_sta"))

  # ---- Q1: identical grouping? ----
  identical_grouping <- all(m$hap_fib == m$hap_sta)
  n_mismatch         <- sum(m$hap_fib != m$hap_sta)

  # assigned individuals only (drop hap == "0" = unassigned)
  m2 <- m[m$hap_fib != "0", ]
  m2$hap <- m2$hap_fib

  ag <- do.call(rbind, lapply(split(m2, m2$hap), function(d) {
    data.frame(
      gene_id      = g,
      annotation   = genes$annotation[i],
      setting      = sprintf("MGmin%d_eps%.2f", genes$MGmin[i], genes$epsilon[i]),
      hap          = d$hap[1],
      n            = nrow(d),
      mean_fiber   = mean(d$Pheno_fib),
      sd_fiber     = sd(d$Pheno_fib),
      median_fiber = median(d$Pheno_fib),
      mean_starch  = mean(d$Pheno_sta),
      sd_starch    = sd(d$Pheno_sta),
      median_starch= median(d$Pheno_sta)
    )
  }))
  ag <- ag[order(-ag$mean_fiber), ]
  group_means[[g]] <- ag

  # ---- Q2: direction ----
  pear <- cor(ag$mean_fiber, ag$mean_starch)
  spear<- cor(ag$mean_fiber, ag$mean_starch, method = "spearman")
  top_fiber_hap   <- ag$hap[1]
  top_fiber_starch_rank <- rank(ag$mean_starch)[1]   # 1 = lowest starch
  n_groups        <- nrow(ag)
  direction_inverse <- (top_fiber_starch_rank == 1)

  # Kruskal-Wallis on the shared partition, each trait (matches crosshap test)
  kw_fib <- kruskal.test(Pheno_fib ~ hap, data = m2)$p.value
  kw_sta <- kruskal.test(Pheno_sta ~ hap, data = m2)$p.value

  summ[[g]] <- data.frame(
    gene_id            = g,
    annotation         = genes$annotation[i],
    setting            = sprintf("MGmin%d_eps%.2f", genes$MGmin[i], genes$epsilon[i]),
    n_individuals      = nrow(m),
    n_assigned         = nrow(m2),
    n_groups           = n_groups,
    grouping_identical = identical_grouping,
    n_mismatch         = n_mismatch,
    kw_p_fiber         = kw_fib,
    kw_p_starch        = kw_sta,
    top_fiber_hap      = top_fiber_hap,
    top_fiber_starch_rank = top_fiber_starch_rank,
    direction_inverse  = direction_inverse,
    pearson_groupmeans = pear,
    spearman_groupmeans= spear
  )
}

gm <- do.call(rbind, group_means)
sm <- do.call(rbind, summ)

# round numeric cols for readable tables
num3 <- function(x) round(x, 3)
gm[, c("mean_fiber","sd_fiber","median_fiber","mean_starch","sd_starch","median_starch")] <-
  lapply(gm[, c("mean_fiber","sd_fiber","median_fiber","mean_starch","sd_starch","median_starch")], num3)
sm$pearson_groupmeans  <- num3(sm$pearson_groupmeans)
sm$spearman_groupmeans <- num3(sm$spearman_groupmeans)
sm$kw_p_fiber  <- signif(sm$kw_p_fiber, 3)
sm$kw_p_starch <- signif(sm$kw_p_starch, 3)

# ---- overall phenotypic correlation anchor ----
fp <- read.table("/mnt/data/shahar/gwas_barley/data/inputs/fiber_corrected_V3.pheno")
sp <- read.table("/mnt/data/shahar/gwas_barley/data/inputs/starch_corrected_V3.pheno")
ph <- merge(fp[, c(1, 3)], sp[, c(1, 3)], by = "V1")
colnames(ph) <- c("id", "fiber", "starch")
ov_p <- cor(ph$fiber, ph$starch)
ov_s <- cor(ph$fiber, ph$starch, method = "spearman")

# ---- write outputs ----
write.table(genes[, c("gene_id","annotation","evidence_tier",
                      "fiber_leadSNP","fiber_class","starch_leadSNP","starch_class",
                      "MGmin","epsilon")],
            file.path(out_dir, "shared_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(gm, file.path(out_dir, "haplotype_group_means.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sm, file.path(out_dir, "tradeoff_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "Overall phenotypic correlation, fiber vs starch BLUP (290 accessions)",
  sprintf("n            = %d", nrow(ph)),
  sprintf("Pearson r    = %.3f", ov_p),
  sprintf("Spearman rho = %.3f", ov_s)
), file.path(out_dir, "overall_phenotypic_correlation.txt"))

cat("\n==== shared_genes.tsv ====\n");            print(genes, row.names = FALSE)
cat("\n==== haplotype_group_means.tsv ====\n");   print(gm, row.names = FALSE)
cat("\n==== tradeoff_summary.tsv ====\n");        print(sm, row.names = FALSE)
cat(sprintf("\n==== overall phenotypic corr: Pearson %.3f / Spearman %.3f (n=%d) ====\n",
            ov_p, ov_s, nrow(ph)))
cat("\nDone. Tables in", out_dir, "\n")
