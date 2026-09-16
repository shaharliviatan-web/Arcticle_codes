#!/usr/bin/env Rscript
# 07d_hit_table_readable.R
# Turn 03_manhattan/v2_vs_v3/hit_counts_compare.tsv into human-readable tables.
#
# Writes into results/comparison_v2_vs_v3/03_manhattan/:
#   hit_counts_readable.md   markdown, one section per trait, configs as rows
#   hit_counts_readable.txt  same content as fixed-width plain text
#
# Column meanings (the three scorings matter and are easy to conflate):
#   v3 top / v3 hits  : v3 run, scored at the v3 threshold 6.0454
#   v2 top / v2 hits  : v2 run, scored at the v2 threshold 6.7712  (as v2 was published)
#   v2 hits @6.045    : v2 run RESCORED at the v3 threshold
#                       -> compare THIS against "v3 hits" to isolate the effect of the
#                          new correction from the effect of the easier threshold.
#
# Created 2026-08-16.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))

PIPE <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
DIR  <- file.path(PIPE, "results", "comparison_v2_vs_v3", "03_manhattan")
h    <- fread(file.path(DIR, "v2_vs_v3", "hit_counts_compare.tsv"))

BONF_V3 <- 6.0454; BONF_V2 <- 6.7712
HEADLINE <- quote(pheno_type == "BLUP" & n_PCs == 3L)

h[, config := sprintf("%s x %d PC", pheno_type, n_PCs)]
h[eval(HEADLINE), config := paste0(config, " (headline)")]
# display order: BLUP 3,5,10 then BLUE 3,5,10
h[, ord := match(pheno_type, c("BLUP", "BLUE")) * 100 + n_PCs]
setorder(h, trait, ord)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

# ---------------- Markdown ----------------
md <- c(
  "# v3 vs v2 — hits and top signal, all traits x all configs",
  "",
  sprintf("Generated %s from `v2_vs_v3/hit_counts_compare.tsv`.", format(Sys.Date())),
  "",
  "**How to read this.** Each version is first scored at its own Bonferroni line",
  sprintf("(v3 = %.4f on 111,017 pruned SNPs; v2 = %.4f on 590,462). The last column", BONF_V3, BONF_V2),
  "rescores the **v2** run at the **v3** threshold — compare it against `v3 hits`",
  "to see what the new correction did, with the threshold change held constant.",
  "")

for (tr in sort(unique(h$trait))) {
  d <- h[trait == tr]
  md <- c(md,
    sprintf("## %s", tools::toTitleCase(tr)), "",
    "| config | v3 top -log10p | v3 hits | v2 top -log10p | v2 hits | v2 hits @6.045 |",
    "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(d))) {
    md <- c(md, sprintf("| %s | %s | %d | %s | %d | %d |",
                        d$config[i], fmt(d$v3_top[i]), d$v3_hits_at_v3thr[i],
                        fmt(d$v2_top[i]), d$v2_hits_at_v2thr[i], d$v2_hits_at_v3thr[i]))
  }
  best <- d[which.max(v3_top)]
  md <- c(md, "",
    sprintf("Strongest v3 signal for %s: **%s** at -log10p = **%s** (%d hits).",
            tr, sub(" \\(headline\\)", "", best$config), fmt(best$v3_top), best$v3_hits_at_v3thr),
    "")
}

tot <- h[, .(v3 = sum(v3_hits_at_v3thr), v2own = sum(v2_hits_at_v2thr),
             v2at3 = sum(v2_hits_at_v3thr))]
md <- c(md,
  "## Totals across all 24 cells", "",
  "| scoring | hits |", "|---|---|",
  sprintf("| v3 at its own threshold (%.4f) | %d |", BONF_V3, tot$v3),
  sprintf("| v2 at its own threshold (%.4f) | %d |", BONF_V2, tot$v2own),
  sprintf("| v2 rescored at the v3 threshold (%.4f) | %d |", BONF_V3, tot$v2at3),
  "",
  sprintf(paste("Note: v2 rescored at the v3 line yields **%d** hits versus v3's **%d**.",
                "The headline jump from %d to %d is therefore driven by the easier",
                "threshold, not by the new correction. See the README for why the",
                "correction is still the better-calibrated model (lambda)."),
          tot$v2at3, tot$v3, tot$v2own, tot$v3),
  "")
writeLines(md, file.path(DIR, "hit_counts_readable.md"))

# ---------------- Plain text ----------------
txt <- c("v3 vs v2 - hits and top signal, all traits x all configs",
         strrep("=", 78), "")
for (tr in sort(unique(h$trait))) {
  d <- h[trait == tr]
  txt <- c(txt, toupper(tr), strrep("-", 78),
    sprintf("%-24s %14s %8s %14s %8s %15s",
            "config", "v3 top", "v3 hits", "v2 top", "v2 hits", "v2 hits@6.045"))
  for (i in seq_len(nrow(d)))
    txt <- c(txt, sprintf("%-24s %14s %8d %14s %8d %15d",
                          d$config[i], fmt(d$v3_top[i]), d$v3_hits_at_v3thr[i],
                          fmt(d$v2_top[i]), d$v2_hits_at_v2thr[i], d$v2_hits_at_v3thr[i]))
  txt <- c(txt, "")
}
txt <- c(txt, strrep("=", 78),
  sprintf("TOTALS   v3@%.4f = %d    v2@%.4f = %d    v2@%.4f = %d",
          BONF_V3, tot$v3, BONF_V2, tot$v2own, BONF_V3, tot$v2at3), "")
writeLines(txt, file.path(DIR, "hit_counts_readable.txt"))

cat(sprintf("[07d] OK: hit_counts_readable.{md,txt} in %s\n", DIR))
cat(sprintf("[07d] totals: v3=%d  v2@own=%d  v2@v3thr=%d\n", tot$v3, tot$v2own, tot$v2at3))
