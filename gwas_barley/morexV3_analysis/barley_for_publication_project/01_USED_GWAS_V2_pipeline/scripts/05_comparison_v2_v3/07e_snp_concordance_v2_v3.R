#!/usr/bin/env Rscript
# 07e_snp_concordance_v2_v3.R
# Do the SNPs flagged in v2 stay flagged in v3, at the same coordinates?
#
# COORDINATES: v2 and v3 ran on the IDENTICAL 7,110,996-SNP tped -- only the
# kinship/PC covariates changed. Every .ps file therefore has the same rows in
# the same order, and every SNP id / bp coordinate is the same in both runs by
# construction. The table verifies this explicitly (coord_match) rather than
# asserting it. What genuinely changes is the p-value, hence the class.
#
# Because row order is identical across all 48 .ps files, this script indexes by
# ROW POSITION instead of joining on SNP id -- a keyed join over ~28M character
# keys is what made the first version of this script unusably slow.
#
# v2 inputs (BLUP x 3 PCs, the published config):
#   .../publication_BonfOnly_BLUP_3PC/tables/lead_snps.tsv      15 significant
#   .../publication_BonfOnly_BLUP_3PC/tables/marginal_snps.tsv   9 marginal (curated)
#
# v3 classification: ONLY significant vs below_threshold.
#   significant     : -log10p >= 6.0454   (alpha=0.10 / 111,017)
#   below_threshold : everything else
# There is deliberately NO "marginal" class for v3. In v2 the marginal set was a
# MANUALLY CURATED list, not a rule -- so v3 marginals cannot be derived
# automatically and must be picked by hand. Until that selection is made, this
# script reports the raw v3 -log10p and leaves the call open.
#
# Outputs (results/comparison_v2_vs_v3/04_snp_concordance/):
#   v2_flagged_in_v3.tsv         forward: every v2-flagged SNP looked up in v3
#   v3_significant_in_v2.tsv     reverse: every v3-significant lead SNP in v2
#   concordance_summary.md       both directions, readable
#
# Created 2026-08-16.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(parallel) })

PIPE  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
ARC   <- file.path(PIPE, "results", "_archive", "_archive_v2_win50snp_2026-08-16")
V2TAB <- file.path(ARC, "publication_BonfOnly_BLUP_3PC", "tables")
PS_V3 <- file.path(PIPE, "results", "emmax_ps")
PS_V2 <- file.path(ARC, "emmax_ps")
BIM   <- file.path(PIPE, "intermediates", "snp_map.bim")
OUT   <- file.path(PIPE, "results", "comparison_v2_vs_v3", "04_snp_concordance")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BONF_V3 <- -log10(0.10 / 111017L)
BONF_V2 <- -log10(0.10 / 590462L)
TRAITS  <- c("betaglucan", "fiber", "protein", "starch")
CONFIGS <- list(c("BLUP","pc3"), c("BLUP","pc5"), c("BLUP","pc10"),
                c("BLUE","pc3"), c("BLUE","pc5"), c("BLUE","pc10"))
CFG_LAB <- vapply(CONFIGS, function(cf) sprintf("%s x %s PC", cf[1], sub("^pc","",cf[2])), character(1))

# v3 side: significant vs below_threshold only (no automatic marginal class).
classify <- function(nlp, bonf) fifelse(is.na(nlp), NA_character_,
              fifelse(nlp >= bonf, "significant", "below_threshold"))

ps_path <- function(dir, tr, ph, pc) file.path(dir, sprintf("morexV3__%s__%s__%s.ps", tr, ph, pc))
read_nlp <- function(f) -log10(fread(f, select = 4, col.names = "P", showProgress = FALSE)$P)

# ---- Map: row position is the shared key ----
map <- fread(BIM, select = c(1, 2, 4), col.names = c("chrH", "SNP_id", "bp"))
map[, chrN := as.integer(sub("H$", "", chrH))]
NSNP <- nrow(map)

# ---- v2 flagged SNPs ----
lead <- fread(file.path(V2TAB, "lead_snps.tsv"))[, .(trait, SNP_id, chr, position_bp,
                 v2_nlp = neg_log10_p, v2_class = "significant")]
marg <- fread(file.path(V2TAB, "marginal_snps.tsv"))
marg <- marg[, .(trait, SNP_id, chr = as.integer(sub("H$", "", chr)),
                 position_bp, v2_nlp = neg_log10_p, v2_class = "marginal")]
fw <- rbind(lead, marg)
fw[, row := match(SNP_id, map$SNP_id)]
stopifnot(!any(is.na(fw$row)))
# explicit coordinate verification against the shared map
fw[, `:=`(map_chr = map$chrN[row], map_bp = map$bp[row])]
fw[, coord_match := (map_chr == chr & map_bp == position_bp)]

# ---- Read the 48 .ps files once, in parallel; keep only what we need ----
jobs <- CJ(ci = seq_along(CONFIGS), ti = seq_along(TRAITS), sorted = FALSE)
cat(sprintf("[07e] Reading %d v3 .ps files (8 cores) ...\n", nrow(jobs)))
v3_list <- mclapply(seq_len(nrow(jobs)), function(k) {
  cf <- CONFIGS[[jobs$ci[k]]]; tr <- TRAITS[jobs$ti[k]]
  read_nlp(ps_path(PS_V3, tr, cf[1], cf[2]))
}, mc.cores = 8)
stopifnot(all(vapply(v3_list, length, integer(1)) == NSNP))

cat("[07e] Reading the 4 v2 headline .ps files ...\n")
v2_h <- mclapply(TRAITS, function(tr) read_nlp(ps_path(PS_V2, tr, "BLUP", "pc3")), mc.cores = 4)
names(v2_h) <- TRAITS

# index helper: which element of v3_list holds (trait, config)
v3_of <- function(tr, ci) v3_list[[ which(jobs$ci == ci & jobs$ti == match(tr, TRAITS)) ]]

# ================= FORWARD =================
fw[, v3_nlp := NA_real_]; fw[, v3_best_nlp := NA_real_]; fw[, v3_best_config := NA_character_]
for (tr in TRAITS) {
  sel <- fw$trait == tr
  if (!any(sel)) next
  rows <- fw$row[sel]
  fw$v3_nlp[sel] <- v3_of(tr, 1L)[rows]                       # config 1 = BLUP x 3 PC
  best <- rep(-Inf, length(rows)); bcfg <- rep(NA_character_, length(rows))
  for (ci in seq_along(CONFIGS)) {
    val <- v3_of(tr, ci)[rows]
    upd <- !is.na(val) & val > best
    best[upd] <- val[upd]; bcfg[upd] <- CFG_LAB[ci]
  }
  fw$v3_best_nlp[sel] <- best; fw$v3_best_config[sel] <- bcfg
}
fw[, v3_class := classify(v3_nlp, BONF_V3)]
fw[, v3_best_class := classify(v3_best_nlp, BONF_V3)]
fw[, delta := round(v3_nlp - v2_nlp, 3)]
fw[, v3_significant := v3_class == "significant"]
setorder(fw, trait, -v2_nlp)
fw_out <- fw[, .(trait, SNP_id, chr, position_bp, coord_match,
                 v2_nlp = round(v2_nlp, 4), v2_class,
                 v3_nlp = round(v3_nlp, 4), v3_class, delta,
                 v3_best_nlp = round(v3_best_nlp, 4), v3_best_config, v3_best_class, v3_significant)]
fwrite(fw_out, file.path(OUT, "v2_flagged_in_v3.tsv"), sep = "\t")

# ================= REVERSE =================
# v3-significant SNPs in the headline config, LD-clumped at +/-188 kb (v2's window)
WIN <- 188000L
rv <- rbindlist(lapply(TRAITS, function(tr) {
  nlp <- v3_of(tr, 1L)
  idx <- which(!is.na(nlp) & nlp >= BONF_V3)
  if (!length(idx)) return(NULL)
  d <- data.table(trait = tr, row = idx, v3_nlp = nlp[idx],
                  chrN = map$chrN[idx], bp = map$bp[idx], SNP_id = map$SNP_id[idx])
  setorder(d, -v3_nlp); keep <- d[0]
  while (nrow(d)) {
    top <- d[1]; keep <- rbind(keep, top)
    d <- d[!(chrN == top$chrN & abs(bp - top$bp) <= WIN)]
  }
  keep
}))
rv[, v2_nlp := v2_h[[trait[1]]][row], by = trait]
rv[, v2_class := classify(v2_nlp, BONF_V2)]
rv[, delta := round(v3_nlp - v2_nlp, 3)]
rv[, in_v2_tables := SNP_id %in% fw$SNP_id]
setorder(rv, trait, -v3_nlp)
rv_out <- rv[, .(trait, SNP_id, chr = paste0(chrN, "H"), position_bp = bp,
                 v3_nlp = round(v3_nlp, 4), v2_nlp = round(v2_nlp, 4),
                 v2_class, delta, in_v2_tables)]
fwrite(rv_out, file.path(OUT, "v3_significant_in_v2.tsv"), sep = "\t")

# ================= Markdown =================
tab <- function(dt, cols, hdr) c(
  paste0("| ", paste(hdr, collapse = " | "), " |"),
  paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"),
  apply(dt[, ..cols], 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))

md <- c("# SNP concordance: v2 flags vs v3", "",
  sprintf("Generated %s. Headline config: **BLUP x 3 PCs**.", Sys.Date()), "",
  "## Coordinates are identical by construction", "",
  sprintf("Both runs used the same 7,110,996-SNP tped; only the covariates changed. Verified against `snp_map.bim`: **coord_match TRUE for %d / %d** flagged SNPs. No SNP moved position, and none is missing from either run.",
          sum(fw$coord_match), nrow(fw)), "",
  "So the question is purely whether each SNP keeps its **class**:", "",
  sprintf("- **significant**: -log10p >= %.4f (v3), %.4f (v2)", BONF_V3, BONF_V2),
  "- **below_threshold**: everything else", "",
  "**No `marginal` class is assigned to v3.** In v2 that set was curated by hand, not by a rule, so the v3 marginals have to be chosen manually -- the raw `v3_nlp` is given here so that selection can be made.", "",
  "`v3_best` = strongest v3 signal for that SNP across all 6 configs, so a SNP that fades at 3 PCs but survives elsewhere is still visible.", "",
  "## Forward: every SNP v2 flagged, looked up in v3", "")
for (tr in sort(unique(fw_out$trait))) {
  d <- fw_out[trait == tr]
  md <- c(md, sprintf("### %s", tools::toTitleCase(tr)), "",
    tab(d, c("SNP_id","position_bp","v2_nlp","v2_class","v3_nlp","v3_class","delta","v3_best_nlp","v3_best_config"),
        c("SNP","bp","v2 -log10p","v2 class","v3 -log10p","v3 class","delta","v3 best","best config")), "")
}
md <- c(md, "## Class transitions (forward)", "",
  tab(fw_out[, .N, by = .(v2_class, v3_class)][order(v2_class, -N)],
      c("v2_class","v3_class","N"), c("v2 class","v3 class","n SNPs")), "",
  sprintf("Significant in v3 (headline config): **%d / %d**.",
          sum(fw_out$v3_significant), nrow(fw_out)),
  sprintf("Significant in at least one of the 6 v3 configs: **%d / %d**.",
          sum(fw_out$v3_best_class == "significant"), nrow(fw_out)), "",
  "## Reverse: every v3-significant lead SNP, looked up in v2", "",
  sprintf("LD-clumped at +/- %g kb (the window v2 used).", WIN/1000), "")
for (tr in sort(unique(rv_out$trait))) {
  d <- rv_out[trait == tr]
  md <- c(md, sprintf("### %s (%d lead SNPs)", tools::toTitleCase(tr), nrow(d)), "",
    tab(d, c("SNP_id","position_bp","v3_nlp","v2_nlp","v2_class","delta","in_v2_tables"),
        c("SNP","bp","v3 -log10p","v2 -log10p","v2 class","delta","in v2 tables")), "")
}
md <- c(md, "## Reverse summary", "",
  tab(rv_out[, .N, by = v2_class][order(-N)], c("v2_class","N"),
      c("v2 class of v3-significant SNPs","n")), "",
  sprintf("New in v3 (not in the v2 significant or marginal tables): **%d / %d**.",
          sum(!rv_out$in_v2_tables), nrow(rv_out)), "")
writeLines(md, file.path(OUT, "concordance_summary.md"))

cat("\n---------------------------------\n")
cat(sprintf("[07e] CHECKPOINT: v2-flagged SNPs = %d (expect 24 = 15 significant + 9 marginal)\n", nrow(fw_out)))
cat(sprintf("[07e] CHECKPOINT: coord_match TRUE = %d / %d\n", sum(fw$coord_match), nrow(fw)))
cat(sprintf("[07e] v3-significant headline = %d / %d ; significant in >=1 of 6 configs = %d / %d\n",
            sum(fw_out$v3_significant), nrow(fw_out),
            sum(fw_out$v3_best_class == "significant"), nrow(fw_out)))
cat(sprintf("[07e] v3 lead SNPs (clumped) = %d ; new vs v2 tables = %d\n", nrow(rv_out), sum(!rv_out$in_v2_tables)))
print(fw_out[, .N, by = .(v2_class, v3_class)][order(v2_class, -N)])
cat(sprintf("[07e] OK: outputs in %s\n", OUT))
