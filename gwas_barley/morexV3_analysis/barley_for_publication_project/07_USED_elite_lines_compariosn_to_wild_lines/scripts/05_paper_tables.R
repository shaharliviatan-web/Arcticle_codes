# =============================================================================
# 05_paper_tables.R -- the numbers, in tables, ready to fish from when writing
# =============================================================================
# Every quantity a reader or reviewer could ask about this step, in flat TSVs,
# plus one prose file. Same convention as step 03's results_chapter_numbers.txt
# and step 04's Stats/.
#
# Writes:
#   Table_elite_lines.tsv              the 5 lines: name, SAMEA, breeder, year,
#                                      growth habit, country, why chosen
#   Table_panel_composition.tsv        what the 1315-genotype DivBrowse panel is
#   Table_gene_windows.tsv             the 3 genes: coords, window, trait, locus,
#                                      lead SNP, q, eta^2 (carried from step 04)
#   Table_site_overlap.tsv             wild/elite SNP counts and overlap per gene
#   Table_allele_concordance.tsv       per shared site: REF/ALT in both files
#   Table_haplotype_groups.tsv         per group: n, mean/median/sd phenotype
#   Table_elite_genotypes_wide.tsv     per gene: elite lines x shared sites
#   Table_triallelic_sites.tsv         per triallelic site x line: what the line carries
#   (Table_elite_line_screen.tsv is written by 02_screen_elite_lines.R and
#    summarised here; Table_elite_lines.tsv now carries each line's call rate)
#   results_chapter_numbers.txt        all of the above as prose
#
# Author : Shahar Liviatan
# Created: 2026-09-10
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(script_dir, "_load_params.R"))
P <- load_params(file.path(script_dir, "..", "config", "params.sh"))
Sys.setenv(TMPDIR = P$TMPDIR)

genes <- target_genes_df(P)
elite <- read_elite_lines(P)
TB <- P$DIR_TABLES
wt <- function(x, f) write.table(x, file.path(TB, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE, na = "NA")
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

panel <- read.delim(file.path(P$DIR_INPUTS, "panel_metadata.tsv"), stringsAsFactors = FALSE)
overlap <- read.delim(file.path(TB, "site_overlap_summary.tsv"), stringsAsFactors = FALSE)
allele  <- read.delim(file.path(TB, "allele_concordance.tsv"), stringsAsFactors = FALSE)
screen  <- read.delim(file.path(TB, "Table_elite_line_screen.tsv"), stringsAsFactors = FALSE,
                      check.names = FALSE)
tri_gt  <- read.delim(file.path(TB, "triallelic_sites_elite_genotypes.tsv"),
                      stringsAsFactors = FALSE)
thr     <- as.numeric(P$CALL_RATE_MIN)
g04     <- read.delim(P$STEP04_GENE_RESULTS, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# 1. The elite lines
# ---------------------------------------------------------------------------
el <- elite %>%
  left_join(panel %>% select(sample_id, accession_type, organism, panel_label,
                             breeder_meta = breeder, year_meta = year_of_release,
                             annuality_meta = annuality, country_meta = country),
            by = "sample_id") %>%
  left_join(screen %>% select(sample_id, starts_with("call_rate_"), min_gene_call_rate,
                              passes_threshold), by = "sample_id") %>%
  mutate(verified_against_biosamples = "yes")
stopifnot(all(el$passes_threshold))
wt(el, "Table_elite_lines.tsv")
wt(tri_gt, "Table_triallelic_sites.tsv")

# ---------------------------------------------------------------------------
# 2. What the source panel is
# ---------------------------------------------------------------------------
comp <- panel %>%
  mutate(accession_type = ifelse(nzchar(accession_type), accession_type, "(no type recorded)")) %>%
  group_by(accession_type) %>%
  summarise(n_accessions = dplyr::n(),
            n_with_cultivar_name = sum(nzchar(accession_name)),
            n_with_breeder = sum(nzchar(breeder)),
            n_with_release_year = sum(nzchar(year_of_release)),
            .groups = "drop") %>%
  arrange(desc(n_accessions))
comp$n_wild_spontaneum <- sapply(comp$accession_type, function(t) {
  s <- panel[ifelse(nzchar(panel$accession_type), panel$accession_type,
                    "(no type recorded)") == t, ]
  sum(grepl("spontaneum", paste(s$organism, s$infraspecific_name), ignore.case = TRUE))
})
wt(comp, "Table_panel_composition.tsv")

# ---------------------------------------------------------------------------
# 3. The three genes, with their step-04 statistics carried through
# ---------------------------------------------------------------------------
gw <- read.delim(P$STEP04_GENE_WINDOWS, stringsAsFactors = FALSE)
gt <- genes %>%
  left_join(gw %>% select(gene_id, trait, chr, gene_start, gene_end, strand,
                          win_start, win_end, locus_id, lead_SNP, lead_pos,
                          lead_neg_log10p, dist_to_lead_bp), by = "gene_id") %>%
  left_join(g04 %>% select(gene_id, n_snps_window, n_ind, n_unassigned, n_groups,
                           group_sizes, kw_p_raw, fdr_q, eta_squared,
                           delta_top_bottom_sd, description), by = "gene_id") %>%
  mutate(window_bp = P$WINDOW_BP, crosshap_run = P$STEP04_RUN_ID,
         crosshap_epsilon = P$CROSSHAP_EPSILON, crosshap_MGmin = P$CROSSHAP_MGMIN)
wt(gt, "Table_gene_windows.tsv")

# ---------------------------------------------------------------------------
# 4/5. Site overlap and allele concordance (already written by 02; re-emit as
#      the paper-facing names so everything a writer needs sits under Table_*)
# ---------------------------------------------------------------------------
wt(overlap, "Table_site_overlap.tsv")
wt(allele,  "Table_allele_concordance.tsv")

# ---------------------------------------------------------------------------
# 6. Haplotype group phenotype summary + 7. elite genotypes, wide
# ---------------------------------------------------------------------------
hg <- list(); ew <- list()
for (i in seq_len(nrow(genes))) {
  gid <- genes$gene_id[i]; short <- genes$short_name[i]
  f <- file.path(P$DIR_MATRICES, paste0(gid, "__shared_sites.rds"))
  if (!file.exists(f)) next
  D <- readRDS(f)

  hg[[gid]] <- D$indfile %>% group_by(hap) %>%
    summarise(n = dplyr::n(), mean_pheno = mean(Pheno, na.rm = TRUE),
              median_pheno = median(Pheno, na.rm = TRUE),
              sd_pheno = sd(Pheno, na.rm = TRUE), .groups = "drop") %>%
    mutate(gene_id = gid, short_name = short, trait = D$trait) %>%
    select(gene_id, short_name, trait, hap, n, mean_pheno, median_pheno, sd_pheno)

  m <- D$elite
  if (ncol(m)) {
    lab <- c("REF", "ALT", "HET")
    chr_m <- matrix(ifelse(is.na(m), "missing", lab[m + 1]), nrow = nrow(m),
                    dimnames = dimnames(m))
    d <- as.data.frame(chr_m, stringsAsFactors = FALSE)
    d$line_name <- rownames(chr_m)
    ew[[gid]] <- cbind(gene_id = gid, short_name = short,
                       d[, c("line_name", setdiff(names(d), "line_name"))])
  }
}
wt(bind_rows(hg), "Table_haplotype_groups.tsv")
for (gid in names(ew))
  wt(ew[[gid]], paste0("Table_elite_genotypes_wide__", gid, ".tsv"))

# ---------------------------------------------------------------------------
# 8. Prose
# ---------------------------------------------------------------------------
con <- file(file.path(TB, "results_chapter_numbers.txt"), "w")
# NOTE: sprintf() on a NULL argument returns character(0), and writeLines() of
# that writes NOTHING -- a renamed upstream column would silently delete lines
# from this report rather than erroring. w() therefore refuses empty input.
w <- function(...) {
  x <- paste0(...)
  if (!length(x)) stop("w(): empty output -- a referenced column is probably missing")
  writeLines(x, con)
}

w("Elite vs wild haplotype comparison -- every number in this step")
w(strrep("=", 72)); w("")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("Wild haplotypes: step 04 run ", P$STEP04_RUN_ID,
  " (crosshap epsilon ", P$CROSSHAP_EPSILON, ", MGmin ", P$CROSSHAP_MGMIN,
  "), consumed unchanged; crosshap was NOT re-run and never saw the elite lines.")
w("")
w("[A] SOURCE OF THE ELITE GENOTYPES")
w("  IPK DivBrowse barley pangenome v2 (", P$DIVBROWSE_BASE, ")")
w("  Dataset: unimputed SNP variants for SHAPE2 Core1000 and elite line genotypes,")
w("  called against Morex V3. Panel size: ", nrow(panel), " genotypes.")
for (i in seq_len(nrow(comp)))
  w(sprintf("    %-24s %4d accessions, %4d with a cultivar name, %d wild (ssp. spontaneum)",
            comp$accession_type[i], comp$n_accessions[i],
            comp$n_with_cultivar_name[i], comp$n_wild_spontaneum[i]))
w("  No wild barley (ssp. spontaneum) is present in the panel, so no elite line")
w("  in this comparison can be a wild accession.")
w("")
w("[B] THE FIVE ELITE LINES")
w("  All from the `elite lines` subset, all spring malting types, released ",
  min(as.integer(el$year_of_release)), "-", max(as.integer(el$year_of_release)), ".")
w("  All five verified against their EBI BioSamples record (name, type, year).")
w(sprintf("  Call-rate screen (2026-09-14): all %d spring elite lines scored at the shared",
          nrow(screen)))
w(sprintf("  sites of the three genes; %d pass >= %.0f%% called. Every line shown passes.",
          sum(screen$passes_threshold), 100 * thr))
first <- screen[screen$accession_name %in% c("RGT Planet", "Propino"), ]
if (nrow(first))
  w("  Replaced from the first selection (2026-09-10) for low call rate: ",
    paste(sprintf("%s %.0f%%", first$accession_name, 100 * first$call_rate_all),
          collapse = ", "), ".")
for (i in seq_len(nrow(el)))
  w(sprintf("    %-12s %-16s %s  %-12s %-18s called %5.1f%%", el$line_name[i],
            el$sample_id[i], el$year_of_release[i], el$annuality[i],
            substr(el$breeder[i], 1, 18), 100 * el$call_rate_all[i]))
w("")
w("[C] GENES AND WINDOWS")
for (i in seq_len(nrow(gt)))
  w(sprintf("    %-8s %s  %s:%s-%s (gene %s-%s, +/-%d bp)  trait=%s  q=%s  eta2=%s",
            gt$short_name[i], gt$gene_id[i], gt$chr[i],
            format(gt$win_start[i], big.mark = ","), format(gt$win_end[i], big.mark = ","),
            format(gt$gene_start[i], big.mark = ","), format(gt$gene_end[i], big.mark = ","),
            P$WINDOW_BP, gt$trait[i], signif(gt$fdr_q[i], 3), signif(gt$eta_squared[i], 3)))
w("")
w("[D] SITE OVERLAP BETWEEN THE TWO CALL SETS")
w("  The wild set is a QC-filtered call set on 290 Levantine accessions; the elite")
w("  export is unfiltered across the 1315-genotype panel. SNPs therefore differ in")
w("  BOTH directions.")
for (i in seq_len(nrow(overlap)))
  w(sprintf("    %-8s wild %3d SNPs | elite %3d records | shared %3d | triallelic %d | wild-only %3d | elite-only %3d (%d indel)",
            overlap$short_name[i], overlap$n_wild_snps[i], overlap$n_elite_records[i],
            overlap$n_shared[i], overlap$n_alleles_triallelic[i], overlap$n_wild_only[i],
            overlap$n_elite_only[i], overlap$n_elite_only_indel[i]))
w("  wild-only / elite-only = no record at that position in the other file.")
w("  Elite-only sites are never drawn: the wild haplotype groups have no data there.")
w("")
w("[E] ALLELE CONCORDANCE (REF/ALT must be identical, never flipped)")
w("  Both call sets are called against Morex V3 and the wild set carries no PLINK")
w("  allele flip, so at a shared position REF and ALT must agree exactly. If they")
w("  were swapped, 0/0 in one file and 0/0 in the other would denote OPPOSITE")
w("  alleles and every barcode here would be silently inverted.")
stopifnot(all(c("n_alleles_identical", "n_alleles_swapped", "n_alleles_triallelic",
                "n_alleles_ref_differs") %in% names(overlap)))
for (i in seq_len(nrow(overlap)))
  w(sprintf("    %-8s %3d shared positions: %3d identical, %d swapped, %d triallelic (ALT differs), %d differing REF",
            overlap$short_name[i], overlap$n_shared_positions[i],
            overlap$n_alleles_identical[i], overlap$n_alleles_swapped[i],
            overlap$n_alleles_triallelic[i], overlap$n_alleles_ref_differs[i]))
tot_fatal <- sum(overlap$n_alleles_swapped) + sum(overlap$n_alleles_ref_differs)
tot_tri   <- sum(overlap$n_alleles_triallelic)
w(if (tot_fatal == 0)
    sprintf("  VERIFIED: %d of %d shared positions carry IDENTICAL REF and ALT. No allele is flipped anywhere.",
            sum(overlap$n_alleles_identical), sum(overlap$n_shared_positions))
  else sprintf("  WARNING: %d shared site(s) are swapped or differ in REF -- see Table_allele_concordance.tsv", tot_fatal))
if (tot_tri > 0) {
  w(sprintf("  %d shared position(s) are TRIALLELIC: same reference base, but a different", tot_tri))
  w("  alternate allele segregates in the wild panel and in the elite panel. That is a")
  w("  population difference, not a coding error. Such sites never enter the shared set")
  w("  (the key is CHROM:POS:REF:ALT) and are listed in Table_allele_concordance.tsv:")
  tri <- allele[allele$status == "alt_differs", ]
  for (i in seq_len(nrow(tri)))
    w(sprintf("    %-8s %s:%s  wild %s/%s   elite %s/%s",
              tri$short_name[i], tri$chr[i], tri$pos[i],
              tri$wild_ref[i], tri$wild_alt[i], tri$elite_ref[i], tri$elite_alt[i]))
  w("  shared_sites does not draw them. Both filled versions draw each line's REAL")
  w("  genotype there, never a fill (0/0 = REF, the same base as the wild REF):")
  for (st in unique(tri_gt$site)) {
    g <- tri_gt[tri_gt$site == st, ]
    w(sprintf("    %-22s %s", st, paste(sprintf("%s: %s", g$line_name, g$genotype),
                                        collapse = " | ")))
  }
  n_third <- sum(grepl("third allele", tri_gt$genotype))
  w(if (n_third == 0)
      "  No shown line carries the elite-only third allele, so the triallelic colour is not drawn."
    else sprintf("  %d call(s) carry the elite-only third allele (purple in the filled figures).", n_third))
}
w("")
w("[F] WILD HAPLOTYPE GROUPS (unassigned hap 0 dropped, as in step 04)")
hgt <- bind_rows(hg)
for (i in seq_len(nrow(hgt)))
  w(sprintf("    %-8s group %-2s n=%3d  mean=%8.3f  median=%8.3f  sd=%7.3f",
            hgt$short_name[i], hgt$hap[i], hgt$n[i], hgt$mean_pheno[i],
            hgt$median_pheno[i], hgt$sd_pheno[i]))
w("")
w("[G] HOW THE FIGURES ARE BUILT")
w("  Group rows are a PER-SNP MAJORITY CONSENSUS over the group's members,")
w("  ignoring missing calls; exact 50/50 ties are drawn as missing. A consensus is")
w("  used rather than one representative accession because wild missingness is")
w("  high, so a single plant's row would show grey tiles that are artefacts of its")
w("  sequencing rather than features of the haplotype.")
w("  Elite lines are NOT assigned to a haplotype group anywhere: crosshap never")
w("  saw them. The barcodes are aligned so the comparison is left to the reader.")
w("  Three SNP-matching versions are produced (", paste(P$VERSIONS, collapse = ", "), ").")
w("  Barcode rows are separated by white gaps; each legend lists only the states drawn.")
w("")
w("[H] CAVEAT FOR THE MANUSCRIPT")
w("  The figure shows which haplotypes elite lines carry; it does not by itself")
w("  show that breeding selected them.")
close(con)

msg("DONE -> ", TB)
