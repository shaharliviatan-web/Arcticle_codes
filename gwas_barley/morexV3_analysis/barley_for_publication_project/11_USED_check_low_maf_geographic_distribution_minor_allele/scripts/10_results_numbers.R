# 10_results_numbers.R — every number of this step, as prose, for the Discussion, the M&M
# and the Online Resource captions. Written from the tables of scripts 02–07 only, so it can
# be regenerated after any re-run. Nothing is computed here that a table does not hold.
#
# Output: results/results_numbers.txt

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))

rd <- function(f) read_tsv(file.path(TAB, f), show_col_types = FALSE)
t02 <- rd("T02_lead_enrichment_tests.tsv"); t04 <- rd("T04_haplotype_group_enrichment_tests.tsv")
t06 <- rd("T06_top_region_vs_background.tsv"); t07 <- rd("T07_background_top_region_by_MAF.tsv")
t08 <- rd("T08_lead_within_site.tsv"); t09 <- rd("T09_lead_carrier_sharing.tsv")
t10 <- rd("T10_accession_minor_allele_load.tsv"); t11 <- rd("T11_kinship_by_region.tsv")
g   <- read_tsv(file.path(INTER, "lead_genotype_long.tsv"), show_col_types = FALSE)

out <- character(); say <- function(...) out <<- c(out, paste0(...))
pct <- function(x) sprintf("%.0f%%", 100 * x)
say("11_USED_check_low_maf_geographic_distribution_minor_allele - results numbers")
say("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by scripts/10_results_numbers.R from results/tables/T01-T11.")
say("Site-permutation tests: ", N_PERM, " permutations of the region labels among the 29 sites, seed ", SEED, ".")
say("")

# ---- 1. leads: Desert and six-region tests ---------------------------------------------------
say("== 1. Lead SNPs (36): origin of minor- vs major-allele carriers (T01, T02) ==")
for (t in TRAIT_ORDER) {
  x <- t02 %>% filter(trait == t)
  say(sprintf("%s: %d loci. Desert test q <= 0.05 at %d (enriched %d, depleted %d); six-region test q <= 0.05 at %d. Top region of the minor-allele carriers = Desert at %d.",
              TRAIT_LAB[t], nrow(x), sum(x$q_desert <= 0.05), sum(x$q_desert <= 0.05 & x$desert_direction == "enriched"),
              sum(x$q_desert <= 0.05 & x$desert_direction == "depleted"), sum(x$q_regions <= 0.05),
              sum(apply(x[, c("minor_North", "minor_Coast", "minor_Desert", "minor_HZ1", "minor_HZ2", "minor_HZ3")], 1, which.max) == 3)))
}
inst <- g %>% filter(allele == "minor") %>% group_by(trait) %>%
  summarise(n = n(), d = mean(region == "Desert"), dh = mean(region %in% c("Desert", "HZ3 (Coast-Desert)")), .groups = "drop")
acc <- g %>% distinct(IID, region)
pd <- mean(acc$region == "Desert"); pdh <- mean(acc$region %in% c("Desert", "HZ3 (Coast-Desert)"))
for (i in seq_len(nrow(inst)))
  say(sprintf("%s: %d minor-allele instances; %s from the Desert, %s from the Desert + HZ3 (panel: %s and %s).",
              TRAIT_LAB[inst$trait[i]], inst$n[i], pct(inst$d[i]), pct(inst$dh[i]), pct(pd), pct(pdh)))
fb <- t02 %>% filter(trait %in% c("fiber", "betaglucan"), desert_pct_minor > desert_pct_major)
say(sprintf("Fiber + beta-glucan leads whose minor-allele carriers are more often from the Desert than the major-allele carriers (any q): %d of 29; of these, desert carriers from >= 2 Desert sites at %d.",
            nrow(fb), sum(fb$n_desert_sites_minor >= 2)))
say("Desert test q <= 0.05, enriched: ", paste(t02$locus_id[t02$q_desert <= 0.05 & t02$desert_direction == "enriched"], collapse = ", "))
say("Desert test q <= 0.05, depleted: ", paste(t02$locus_id[t02$q_desert <= 0.05 & t02$desert_direction == "depleted"], collapse = ", "))
for (s in c("fiber_L17", "starch_L05")) {
  x <- t02 %>% filter(locus_id == s)
  say(sprintf("7H %s (%s): %d carriers from %d sites in %d regions (N/C/D/HZ1/HZ2/HZ3 = %d/%d/%d/%d/%d/%d); Desert %.1f%% of minor vs %.1f%% of major carriers, q Desert = %s, q regions = %s.",
              s, x$lead_SNP, x$n_minor, x$n_sites_minor, x$n_regions_minor, x$minor_North, x$minor_Coast, x$minor_Desert,
              x$minor_HZ1, x$minor_HZ2, x$minor_HZ3, x$desert_pct_minor, x$desert_pct_major, x$q_desert, x$q_regions))
}
say("")

# ---- 2. matched background ------------------------------------------------------------------
say("== 2. Matched genome-wide background (T05-T07) ==")
for (grp in c("fiber", "betaglucan", "fiber+betaglucan", "starch", "all_but_protein")) {
  x <- t06 %>% filter(leads == grp)
  d <- x %>% filter(top_region == "Desert"); n <- x %>% filter(top_region == "North")
  say(sprintf("%s (%d leads): top region Desert %d vs %.2f expected (p = %s); North %d vs %.2f expected (p fewer = %s).",
              grp, d$n_leads, d$observed, d$expected, d$p_more, n$observed, n$expected, n$p_fewer))
}
r <- t07 %>% filter(maf_class == "MAF < 0.10")
say("Background SNPs with MAF < 0.10, top region: ", paste0(REGION_SHORT[match(r$top_region, REGION_ORDER)], " ", r$pct, "%", collapse = ", "), ".")
say("")

# ---- 3. within site ------------------------------------------------------------------------------
say("== 3. Carriers vs their own site-mates (T08) ==")
for (t in TRAIT_ORDER) {
  x <- t08 %>% filter(trait == t)
  say(sprintf("%s: median share of the raw carrier difference predicted by the carriers' sites %.0f%%; median carriers beyond their site-mates %.0f%%; median raw difference kept without the 2 most extreme accessions %.0f%%.",
              TRAIT_LAB[t], 100 * median(x$frac_expected_site), median(x$pct_carriers_beyond_sitemates), median(x$pct_diff_kept_without_top2)))
}
few <- t08 %>% filter(n_minor_with_sitemates / n_minor < 0.6)
say(sprintf("Loci where < 60%% of carriers have a non-carrier site-mate (%d): %s.", nrow(few),
            paste0(few$locus_id, ifelse(!is.na(few$gene), paste0(" [", few$gene, "]"), ""), " ", few$n_minor_with_sitemates, "/", few$n_minor, collapse = "; ")))
say("Raw carrier difference opposite in sign to the GWAS beta: ", paste(t08$locus_id[!t08$raw_sign_agrees_with_beta], collapse = ", "), ".")
say("")

# ---- 4. carrier sharing -----------------------------------------------------------------------
say("== 4. Shared carriers across loci (T09, T10) ==")
x <- t09 %>% filter(!same_chr) %>% group_by(trait) %>% summarise(n = n(), s = sum(p_hyper < 0.001), .groups = "drop")
say("Pairs of loci on different chromosomes sharing more carriers than chance (hypergeometric p < 0.001): ",
    paste0(x$trait, " ", x$s, "/", x$n, collapse = ", "), ".")
top <- t10 %>% slice_max(minor_fiber, n = 3, with_ties = FALSE)
for (i in seq_len(nrow(top)))
  say(sprintf("%s (%s, %s): minor allele at %d/18 fiber, %d/11 beta-glucan, %d/5 starch loci; fiber rank %d (high), beta-glucan rank %d, starch rank %d (low); genome-wide rare-allele load rank %d of 290.",
              top$IID[i], top$location[i], top$region[i], top$minor_fiber[i], top$minor_betaglucan[i], top$minor_starch[i],
              top$fiber_rank_high[i], top$betaglucan_rank_high[i], top$starch_rank_low[i], top$gw_rare_rate_rank[i]))
say(sprintf("Accessions carrying fiber minor alleles at >= 3 loci: %d of 50 Desert, %d of 240 others.",
            sum(t10$region == "Desert" & t10$minor_fiber >= 3), sum(t10$region != "Desert" & t10$minor_fiber >= 3)))
say("")

# ---- 5. haplotype groups ------------------------------------------------------------------------
say("== 5. Haplotype groups of the four presented genes (T03, T04) ==")
for (gn in unique(t04$gene)) {
  x <- t04 %>% filter(gene == gn)
  say(sprintf("%s (%s): %d groups over %d assigned accessions; groups x regions q = %s.", gn, x$trait[1], nrow(x),
              x$n_assigned_gene[1], x$q_regions_gene[1]))
  for (i in seq_len(nrow(x)))
    say(sprintf("   %s n = %d, mean BLUP %s: N/C/D/HZ1/HZ2/HZ3 = %d/%d/%d/%d/%d/%d, %d sites; Desert %.1f%% vs %.1f%% of the other groups (%s, q = %s).",
                x$group[i], x$n[i],
                paste(na.omit(c(if (!is.na(x$mean_BLUP_fiber[i])) paste0("fiber ", x$mean_BLUP_fiber[i]),
                                if (!is.na(x$mean_BLUP_starch[i])) paste0("starch ", x$mean_BLUP_starch[i]))), collapse = ", "),
                x$North[i], x$Coast[i], x$Desert[i], x$HZ1[i], x$HZ2[i], x$HZ3[i], x$n_sites[i],
                x$desert_pct_group[i], x$desert_pct_other_groups[i], x$desert_direction[i], x$q_desert[i]))
}
say("")

# ---- 6. kinship --------------------------------------------------------------------------------
say("== 6. Relatedness within and between regions, aIBS (T11) ==")
w <- t11 %>% filter(pair_type == "different sites, same region")
s <- t11 %>% filter(pair_type == "same site")
say("Mean kinship between accessions of different sites of the same region: ",
    paste0(REGION_SHORT[match(w$region_1, REGION_ORDER)], " ", w$mean_kinship, collapse = ", "), ".")
say("Mean kinship within sites: ", paste0(REGION_SHORT[match(s$region_1, REGION_ORDER)], " ", s$mean_kinship, collapse = ", "), ".")
say("All pairs: ", t11$mean_kinship[t11$pair_type == "all pairs"], ".")

writeLines(out, file.path(HERE, "results", "results_numbers.txt"))
cat("[10] wrote results/results_numbers.txt (", length(out), " lines)\n", sep = "")
