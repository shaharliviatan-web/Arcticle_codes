# 05_lead_within_site.R — do the minor-allele carriers differ from their own site-mates?
# A descriptive decomposition of the raw carrier difference; no association is re-tested
# and no covariate is added.
#
# Per lead (phenotype = the BLUP the GWAS used, in trait SD):
#   observed   = mean(minor-allele carriers) - mean(major-allele carriers)
#   expected   = what the carriers' origin alone predicts: each carrier is replaced by the mean
#                of the major-allele accessions of its own site (its non-carrier site-mates) or
#                of its own region; a site with no major-allele accession falls back to its region.
#   frac_expected_site          = expected (site) / observed
#   pct_carriers_beyond_sitemates = carriers whose value exceeds their site-mates' mean, in the
#                                   direction of the minor allele's effect
#   n_minor_with_sitemates      = carriers whose site holds at least one major-allele accession
#                                 (the others sit in sites fixed for the minor allele, where
#                                 allele and population cannot be separated)
# and the dependence on the panel's most extreme accessions, ranked on the phenotype alone in
# the direction of the minor allele's effect: carriers among the top 10, and the observed
# difference recomputed without the top 2.
#
# Also: raw_sign_agrees_with_beta — the lead_beta is estimated with 3 PCs + kinship; where the
# raw carrier difference has the other sign, the association exists only after that correction.
#
# Inputs:  intermediates/lead_genotype_long.tsv (script 02)
# Outputs: results/tables/T08_lead_within_site.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))

g <- read_tsv(file.path(INTER, "lead_genotype_long.tsv"), show_col_types = FALSE)
lead_order <- read_tsv(LOCI_FILE, show_col_types = FALSE)$lead_SNP

t08 <- bind_rows(lapply(split(g, g$lead_SNP), function(d) {
  sd_y <- sd(d$y); dir <- sign(d$lead_beta[1])
  called <- d %>% filter(allele != "missing")
  M <- called %>% filter(allele == "minor"); J <- called %>% filter(allele == "major")
  mu_J <- mean(J$y)
  bg_reg  <- J %>% group_by(region) %>% summarise(bg_region = mean(y), .groups = "drop")
  bg_site <- J %>% group_by(site) %>% summarise(bg_site = mean(y), n_sitemates = n(), .groups = "drop")
  Mb <- M %>% left_join(bg_reg, by = "region") %>% left_join(bg_site, by = "site") %>%
    mutate(n_sitemates = coalesce(n_sitemates, 0L), bg_site_fb = coalesce(bg_site, bg_region))
  obs <- mean(M$y) - mu_J
  exp_reg <- mean(Mb$bg_region) - mu_J; exp_site <- mean(Mb$bg_site_fb) - mu_J
  within <- Mb$y - Mb$bg_site_fb
  ext <- d$IID[order(-dir * d$y)]; top10 <- ext[1:10]; top2 <- ext[1:2]
  obs_wo <- mean(M$y[!M$IID %in% top2]) - mean(J$y[!J$IID %in% top2])
  tibble(trait = d$trait[1], locus_id = d$locus_id[1], lead_SNP = d$lead_SNP[1],
         gene = unname(ifelse(d$lead_SNP[1] %in% names(GENE_TAG), GENE_TAG[d$lead_SNP[1]], "")),
         minor_effect = ifelse(dir > 0, "raises", "lowers"), n_minor = nrow(M),
         obs_diff_SD = round(obs / sd_y, 3), raw_sign_agrees_with_beta = sign(obs) == dir,
         exp_diff_region_SD = round(exp_reg / sd_y, 3), exp_diff_site_SD = round(exp_site / sd_y, 3),
         frac_expected_region = round(exp_reg / obs, 3), frac_expected_site = round(exp_site / obs, 3),
         n_minor_with_sitemates = sum(Mb$n_sitemates > 0),
         pct_carriers_beyond_sitemates = round(100 * mean(sign(within) == dir), 1),
         top2_extreme = paste(top2, collapse = ","), n_top2_extreme_carriers = sum(M$IID %in% top2),
         n_top10_extreme_carriers = sum(M$IID %in% top10),
         obs_diff_SD_without_top2 = round(obs_wo / sd_y, 3),
         pct_diff_kept_without_top2 = round(100 * obs_wo / obs, 1))
})) %>% arrange(factor(trait, levels = TRAIT_ORDER), match(lead_SNP, lead_order))
write_tab(t08, file.path(TAB, "T08_lead_within_site.tsv"))
cat("[05] leads:", nrow(t08), "| raw sign disagrees with beta:",
    paste(t08$locus_id[!t08$raw_sign_agrees_with_beta], collapse = ", "), "\n")
