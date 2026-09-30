# 06_lead_carrier_sharing.R — are the minor-allele carriers of different loci the same
# accessions?
#
#   T09 — every pair of loci of one trait: minor-allele carriers shared among the accessions
#         called at both, the number expected by chance, and a hypergeometric upper-tail p.
#         Pairs on different chromosomes cannot share carriers through linkage.
#   T10 — every accession: the number of lead minor alleles it carries per trait, its BLUPs, and
#         its genome-wide rare-allele load (background SNPs with MAF < 0.10, script 04), so an
#         accession that recurs at the leads can be compared with its load elsewhere.
#
# Inputs:  intermediates/lead_genotype_long.tsv (02), accession_genomewide_rare_rate.tsv (04)
# Outputs: results/tables/T09_lead_carrier_sharing.tsv, T10_accession_minor_allele_load.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))

g  <- read_tsv(file.path(INTER, "lead_genotype_long.tsv"), show_col_types = FALSE)
L  <- g %>% distinct(trait, locus_id, lead_SNP, chr, lead_bp)
pheno <- load_pheno()

t09 <- bind_rows(lapply(TRAIT_ORDER, function(t) {
  l <- L %>% filter(trait == t); if (nrow(l) < 2) return(NULL)
  cmb <- combn(l$lead_SNP, 2)
  bind_rows(lapply(seq_len(ncol(cmb)), function(k) {
    a <- g %>% filter(lead_SNP == cmb[1, k]); b <- g %>% filter(lead_SNP == cmb[2, k])
    both <- intersect(a$IID[a$allele != "missing"], b$IID[b$allele != "missing"])
    ca <- intersect(a$IID[a$allele == "minor"], both); cb <- intersect(b$IID[b$allele == "minor"], both)
    sh <- length(intersect(ca, cb)); la <- l[l$lead_SNP == cmb[1, k], ]; lb <- l[l$lead_SNP == cmb[2, k], ]
    tibble(trait = t, locus_a = la$locus_id, lead_a = la$lead_SNP, locus_b = lb$locus_id, lead_b = lb$lead_SNP,
           same_chr = la$chr == lb$chr,
           distance_Mb = ifelse(la$chr == lb$chr, round(abs(la$lead_bp - lb$lead_bp) / 1e6, 3), NA),
           n_called_both = length(both), n_minor_a = length(ca), n_minor_b = length(cb),
           n_shared = sh, expected_shared = round(length(ca) * length(cb) / length(both), 2),
           jaccard = round(sh / length(union(ca, cb)), 3),
           p_hyper = signif(phyper(sh - 1, length(ca), length(both) - length(ca), length(cb), lower.tail = FALSE), 3))
  }))
}))
write_tab(t09, file.path(TAB, "T09_lead_carrier_sharing.tsv"))

rr <- read_tsv(file.path(INTER, "accession_genomewide_rare_rate.tsv"), show_col_types = FALSE)
t10 <- g %>% group_by(IID, short_Tag, site, location, region, trait) %>%
  summarise(n = sum(allele == "minor"), .groups = "drop") %>%
  pivot_wider(names_from = trait, values_from = n, names_prefix = "minor_") %>%
  mutate(minor_total = minor_betaglucan + minor_fiber + minor_starch + minor_protein) %>%
  left_join(pheno %>% pivot_wider(names_from = trait, values_from = y, names_prefix = "BLUP_"), by = "IID") %>%
  left_join(rr, by = "IID") %>%
  group_by() %>% mutate(fiber_rank_high = rank(-BLUP_fiber), betaglucan_rank_high = rank(-BLUP_betaglucan),
                        starch_rank_low = rank(BLUP_starch)) %>% ungroup() %>%
  arrange(desc(minor_total), region, site)
write_tab(t10, file.path(TAB, "T10_accession_minor_allele_load.tsv"))
x <- t09 %>% filter(!same_chr)
cat("[06] pairs on different chromosomes with p_hyper < 0.001:",
    paste(x %>% group_by(trait) %>% summarise(s = paste0(trait[1], " ", sum(p_hyper < 0.001), "/", n())) %>% pull(s),
          collapse = "; "), "\n")
