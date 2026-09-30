# 02_lead_carrier_origin.R — regional origin of the minor- and major-allele carriers of the
# 36 GWAS lead SNPs, and the Desert and six-region enrichment tests.
#
# Unit: the lead SNP of each locus. At each lead the 290 accessions are split into
# minor-allele carriers (homozygous A1; the call set has no heterozygotes), major-allele
# carriers (homozygous A2) and missing calls (excluded from the tests).
#
# Tests (site permutation, 00_config.R: perm_tests):
#   Desert  — is the Desert share of the minor-allele carriers different from that of the
#             major-allele carriers? two-sided.
#   regions — are minor- and major-allele carriers distributed differently over the six
#             regions? chi-square statistic, upper tail.
# BH within trait for each test. The accession-level Fisher odds ratio (Desert, minor vs
# major) is reported as an effect size only: it treats accessions of one site as
# independent, which they are not.
#
# Inputs:  intermediates/lead_genotypes.raw (script 01), Table_loci_master.tsv, region table, BLUPs
# Outputs: intermediates/lead_genotype_long.tsv (accession x lead)
#          results/tables/T01_lead_carriers_by_region.tsv
#          results/tables/T02_lead_enrichment_tests.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))

loci <- read_tsv(LOCI_FILE, show_col_types = FALSE)
stopifnot(nrow(loci) == 36)
acc   <- load_accessions()
pheno <- load_pheno()

raw <- read.table(file.path(INTER, "lead_genotypes.raw"), header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
gcols <- names(raw)[-(1:6)]
csnp <- sub("_[ACGT]+$", "", gcols); ca1 <- sub("^.*_", "", gcols)
# the counted allele of --recode A must be the table's lead_A1 (the minor allele)
stopifnot(setequal(csnp, loci$lead_SNP), all(ca1 == loci$lead_A1[match(csnp, loci$lead_SNP)]))
names(raw)[-(1:6)] <- csnp
stopifnot(setequal(raw$IID, acc$IID))

g <- raw %>% select(IID, all_of(loci$lead_SNP)) %>%
  pivot_longer(-IID, names_to = "lead_SNP", values_to = "dosage")
stopifnot(sum(g$dosage == 1, na.rm = TRUE) == 0)          # no heterozygous calls
g <- g %>%
  mutate(allele = factor(case_when(is.na(dosage) ~ "missing", dosage == 2 ~ "minor", dosage == 0 ~ "major"),
                         levels = c("minor", "major", "missing"))) %>%
  left_join(select(loci, locus_id, trait, chr, lead_SNP, lead_bp, lead_MAF, lead_beta, lead_neg_log10_p),
            by = "lead_SNP") %>%
  left_join(select(acc, IID, short_Tag, site, location, region), by = "IID") %>%
  left_join(pheno, by = c("trait", "IID"))
stopifnot(!anyNA(g$region), !anyNA(g$y))
write_tsv(g, file.path(INTER, "lead_genotype_long.tsv"))

# ---- T01: locus x region counts ------------------------------------------------------------
t01 <- g %>% count(trait, locus_id, lead_SNP, region, allele, .drop = FALSE) %>%
  filter(!is.na(locus_id)) %>%
  pivot_wider(names_from = allele, values_from = n, values_fill = 0) %>%
  rename(n_minor = minor, n_major = major, n_missing = missing) %>%
  group_by(locus_id) %>%
  mutate(pct_minor_within_region = round(100 * n_minor / (n_minor + n_major), 1),
         pct_of_minor_carriers = round(100 * n_minor / sum(n_minor), 1),
         pct_of_major_carriers = round(100 * n_major / sum(n_major), 1)) %>% ungroup() %>%
  arrange(factor(trait, levels = TRAIT_ORDER), match(lead_SNP, loci$lead_SNP), region)
write_tab(t01, file.path(TAB, "T01_lead_carriers_by_region.tsv"))

# ---- T02: one row per lead, with the tests --------------------------------------------------
SP <- site_perm_matrix(acc)
t02 <- bind_rows(lapply(loci$lead_SNP, function(snp) {
  d  <- g %>% filter(lead_SNP == snp)
  cd <- d %>% filter(allele != "missing")
  M  <- cd %>% filter(allele == "minor")
  focal <- cd$allele == "minor"
  pt <- perm_tests(cd$site, cd$region, focal, factor(focal), SP)
  rc <- region_counts(M$region)
  st <- sort(table(M$site), decreasing = TRUE)
  tibble(
    trait = d$trait[1], locus_id = d$locus_id[1], lead_SNP = snp,
    gene = unname(ifelse(snp %in% names(GENE_TAG), GENE_TAG[snp], "")),
    lead_neg_log10_p = d$lead_neg_log10_p[1], MAF = round(nrow(M) / nrow(cd), 4),
    minor_effect = ifelse(d$lead_beta[1] > 0, "raises", "lowers"),
    n_minor = nrow(M), n_major = sum(!focal), n_missing = sum(d$allele == "missing"),
    minor_North = rc[1], minor_Coast = rc[2], minor_Desert = rc[3],
    minor_HZ1 = rc[4], minor_HZ2 = rc[5], minor_HZ3 = rc[6],
    n_regions_minor = sum(rc > 0), n_sites_minor = length(st),
    n_desert_sites_minor = n_distinct(M$site[M$region == "Desert"]),
    top_site = d$location[match(names(st)[1], d$site)], top_site_n = as.integer(st[1]),
    desert_pct_minor = round(100 * mean(M$region == "Desert"), 1),
    desert_pct_major = round(100 * mean(cd$region[!focal] == "Desert"), 1),
    desert_OR = round(desert_or(focal, cd$region), 2),
    desert_direction = ifelse(pt$desert_diff > 0, "enriched", ifelse(pt$desert_diff < 0, "depleted", "equal")),
    p_desert = signif(pt$p_desert, 3), p_regions = signif(pt$p_regions, 3))
}))
t02 <- t02 %>% group_by(trait) %>%
  mutate(q_desert = signif(p.adjust(p_desert, "BH"), 3),
         q_regions = signif(p.adjust(p_regions, "BH"), 3)) %>% ungroup() %>%
  relocate(q_desert, .after = p_desert) %>%
  arrange(factor(trait, levels = TRAIT_ORDER), match(lead_SNP, loci$lead_SNP))
write_tab(t02, file.path(TAB, "T02_lead_enrichment_tests.tsv"))
cat("[02] leads:", nrow(t02), "| Desert q<=0.05:", sum(t02$q_desert <= 0.05),
    "| regions q<=0.05:", sum(t02$q_regions <= 0.05), "\n")
