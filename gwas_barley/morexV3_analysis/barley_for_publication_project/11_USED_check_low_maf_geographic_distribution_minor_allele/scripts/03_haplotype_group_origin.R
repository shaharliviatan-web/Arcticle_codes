# 03_haplotype_group_origin.R — do the haplotype groups of the four presented genes come
# from different regions and sites?
#
# Genes and groups (read-only, never re-run):
#   GPAT6, GH17, PHT4;3 — step 04 crosshap cache, run loci_LDspan_eps06_V4 (MGmin 2, eps 0.6),
#                         HapObject Indfile (the groups of Fig. 4 and Table 2)
#   GDSL (7H)           — 03_01 branch, mgmin3_haplotype_assignment.tsv (MGmin 3, eps 0.9;
#                         the groups of Fig. 5; identical for fiber and starch)
# Group "0" = unassigned: kept in the counts, left out of every test (as in the haplotype test).
#
# Tests (site permutation, 00_config.R: perm_tests), on the assigned accessions:
#   gene level  — do the groups of a gene differ in their regional composition?
#                 groups x six regions, chi-square statistic, upper tail; BH over the 4 genes.
#   group level — each group against the other assigned accessions of the same gene:
#                 Desert share (two-sided) and six regions (upper tail); BH within gene.
#
# Outputs: intermediates/haplotype_group_long.tsv (gene x accession)
#          results/tables/T03_haplotype_group_members_by_region.tsv
#          results/tables/T04_haplotype_group_enrichment_tests.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))

acc   <- load_accessions()
pheno <- load_pheno()

read_step04 <- function(gene_id, trait) {
  ho  <- readRDS(file.path(STEP04_CACHE, trait, gene_id, "MGmin_2", "HapObject.rds"))
  ind <- ho$HapObject[["Haplotypes_MGmin2_E0.6"]]$Indfile
  # Indfile Pheno must be the same BLUPs as the GWAS phenotype files
  chk <- pheno %>% filter(trait == !!trait)
  stopifnot(nrow(ind) == 290, isTRUE(all.equal(ind$Pheno, chk$y[match(ind$Ind, chk$IID)])))
  tibble(IID = ind$Ind, hap = as.character(ind$hap))
}
read_gdsl <- function() {
  a  <- read_tsv(GDSL_ASSIGN, show_col_types = FALSE)
  fb <- a %>% filter(trait == "fiber"); sc <- a %>% filter(trait == "starch")
  stopifnot(nrow(fb) == 290, all(fb$MGmin == 3), all(fb$epsilon == 0.9),
            identical(fb$hap[order(fb$Ind)], sc$hap[order(sc$Ind)]))   # one grouping for both traits
  tibble(IID = fb$Ind, hap = as.character(fb$hap))
}

h <- bind_rows(lapply(seq_len(nrow(GENES)), function(i) {
  x <- if (GENES$source[i] == "step04") read_step04(GENES$gene_id[i], GENES$trait[i]) else read_gdsl()
  x %>% mutate(gene = GENES$gene[i], gene_id = GENES$gene_id[i], trait = GENES$trait[i])
})) %>% left_join(select(acc, IID, site, location, region), by = "IID")
stopifnot(!anyNA(h$region))
h <- h %>%
  left_join(pheno %>% filter(trait == "fiber")  %>% select(IID, BLUP_fiber = y),  by = "IID") %>%
  left_join(pheno %>% filter(trait == "starch") %>% select(IID, BLUP_starch = y), by = "IID")
write_tsv(h, file.path(INTER, "haplotype_group_long.tsv"))

# ---- T03: gene x group x region counts (unassigned = "0") ------------------------------------
t03 <- h %>% mutate(group = ifelse(hap == "0", "unassigned", hap), region = factor(region, levels = REGION_ORDER)) %>%
  count(gene, gene_id, trait, group, region, .drop = FALSE) %>%
  group_by(gene, group) %>% mutate(pct_of_group = round(100 * n / sum(n), 1)) %>% ungroup() %>%
  arrange(match(gene, GENES$gene), group, region)
write_tab(t03, file.path(TAB, "T03_haplotype_group_members_by_region.tsv"))

# ---- T04: group-level and gene-level tests -----------------------------------------------------
SP <- site_perm_matrix(acc)
t04 <- bind_rows(lapply(GENES$gene, function(gn) {
  d  <- h %>% filter(gene == gn, hap != "0")
  gr <- sort(unique(d$hap))
  gene_p <- perm_tests(d$site, d$region, d$hap == gr[1], factor(d$hap), SP)$p_regions
  bind_rows(lapply(gr, function(k) {
    focal <- d$hap == k; M <- d[focal, ]
    pt <- perm_tests(d$site, d$region, focal, factor(focal), SP)
    rc <- region_counts(M$region); st <- sort(table(M$site), decreasing = TRUE)
    tibble(
      gene = gn, gene_id = d$gene_id[1], trait = d$trait[1], group = k, n = nrow(M),
      n_assigned_gene = nrow(d), n_groups_gene = length(gr),
      mean_BLUP_fiber  = if (grepl("fiber",  d$trait[1])) round(mean(M$BLUP_fiber), 3)  else NA_real_,
      mean_BLUP_starch = if (grepl("starch", d$trait[1])) round(mean(M$BLUP_starch), 3) else NA_real_,
      North = rc[1], Coast = rc[2], Desert = rc[3], HZ1 = rc[4], HZ2 = rc[5], HZ3 = rc[6],
      n_regions = sum(rc > 0), n_sites = length(st),
      n_desert_sites = n_distinct(M$site[M$region == "Desert"]),
      top_site = M$location[match(names(st)[1], M$site)], top_site_n = as.integer(st[1]),
      desert_pct_group = round(100 * mean(M$region == "Desert"), 1),
      desert_pct_other_groups = round(100 * mean(d$region[!focal] == "Desert"), 1),
      desert_OR = round(desert_or(focal, d$region), 2),
      desert_direction = ifelse(pt$desert_diff > 0, "enriched", ifelse(pt$desert_diff < 0, "depleted", "equal")),
      p_desert = signif(pt$p_desert, 3), p_regions_group = signif(pt$p_regions, 3),
      p_regions_gene = signif(gene_p, 3))
  }))
}))
t04 <- t04 %>% group_by(gene) %>%
  mutate(q_desert = signif(p.adjust(p_desert, "BH"), 3),
         q_regions_group = signif(p.adjust(p_regions_group, "BH"), 3)) %>% ungroup()
gq <- t04 %>% distinct(gene, p_regions_gene) %>% mutate(q_regions_gene = signif(p.adjust(p_regions_gene, "BH"), 3))
t04 <- t04 %>% left_join(select(gq, gene, q_regions_gene), by = "gene") %>%
  relocate(q_desert, .after = p_desert) %>% relocate(q_regions_group, .after = p_regions_group) %>%
  arrange(match(gene, GENES$gene), group)
write_tab(t04, file.path(TAB, "T04_haplotype_group_enrichment_tests.tsv"))
cat("[03] genes:", nrow(GENES), "| groups:", nrow(t04), "| gene-level q:",
    paste(gq$gene, gq$q_regions_gene, collapse = "; "), "\n")
