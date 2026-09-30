# 07_kinship_by_region.R — how related are the accessions within and between regions?
# Backs the Discussion statement that the desert populations are closely related to one
# another (TODO 405). Descriptive, from the kinship matrix the GWAS used (EMMAX aIBS).
#
# Every pair of the 290 accessions is classed by the region of each member and by whether
# the two come from the same site. Reported per region pair: mean, median and number of
# pairs, for same-site pairs and for pairs from different sites. The different-site mean
# within a region is the relatedness between its populations.
#
# Inputs:  01_.../intermediates/morexV3_kinship.aIBS.kinf (rows/cols in .fam order)
# Outputs: results/tables/T11_kinship_by_region.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))

acc <- load_accessions()
K   <- as.matrix(read.table(KIN_FILE))
ids <- read.table(paste0(BFILE, ".fam"))$V1
stopifnot(nrow(K) == 290, ncol(K) == 290, setequal(ids, acc$IID))
reg  <- as.character(acc$region[match(ids, acc$IID)]); site <- acc$site[match(ids, acc$IID)]

ut <- which(upper.tri(K), arr.ind = TRUE)
pairs <- tibble(k = K[ut], r1 = reg[ut[, 1]], r2 = reg[ut[, 2]], same_site = site[ut[, 1]] == site[ut[, 2]]) %>%
  mutate(a = pmin(match(r1, REGION_ORDER), match(r2, REGION_ORDER)),
         b = pmax(match(r1, REGION_ORDER), match(r2, REGION_ORDER)),
         region_1 = REGION_ORDER[a], region_2 = REGION_ORDER[b])

t11 <- pairs %>%
  mutate(pair_type = ifelse(same_site, "same site", ifelse(region_1 == region_2, "different sites, same region",
                                                            "different regions"))) %>%
  group_by(region_1, region_2, pair_type) %>%
  summarise(n_pairs = n(), mean_kinship = round(mean(k), 4), median_kinship = round(median(k), 4), .groups = "drop") %>%
  bind_rows(pairs %>% summarise(region_1 = "all", region_2 = "all", pair_type = "all pairs", n_pairs = n(),
                                mean_kinship = round(mean(k), 4), median_kinship = round(median(k), 4))) %>%
  arrange(match(region_1, c(REGION_ORDER, "all")), match(region_2, c(REGION_ORDER, "all")), pair_type)
write_tab(t11, file.path(TAB, "T11_kinship_by_region.tsv"))
w <- t11 %>% filter(pair_type == "different sites, same region")
cat("[07] mean aIBS between sites of the same region:", paste(REGION_SHORT[match(w$region_1, REGION_ORDER)], w$mean_kinship, collapse = "; "), "\n")
