# 04_matched_background.R — is the Desert concentration of the lead minor alleles special,
# or simply where rare alleles sit in this panel? A genome-wide background of SNPs matched
# to each lead on minor-allele frequency and call rate.
#
# For each of the 36 leads, N_BG SNPs are drawn from the 7.1 M set with
#   |MAF - lead MAF| <= 0.005 and |called chromosomes - lead's| <= 20 (+/- 10 accessions),
# excluding everything within 2 Mb of any lead. The carrier-origin metrics of script 02 are
# computed for every background SNP, and each lead is placed in its own matched distribution.
#
# Inputs:  01_.../intermediates/morexV3_290_freq.frq, morexV3_290.{bed,bim,fam},
#          intermediates/lead_genotype_long.tsv (script 02)
# Outputs: intermediates/background_snps.tsv, background_snp_ids.txt, background_genotypes.raw
#          intermediates/accession_genomewide_rare_rate.tsv (used by script 06)
#          results/tables/T05_lead_matched_background.tsv
#          results/tables/T06_top_region_vs_background.tsv
#          results/tables/T07_background_top_region_by_MAF.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))
suppressPackageStartupMessages(library(data.table))
set.seed(SEED)
N_BG <- 500; MAF_TOL <- 0.005; NCHR_TOL <- 20; EXCL_BP <- 2e6

g   <- read_tsv(file.path(INTER, "lead_genotype_long.tsv"), show_col_types = FALSE)
acc <- g %>% distinct(IID, site, region)
leads <- g %>% filter(allele != "missing") %>%
  group_by(locus_id, trait, lead_SNP, chr, lead_bp) %>%
  summarise(n_called = n(), n_minor = sum(allele == "minor"), .groups = "drop") %>%
  mutate(nchrobs = 2 * n_called, maf = n_minor / n_called)

# ---- draw the matched background ---------------------------------------------------------
frq <- fread(FRQ_FILE, select = c("CHR", "SNP", "A1", "MAF", "NCHROBS"))[MAF > 0]
bp  <- as.numeric(sub("^.*:", "", frq$SNP))
excl <- rep(FALSE, nrow(frq))
for (i in seq_len(nrow(leads)))
  excl <- excl | (frq$CHR == leads$chr[i] & abs(bp - leads$lead_bp[i]) <= EXCL_BP)
frq <- frq[!excl]
bg <- bind_rows(lapply(seq_len(nrow(leads)), function(i) {
  cand <- frq[abs(MAF - leads$maf[i]) <= MAF_TOL & abs(NCHROBS - leads$nchrobs[i]) <= NCHR_TOL]
  stopifnot(nrow(cand) >= N_BG)
  s <- cand[sample.int(nrow(cand), N_BG)]
  tibble(lead_SNP = leads$lead_SNP[i], bg_SNP = s$SNP, bg_A1 = s$A1, n_candidates = nrow(cand))
}))
write_tsv(bg, file.path(INTER, "background_snps.tsv"))
writeLines(unique(bg$bg_SNP), file.path(INTER, "background_snp_ids.txt"))

system2(PLINK, c("--bfile", BFILE, "--extract", file.path(INTER, "background_snp_ids.txt"),
                 "--keep-allele-order", "--allow-extra-chr", "--memory", "4000",
                 "--recode", "A", "--out", file.path(INTER, "background_genotypes")),
        stdout = file.path(LOGS, "04_plink_background.log"), stderr = file.path(LOGS, "04_plink_background.log"))
invisible(file.remove(list.files(INTER, "^background_genotypes\\.(nosex|log)$", full.names = TRUE)))

raw <- fread(file.path(INTER, "background_genotypes.raw"))
snp_cols <- names(raw)[-(1:6)]; col_snp <- sub("_[ACGT]+$", "", snp_cols)
stopifnot(all(sub("^.*_", "", snp_cols) == frq$A1[match(col_snp, frq$SNP)]))   # counted allele = minor
setnames(raw, snp_cols, col_snp)
G <- as.matrix(raw[, ..col_snp])
stopifnot(sum(G == 1, na.rm = TRUE) == 0)                                          # no heterozygotes
reg  <- as.character(acc$region[match(raw$IID, acc$IID)]); site <- acc$site[match(raw$IID, acc$IID)]
stopifnot(!anyNA(reg))

# ---- carrier-origin metrics per background SNP ------------------------------------------------
metrics <- function(x) {
  minor <- !is.na(x) & x == 2
  rt <- table(factor(reg[minor], levels = REGION_ORDER))
  c(desert_share = mean(reg[minor] == "Desert"), top_region_share = max(rt) / sum(minor),
    top_region = unname(which.max(rt)), n_sites = length(unique(site[minor])))
}
bgm <- as.data.frame(t(apply(G, 2, metrics)))
bgm$bg_SNP <- rownames(bgm); bgm$top_region <- REGION_ORDER[bgm$top_region]
bgm <- bg %>% left_join(bgm, by = "bg_SNP")

lead_m <- g %>% filter(allele == "minor") %>% group_by(lead_SNP) %>%
  summarise(lead_desert_share = mean(region == "Desert"), lead_n_sites = n_distinct(site),
            lead_top_region = names(which.max(table(factor(region, levels = REGION_ORDER)))),
            lead_top_region_share = max(table(factor(region, levels = REGION_ORDER))) / n(),
            .groups = "drop")

t05 <- bgm %>% left_join(lead_m, by = "lead_SNP") %>% group_by(lead_SNP) %>%
  summarise(n_bg = n(), n_candidates = first(n_candidates),
            lead_desert_pct = round(100 * first(lead_desert_share), 1),
            bg_desert_pct_median = round(100 * median(desert_share), 1),
            bg_desert_pct_q95 = round(100 * quantile(desert_share, 0.95), 1),
            p_bg_desert = round((1 + sum(desert_share >= first(lead_desert_share) - 1e-12)) / (n() + 1), 4),
            lead_top_region = first(lead_top_region),
            bg_frac_top_region_desert = round(mean(top_region == "Desert"), 3),
            bg_frac_top_region_north = round(mean(top_region == "North"), 3),
            lead_n_sites = first(lead_n_sites), bg_n_sites_median = median(n_sites),
            p_bg_fewer_sites = round((1 + sum(n_sites <= first(lead_n_sites))) / (n() + 1), 4),
            lead_top_region_pct = round(100 * first(lead_top_region_share), 1),
            bg_top_region_pct_median = round(100 * median(top_region_share), 1), .groups = "drop") %>%
  left_join(select(leads, locus_id, trait, lead_SNP, maf), by = "lead_SNP") %>%
  mutate(maf = round(maf, 4)) %>% relocate(trait, locus_id, lead_SNP, maf) %>%
  arrange(factor(trait, levels = TRAIT_ORDER), match(lead_SNP, leads$lead_SNP))
write_tab(t05, file.path(TAB, "T05_lead_matched_background.tsv"))

# ---- T06: leads per top region, observed vs expected from each lead's own background -----------
# p by simulation of the Poisson-binomial (sum of per-lead Bernoulli draws). The leads of a trait
# share carriers (T09), so they are not independent draws and these p-values are optimistic.
bg_frac <- bgm %>% count(lead_SNP, top_region) %>% group_by(lead_SNP) %>% mutate(f = n / sum(n)) %>%
  ungroup() %>% complete(lead_SNP, top_region = REGION_ORDER, fill = list(n = 0, f = 0))
lt <- lead_m %>% left_join(distinct(leads, lead_SNP, trait), by = "lead_SNP")
t06 <- bind_rows(lapply(c(TRAIT_ORDER, "fiber+betaglucan", "all_but_protein"), function(t) {
  sel <- switch(t, "fiber+betaglucan" = lt$trait %in% c("fiber", "betaglucan"),
                "all_but_protein" = lt$trait != "protein", lt$trait == t)
  l <- lt[sel, ]
  bind_rows(lapply(REGION_ORDER, function(r) {
    f <- bg_frac %>% filter(lead_SNP %in% l$lead_SNP, top_region == r) %>% pull(f)
    obs <- sum(l$lead_top_region == r)
    sims <- colSums(matrix(runif(length(f) * 1e5) < f, nrow = length(f)))
    tibble(leads = t, n_leads = nrow(l), top_region = r, observed = obs, expected = round(sum(f), 2),
           p_more = signif((1 + sum(sims >= obs)) / (1e5 + 1), 3),
           p_fewer = signif((1 + sum(sims <= obs)) / (1e5 + 1), 3))
  }))
}))
write_tab(t06, file.path(TAB, "T06_top_region_vs_background.tsv"))

# ---- T07: pooled background, which region holds most carriers, by MAF class --------------------
t07 <- bgm %>% distinct(bg_SNP, .keep_all = TRUE) %>%
  left_join(select(leads, lead_SNP, maf), by = "lead_SNP") %>%
  mutate(maf_class = ifelse(maf < 0.10, "MAF < 0.10", "MAF >= 0.10")) %>%
  count(maf_class, top_region) %>% group_by(maf_class) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>% arrange(maf_class, factor(top_region, levels = REGION_ORDER))
write_tab(t07, file.path(TAB, "T07_background_top_region_by_MAF.tsv"))

# ---- genome-wide rare-allele load per accession (background SNPs with MAF < 0.10) --------------
R <- G[, colMeans(G, na.rm = TRUE) / 2 < 0.10, drop = FALSE]
rr <- tibble(IID = raw$IID, gw_rare_minor = rowSums(R == 2, na.rm = TRUE), gw_called = rowSums(!is.na(R))) %>%
  mutate(gw_rare_minor_rate = round(gw_rare_minor / gw_called, 4),
         gw_missing_rate = round(1 - gw_called / ncol(R), 4),
         gw_rare_rate_rank = rank(-gw_rare_minor_rate, ties.method = "min"))
write_tsv(rr, file.path(INTER, "accession_genomewide_rare_rate.tsv"))
cat("[04] background SNPs:", length(unique(bg$bg_SNP)), "| rare background SNPs:", ncol(R), "\n")
