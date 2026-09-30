#!/usr/bin/env Rscript
# 04_genotype_phenotype_split.R
# THE MAIN RESULT OF THIS BRANCH.
# Split the 290 accessions directly on the two lead SNPs (no haplotype clustering)
# and compare centred BLUP means for fiber and starch.
#
# WHY DIRECT AND NOT CROSSHAP: crosshap at any epsilon discards a large share of the
# panel as unassigned (see 05/06). At eps=1.5 its haplotype "C" (n=11) is a strict
# SUBSET of the true minor-allele class -- 34 accessions are ALT at both leads, but
# 20+ of them were dumped into the unassigned bin over missing data elsewhere in the
# window. The direct split uses every accession with a call at both leads.
# Outputs -> results/tables/{genotype_phenotype_per_accession.tsv, genotype_group_summary.tsv}
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("BRANCH"); OUT <- file.path(B,"results","tables")
FL <- Sys.getenv("FIBER_LEAD"); SL <- Sys.getenv("STARCH_LEAD")
V <- file.path(B,"intermediates","gene_window_raw.vcf.gz"); BCF <- Sys.getenv("BCFTOOLS")

sm <- system(paste(BCF,"query -l",V), intern=TRUE)
gt <- fread(cmd=sprintf("%s query -f '%%POS[\\t%%GT]\\n' -r 7H:%s-%s %s", BCF, FL, SL, V), header=FALSE)
pos <- as.character(gt$V1); m <- as.matrix(gt[,-1]); rownames(m)<-pos; colnames(m)<-sm
dose <- function(x) fifelse(x %in% c("1/1","1|1"),2L, fifelse(x %in% c("0/1","1/0","0|1","1|0"),1L,
                     fifelse(x %in% c("0/0","0|0"),0L, NA_integer_)))
D <- as.data.table(apply(m,1,dose)); D[, Ind := sm]

ph <- rbindlist(lapply(c("fiber","starch"), function(t){
  p <- fread(file.path(Sys.getenv("PHENO_ROOT"), paste0(t, Sys.getenv("PHENO_SUFFIX"))), header=FALSE)
  data.table(Ind=as.character(p$V2), trait=t, val=as.numeric(p$V3))}))
X <- merge(dcast(ph, Ind~trait, value.var="val"),
           D[, .(Ind, fiber_lead_dose=get(FL), starch_lead_dose=get(SL))], by="Ind")
X[, genotype_group := fifelse(is.na(fiber_lead_dose)|is.na(starch_lead_dose), "missing_call",
                       fifelse(fiber_lead_dose>0 & starch_lead_dose>0, "ALT_ALT_minor",
                        fifelse(fiber_lead_dose==0 & starch_lead_dose==0, "REF_REF_major", "discordant")))]
setorder(X, genotype_group, Ind)
fwrite(X, file.path(OUT,"genotype_phenotype_per_accession.tsv"), sep="\t")

S <- X[, .(n=.N,
  fiber_mean=mean(fiber,na.rm=TRUE), fiber_sd=sd(fiber,na.rm=TRUE), fiber_median=median(fiber,na.rm=TRUE),
  starch_mean=mean(starch,na.rm=TRUE), starch_sd=sd(starch,na.rm=TRUE), starch_median=median(starch,na.rm=TRUE)),
  by=genotype_group][order(-n)]
fwrite(S, file.path(OUT,"genotype_group_summary.tsv"), sep="\t")

# STATISTICS: deliberately the SAME test and the SAME effect sizes as step 04, so this
# branch adds no new method to the paper's M&M. The only difference from step 04 is how
# the groups are DEFINED -- observed genotype at the two lead SNPs, instead of DBSCAN
# haplotype clustering. Formulas are copied from 03_run_crosshap_pipeline.R:
#   eta_squared         = (H - k + 1) / (n - k)
#   delta_top_bottom_sd = (mean of highest group - mean of lowest) / SD of all values used
# p-values are RAW: a single pre-specified gene x 2 traits, so no BH correction is applied
# (step 04 corrects within trait ACROSS genes, which does not apply here).
TWO <- X[genotype_group %in% c("ALT_ALT_minor","REF_REF_major")]
tests <- rbindlist(lapply(c("fiber","starch"), function(tt){
  d  <- TWO[!is.na(get(tt))]
  kt <- kruskal.test(d[[tt]] ~ d$genotype_group)
  n  <- nrow(d); k <- 2L; H <- unname(kt$statistic)
  g  <- d[, .(m=mean(get(tt))), by=genotype_group][order(-m)]
  sd_all <- sd(d[[tt]])
  a <- d[genotype_group=="ALT_ALT_minor"]; b <- d[genotype_group=="REF_REF_major"]
  data.table(trait=tt, n_ALT=nrow(a), n_REF=nrow(b),
    mean_ALT=mean(a[[tt]]), mean_REF=mean(b[[tt]]),
    diff=mean(a[[tt]])-mean(b[[tt]]),
    kw_H=H, kw_df=unname(kt$parameter), kw_p_raw=kt$p.value,
    eta_squared=max(0, min(1, (H-k+1)/(n-k))),
    delta_top_bottom_sd=(g$m[1]-g$m[nrow(g)])/sd_all)}))
fwrite(tests, file.path(OUT,"genotype_group_tests.tsv"), sep="\t")
print(S); print(tests)
cat("04 OK\n")
