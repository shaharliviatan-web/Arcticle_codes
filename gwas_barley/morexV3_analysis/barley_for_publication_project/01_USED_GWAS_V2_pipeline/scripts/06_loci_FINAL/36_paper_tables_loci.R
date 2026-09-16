#!/usr/bin/env Rscript
# 36_paper_tables_loci.R
# Build the publication-ready table set for the FINAL locus definition, so the
# paper can be written from tables alone without re-deriving anything.
#
# Locus definition (values are read at runtime from locus_definition_params.tsv):
#   leads      PLINK --clump --clump-p1 9.008e-07  (Bonferroni alpha=0.10 over the
#              111,017 LD-pruned SNPs) -- only genome-wide-significant SNPs lead
#   members    --clump-p2 1  => membership on LD alone, NO p-value condition,
#              because a causal variant need not itself be significant
#   LD         --clump-r2 0.5
#   reach      --clump-kb 2000  (+/- 2 Mb, so 4 Mb maximum span)
#   contiguity severed at the first gap > 60 kb between consecutive members,
#              walking outward from the lead in each direction
#
# Outputs -> 02_loci_FINAL/tables/
#   Table_loci_master.tsv        one row per locus: coordinates, span, members,
#                                lead stats (beta/SE/p/MAF), severing diagnostics
#   Table_loci_members_full.tsv  every member SNP with its position and p-value
#   Table_loci_per_trait.tsv     per-trait locus counts and span distribution
#   Table_analysis_parameters.tsv every constant, one row each
#   Table_extra_peaks.tsv        the 2 sub-threshold 3H protein peaks
#   Table_loci_for_gene_search.tsv  HANDOFF: the minimal coordinate table consumed
#                                by 03_USED_candidate_genes_around_leading_snps --
#                                one row per search interval, loci and extra peaks
#                                together, with an explicit include_in_gene_search
#                                flag so the two sub-threshold peaks can be
#                                switched off without editing the file.
# Created 2026-09-08.
Sys.setenv(TMPDIR="/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))
PIPE <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
SET  <- file.path(PIPE,"results","00_FINAL_BLUP_3PC")
LOC  <- file.path(SET,"02_loci_FINAL"); TAB <- file.path(LOC,"tables")
INTER<- file.path(PIPE,"intermediates")

N_PRUNED <- length(readLines(file.path(INTER,"morexV3_pruned_for_covs.prune.in")))
ALPHA <- 0.10; BONF <- -log10(ALPHA/N_PRUNED)
PAR <- fread(file.path(TAB,"locus_definition_params.tsv"))
gp  <- function(k) PAR[parameter==k, value][1]
R2 <- as.numeric(gp("clump_r2")); KB <- as.numeric(gp("clump_kb"))
MAX_GAP_KB <- as.numeric(gp("gap_rule_kb")); P1 <- ALPHA/N_PRUNED

S  <- fread(file.path(TAB,"loci_summary.tsv"))
M  <- fread(file.path(TAB,"loci_members.tsv"))
# The two sub-threshold protein peaks were dropped from the analysis set on
# 2026-09-09: protein keeps only its 2 Bonferroni-significant loci. The blocks below
# stay conditional so the tables build whether or not extra-peak files are present.
fe <- file.path(TAB,"loci_summary_extra.tsv"); fm <- file.path(TAB,"loci_members_extra.tsv")
HAS_EXTRA <- file.exists(fe) && file.exists(fm)
SE <- if (HAS_EXTRA) fread(fe) else NULL
ME <- if (HAS_EXTRA) fread(fm) else NULL
maf <- fread(file.path(INTER,"morexV3_290_freq.frq"))[, .(SNP, MAF)]; setkey(maf, SNP)
bim <- fread(file.path(INTER,"morexV3_290.bim"), header=FALSE,
             col.names=c("chr","SNP","cm","bp","A1","A2"))[, .(SNP, A1, A2)]; setkey(bim, SNP)

ps <- rbindlist(lapply(c("betaglucan","fiber","protein","starch"), function(tr)
  fread(file.path(PIPE,"results","emmax_ps",sprintf("morexV3__%s__BLUP__pc3.ps", tr)),
        header=FALSE, col.names=c("SNP","beta","SE","P"))[, trait := tr]))
setkey(ps, trait, SNP)

add_lead <- function(D) {
  D <- copy(D)
  D[, `:=`(lead_MAF   = maf[lead_SNP, MAF],
           lead_A1    = bim[lead_SNP, A1],
           lead_A2    = bim[lead_SNP, A2])]
  st <- ps[.(D$trait, D$lead_SNP), .(beta, SE, P)]
  D[, `:=`(lead_beta = round(st$beta,6), lead_SE = round(st$SE,6),
           lead_p = signif(st$P,4), lead_neg_log10_p = round(-log10(st$P),4))]
  D[]
}
S <- add_lead(S); if (HAS_EXTRA) SE <- add_lead(SE)

master <- S[, .(locus_id, trait, chr, lead_SNP, lead_bp,
                locus_start = start, locus_end = end, span_kb = round(span_kb,3),
                n_member_SNPs = n_members,
                lead_A1, lead_A2, lead_MAF = round(lead_MAF,4),
                lead_beta, lead_SE, lead_p, lead_neg_log10_p,
                max_internal_gap_kb = round(max_internal_gap_kb,3),
                clump_span_before_severing_kb = round(raw_clump_span_kb,3),
                n_clump_members_before_severing = n_clump_members,
                n_members_severed = n_severed)]
setorder(master, trait, chr, lead_bp)
fwrite(master, file.path(TAB,"Table_loci_master.tsv"), sep="\t")

MM <- if (HAS_EXTRA) rbind(M, ME, fill=TRUE) else M
key <- if (HAS_EXTRA) rbind(S[, .(locus_id, trait, lead_SNP)], SE[, .(locus_id, trait, lead_SNP)]) else S[, .(locus_id, trait, lead_SNP)]
MM <- merge(MM, key, by=c("locus_id","trait"), all.x=TRUE)
MM[, `:=`(MAF = round(maf[SNP, MAF],4))]
st <- ps[.(MM$trait, MM$SNP), .(P)]
MM[, `:=`(p_value = signif(st$P,4), neg_log10_p = round(-log10(st$P),4))]
setorder(MM, trait, locus_id, BP)
fwrite(MM[, .(locus_id, trait, lead_SNP, SNP, position_bp = BP, is_index,
              MAF, p_value, neg_log10_p)],
       file.path(TAB,"Table_loci_members_full.tsv"), sep="\t")

per <- S[, .(n_loci = .N,
             n_lead_only = sum(n_members == 1L),
             total_member_SNPs = sum(n_members),
             span_min_kb = round(min(span_kb),1), span_median_kb = round(median(span_kb),1),
             span_max_kb = round(max(span_kb),1),
             largest_internal_gap_kb = round(max(max_internal_gap_kb),1),
             total_severed = sum(n_severed)), by = trait]
fwrite(per, file.path(TAB,"Table_loci_per_trait.tsv"), sep="\t")

if (HAS_EXTRA) fwrite(SE[, .(peak_id = locus_id, trait, chr, lead_bp, lead_A1, lead_A2,
              lead_MAF = round(lead_MAF,4), lead_beta, lead_SE, lead_p, lead_neg_log10_p,
              Bonferroni_threshold_neg_log10_p = round(BONF,4),
              locus_start = start, locus_end = end, span_kb = round(span_kb,3),
              n_member_SNPs = n_members, max_internal_gap_kb = round(max_internal_gap_kb,3),
              n_members_severed = n_severed,
              note = "sub-threshold; below the lead threshold so --clump cannot use it as an index SNP; anchored directly under identical rules")],
       file.path(TAB,"Table_extra_peaks.tsv"), sep="\t")

fwrite(data.table(
  parameter = c("phenotype_value","n_PCs_covariates","kinship","GWAS_software","n_samples",
                "n_SNPs_tested","n_SNPs_LD_pruned","LD_pruning_call","alpha",
                "Bonferroni_neg_log10_p","Bonferroni_p_value",
                "clump_p1_lead_threshold","clump_p2_member_threshold","clump_r2","clump_kb",
                "max_span_kb","gap_rule_kb","n_loci","n_lead_only_loci","median_span_kb",
                "n_extra_subthreshold_peaks"),
  value = c("BLUP", 3, "EMMAX aIBS 290x290", "EMMAX", 290,
            format(nrow(ps)/4, big.mark=","), format(N_PRUNED, big.mark=","),
            "--indep-pairwise 1000kb 1 0.2", ALPHA,
            sprintf("%.4f", BONF), format(signif(P1,4), scientific=TRUE),
            format(signif(P1,4), scientific=TRUE),
            "1 (no p-value condition on members)", R2, KB, 2*KB, MAX_GAP_KB,
            nrow(S), sum(S$n_members==1L), round(median(S$span_kb),1), if (HAS_EXTRA) nrow(SE) else 0L)),
  file.path(TAB,"Table_analysis_parameters.tsv"), sep="\t")

# ---- handoff table for the candidate-gene search ----
hand <- S[, .(search_id = locus_id, trait, chr, lead_SNP, lead_bp,
              start, end, span_kb = round(span_kb,3), n_member_SNPs = n_members,
              lead_neg_log10_p, lead_MAF = round(lead_MAF,4),
              class = "significant_locus", include_in_gene_search = TRUE)]
if (HAS_EXTRA) hand <- rbind(hand,
  SE[, .(search_id = locus_id, trait, chr, lead_SNP, lead_bp,
         start, end, span_kb = round(span_kb,3), n_member_SNPs = n_members,
         lead_neg_log10_p, lead_MAF = round(lead_MAF,4),
         class = "subthreshold_peak", include_in_gene_search = TRUE)])
setorder(hand, trait, chr, lead_bp)
fwrite(hand, file.path(TAB,"Table_loci_for_gene_search.tsv"), sep="\t")

cat(sprintf("[36] %d loci + %d extra peaks | %d member SNPs | median span %.1f kb\n",
            nrow(S), if (HAS_EXTRA) nrow(SE) else 0L, nrow(MM), median(S$span_kb)))
cat("[36] tables written:\n"); for (f in list.files(TAB, "^Table_")) cat("   ", f, "\n")
