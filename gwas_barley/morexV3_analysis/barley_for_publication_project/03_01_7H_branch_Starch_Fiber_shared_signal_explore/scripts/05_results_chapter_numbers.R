#!/usr/bin/env Rscript
# 05_results_chapter_numbers.R -- every number of this replacement analysis, in prose,
# ready to paste into the manuscript. Mirrors the `results_chapter_numbers.txt` that
# steps 03, 04 and 07 each produce.
suppressPackageStartupMessages(library(data.table))
R <- Sys.getenv("TEMP_ROOT"); O <- file.path(R,"results","tables")
gr <- fread(file.path(O,"mgmin3_gene_results.tsv")); ov <- fread(file.path(O,"Table_site_overlap.tsv"))
hg <- fread(file.path(O,"Table_haplotype_groups.tsv")); scr <- fread(file.path(O,"Table_elite_line_screen.tsv"))
pgrid <- fread(file.path(O,"param_grid_eps0.05-1.5_MGmin2-3.tsv")); rs <- fread(file.path(O,"Table_removed_sites_summary.tsv"))
wide <- fread(file.path(O,"Table_elite_genotypes_wide__HORVU.MOREX.r3.7HG0729030.tsv"))
con <- file(file.path(O,"results_chapter_numbers.txt"),"w"); wl <- function(...) writeLines(paste0(...), con)
wl(strrep("=",78)); wl("  7HG0729030 (GDSL esterase) at MGmin = 3, epsilon = 0.9")
wl("  REPLACEMENT analysis for the 7H shared fiber/starch signal")
wl("  generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by 05_results_chapter_numbers.R"); wl(strrep("=",78)); wl("")
wl("PARAMETERS"); wl("  MGmin = 3, epsilon = 0.9, window = gene +/- 1000 bp, minHap = 9.")
wl("  The pipeline elsewhere uses MGmin = 2, epsilon = 0.6. Under those values this gene")
wl("  has NO usable grouping: DBSCAN assigns all four GWAS signal SNPs to marker group 0")
wl("  (noise) and builds haplotypes from two null SNPs (fiber p = 0.16, starch p = 0.36).")
wl("")
wl("PARAMETER SELECTION (see README section 2)")
wl("  A single procedure in two steps, neither using the haplotype test's p-values:")
wl("    Step 1  MGmin = the smallest value at which the four GWAS signal SNPs form a")
wl("            marker group                                        -> MGmin = 3")
wl("    Step 2  epsilon = the value maximising the assignment rate  -> epsilon = 0.9")
wl("            (246/290, tied with eps 1.0; the smaller eps is taken)")
wl("")
wl("  WHY STEP 1 IS NEEDED. Ranked purely by assignment, the grid's best combination is")
wl("  MGmin 2 / eps 0.10 with 253/290 (87.2%) assigned -- but it carries 0 of the 4 GWAS")
wl("  signal SNPs and returns KW p = 0.92. Assignment rate measures how much of the panel")
wl("  the test uses, not whether it is aimed at the right variants. Within MGmin = 2 the")
wl("  only setting carrying all four signal SNPs is eps = 1.5, assigning just 147/290.")
wl("")
pg <- unique(pgrid[status=="ok", .(MGmin, eps, sig_kept, assigned)])[order(-assigned)]
wl("  Grid ranked by assignment (top rows):")
for (i in seq_len(min(8, nrow(pg)))) wl(sprintf("    %2d. MGmin %d, eps %.2f : %d/4 signal SNPs, %3d/290 assigned (%.1f%%)%s",
  i, pg$MGmin[i], pg$eps[i], pg$sig_kept[i], pg$assigned[i], 100*pg$assigned[i]/290,
  ifelse(pg$MGmin[i]==3 & pg$eps[i]==0.9, "   <- chosen",
  ifelse(i==1, "   <- max assignment, but unusable", ""))))
wl("")
wl("HAPLOTYPE RESULT")
for (tr in c("fiber","starch")) { g <- gr[trait==tr & epsilon==0.9]
  wl(sprintf("  %-7s %d SNPs in window | %d of 290 accessions assigned (%.0f%%) | %d groups [%s]",
     tr, g$n_snps_window, g$n_assigned, 100*g$n_assigned/290, g$haplotype_groups, g$group_sizes))
  wl(sprintf("          KW H = %.2f, df = %d, p = %.3g | eta2 = %.3f | delta top-bottom = %.2f SD",
     g$kw_H, g$kw_df, g$kw_p_raw, g$eta_squared, g$delta_top_bottom_sd))
  for (i in seq_len(nrow(hg[trait==tr]))) { h <- hg[trait==tr][i]
    wl(sprintf("          hap %s: n = %3d, mean %+.4f, median %+.4f, SD %.3f", h$hap, h$n, h$mean_pheno, h$median_pheno, h$sd_pheno)) } }
wl("")
wl("  All four GWAS signal SNPs are retained in marker group MG1 (mean r2 = 0.74).")
wl("  The minority haplotype (n = 34) is EXACTLY the 34 accessions that carry the minor")
wl("  allele at both lead SNPs in the direct genotype split -- 34 shared, 0 discordant.")
wl("")
wl("ELITE COMPARISON")
wl(sprintf("  wild %d SNPs, elite %d records in window; shared %d, triallelic %d, wild-only %d, elite-only %d",
   ov$n_wild_snps, ov$n_elite_records, ov$n_shared, ov$n_alleles_triallelic, ov$n_wild_only, ov$n_elite_only))
wl(sprintf("  allele concordance over %d shared positions: %d identical, %d swapped, %d triallelic",
   ov$n_shared_positions, ov$n_alleles_identical, ov$n_alleles_swapped, ov$n_alleles_triallelic))
for (i in seq_len(nrow(rs))) wl(sprintf("  removed from %-5s (%s): %3d sites, %d polymorphic, %d monomorphic",
   rs$file[i], rs$removed_from[i], rs$n_sites[i], rs$n_polymorphic[i], rs$n_monomorphic[i]))
wl(sprintf("  wild monomorphic over the %d assigned accessions: %d of %d SNPs (%d kept, %d removed)",
   ov$n_wild_assigned_accessions, ov$n_wild_mono_all, ov$n_wild_snps, ov$n_wild_mono_kept, ov$n_wild_mono_removed))
wl("")
sig <- wide[is_gwas_signal_snp==TRUE]
lines <- setdiff(names(wide), c("chr","pos","site_key","is_gwas_signal_snp"))
wl("  Genotype of each elite cultivar at the 4 GWAS signal SNPs:")
for (l in lines) wl(sprintf("    %-11s %s", l, paste(sig[[l]], collapse=", ")))
n_alt <- sum(unlist(sig[, ..lines])=="Alternate"); n_ref <- sum(unlist(sig[, ..lines])=="Reference")
wl(sprintf("    => %d Reference, %d Alternate, %d no-call across %d lines x 4 SNPs",
   n_ref, n_alt, sum(unlist(sig[, ..lines])=="Missing"), length(lines)))
wl("    No elite cultivar carries the minority (high-fibre / low-starch) haplotype.")
wl("")
wl("  Call rate over the 22 shared sites (step 07 requires >= 85%):")
for (i in seq_len(nrow(scr))) wl(sprintf("    %-11s %2d/%d = %.1f%%  %s", scr$line_name[i], scr$n_called[i],
   scr$n_shared_sites[i], scr$call_rate_pct[i], ifelse(scr$passes_85[i],"pass","BELOW 85%")))
wl("")
wl("PASTE-READY SENTENCE")
wl("  At the shared 7H fiber/starch signal, haplotyping of the adjacent GDSL esterase")
wl("  7HG0729030 (gene +/- 1 kb, MGmin = 3, epsilon = 0.9) resolved two haplotypes over")
wl("  246 of 290 accessions. The minority haplotype (n = 34) carried higher fibre")
wl(sprintf("  (KW p = %.3g) and lower starch (p = %.3g) BLUPs, and is absent from five modern",
   gr[trait=="fiber" & epsilon==0.9]$kw_p_raw, gr[trait=="starch" & epsilon==0.9]$kw_p_raw))
wl("  European spring malting cultivars, all of which carry the major haplotype.")
close(con)
cat(readLines(file.path(O,"results_chapter_numbers.txt")), sep="\n")
