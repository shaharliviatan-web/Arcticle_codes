#!/usr/bin/env Rscript
# 04_paper_tables.R -- the full step-07 table set, for 7HG0729030 at MGmin=3, eps=0.9.
# Column names and definitions are copied from 07_.../scripts/05_paper_tables.R and
# 06_removed_sites_table.R so a future session can join these to the published tables.
#
# "polymorphic" is judged over the SAME reference sets step 07 uses:
#   wild  -- only the accessions crosshap ASSIGNED to a haplotype group (246 here)
#   elite -- only the 5 CONFIGURED lines, not the 136-line pool
suppressPackageStartupMessages({library(vcfR); library(data.table)})
R <- Sys.getenv("TEMP_ROOT"); P <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
E <- file.path(P,"07_USED_elite_lines_compariosn_to_wild_lines")
GENE<-"HORVU.MOREX.r3.7HG0729030"; SHORT<-"GDSL"; CHR<-"7H"
GS<-573604051; GE<-573605728; WS<-573603051; WE<-573606728; WIN<-1000L
MG<-3L; EPS<-0.9; SIG<-c(573606282,573606306,573606460,573606491)
LOCI <- list(fiber=list(locus="fiber_L17", lead="7H:573606306", lead_pos=573606306, p=7.2451),
             starch=list(locus="starch_L05", lead="7H:573606460", lead_pos=573606460, p=6.4367))
OUT <- file.path(R,"results","tables"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
elite_cfg <- fread(cmd=paste("grep -v '^#'", shQuote(file.path(E,"config/elite_lines.tsv"))))
gt_state <- function(g){o<-rep("Missing",length(g)); o[g %in% c("0/0","0|0")]<-"Reference"
  o[g %in% c("1/1","1|1")]<-"Alternate"; o[g %in% c("0/1","1/0","0|1","1|0")]<-"Heterozygous"; o}
read_v <- function(p){ v<-suppressMessages(read.vcfR(p,verbose=FALSE))
  fx<-as.data.frame(getFIX(v),stringsAsFactors=FALSE); g<-extract.gt(v,"GT",as.numeric=FALSE)
  fx$CHROM<-sub("^chr","",fx$CHROM); fx$POS<-as.integer(fx$POS)
  list(fx=fx, gt=g, snv=nchar(fx$REF)==1 & nchar(fx$ALT)==1,
       key=paste(fx$CHROM,fx$POS,fx$REF,fx$ALT,sep=":"), pk=paste(fx$CHROM,fx$POS,sep=":")) }
W <- read_v(file.path(R,"work/raw/fiber",paste0(GENE,".vcf.gz")))
L <- read_v(file.path(R,"work/elite/7HG0729030.elite.vcf.gz"))
ho <- readRDS(file.path(R,"Cache","fiber",GENE,"MGmin_3","eps_0.9","HapObject.rds"))
ind <- as.data.table(ho$HapObject[[paste0("Haplotypes_MGmin",MG,"_E",EPS)]]$Indfile)
ind[,`:=`(Ind=as.character(Ind),hap=as.character(hap),Pheno=as.numeric(Pheno))]
assigned <- ind[hap!="0" & !is.na(Pheno)]
Wsnv <- which(W$snv); Lsnv <- which(L$snv)
shared_keys <- intersect(W$key[Wsnv], L$key[Lsnv])
shared_pos  <- intersect(W$pk[Wsnv], L$pk[Lsnv])
# --- allele concordance at every shared POSITION -----------------------------
ac <- rbindlist(lapply(shared_pos, function(pp){
  wi<-Wsnv[W$pk[Wsnv]==pp][1]; li<-Lsnv[L$pk[Lsnv]==pp][1]
  st <- if (W$fx$REF[wi]!=L$fx$REF[li]) "ref_differs"
        else if (W$fx$ALT[wi]==L$fx$ALT[li]) "identical"
        else if (W$fx$ALT[wi]==L$fx$REF[li] && W$fx$REF[wi]==L$fx$ALT[li]) "swapped" else "triallelic"
  data.table(gene_id=GENE, short_name=SHORT, chr=CHR, pos=as.integer(sub(".*:","",pp)),
    wild_ref=W$fx$REF[wi], wild_alt=W$fx$ALT[wi], elite_ref=L$fx$REF[li], elite_alt=L$fx$ALT[li], status=st)}))
setorder(ac, pos); fwrite(ac, file.path(OUT,"Table_allele_concordance.tsv"), sep="\t")
stopifnot(!any(ac$status=="swapped"))          # step 07's hard check
tri <- ac[status=="triallelic"]
fwrite(tri, file.path(OUT,"Table_triallelic_sites.tsv"), sep="\t")
# --- polymorphism over the two reference sets --------------------------------
poly_w <- function(i){ s<-gt_state(W$gt[i, colnames(W$gt) %in% assigned$Ind]); s<-s[s!="Missing"]
  if(!length(s)) "undetermined" else if(length(unique(s))>1) "polymorphic" else "monomorphic" }
elite_ids <- elite_cfg$sample_id
poly_l <- function(i){ s<-gt_state(L$gt[i, colnames(L$gt) %in% elite_ids]); s<-s[s!="Missing"]
  if(!length(s)) "undetermined" else if(length(unique(s))>1) "polymorphic" else "monomorphic" }
wild_only  <- setdiff(W$key[Wsnv], shared_keys); elite_only <- setdiff(L$key[Lsnv], shared_keys)
rm_rows <- rbindlist(list(
  if(length(wild_only)) data.table(file="wild", removed_from="wild_only", site_key=wild_only,
     pos=as.integer(sub("^[^:]+:([0-9]+):.*","\\1",wild_only)),
     poly=vapply(match(wild_only,W$key), poly_w, character(1)), is_indel=FALSE),
  if(length(elite_only)) data.table(file="elite", removed_from="elite_only", site_key=elite_only,
     pos=as.integer(sub("^[^:]+:([0-9]+):.*","\\1",elite_only)),
     poly=vapply(match(elite_only,L$key), poly_l, character(1)), is_indel=FALSE),
  if(sum(!L$snv)) data.table(file="elite", removed_from="elite_only", site_key=L$key[!L$snv],
     pos=L$fx$POS[!L$snv], poly="undetermined", is_indel=TRUE)), fill=TRUE)
rm_rows[, `:=`(gene_id=GENE, short_name=SHORT)]
setorder(rm_rows, file, pos)
fwrite(rm_rows[, .(gene_id,short_name,file,removed_from,chr=CHR,pos,site_key,polymorphic=poly,is_indel)],
       file.path(OUT,"Table_removed_sites.tsv"), sep="\t")
summ <- rm_rows[, .(n_sites=.N, n_polymorphic=sum(poly=="polymorphic"),
                    n_monomorphic=sum(poly=="monomorphic"), n_undetermined=sum(poly=="undetermined"),
                    n_indel=sum(is_indel)), by=.(short_name, file, removed_from)]
summ[, trait := "fiber+starch (same window)"]
fwrite(summ[, .(short_name, trait, file, removed_from, n_sites, n_polymorphic, n_monomorphic, n_undetermined, n_indel)],
       file.path(OUT,"Table_removed_sites_summary.tsv"), sep="\t")
# --- site overlap ------------------------------------------------------------
mono_w_all <- sum(vapply(Wsnv, poly_w, character(1))=="monomorphic")
mono_w_kept<- sum(vapply(match(shared_keys,W$key), poly_w, character(1))=="monomorphic")
mono_l_all <- sum(vapply(Lsnv, poly_l, character(1))=="monomorphic")
mono_l_kept<- sum(vapply(match(shared_keys,L$key), poly_l, character(1))=="monomorphic")
ov <- data.table(gene_id=GENE, short_name=SHORT, trait="fiber+starch", chr=CHR, win_start=WS, win_end=WE,
  n_wild_snps=length(Wsnv), n_elite_records=nrow(L$fx), n_shared=length(shared_keys),
  n_wild_only=length(wild_only), n_elite_only=length(elite_only), n_elite_only_indel=sum(!L$snv),
  n_shared_positions=length(shared_pos), n_alleles_identical=sum(ac$status=="identical"),
  n_alleles_swapped=sum(ac$status=="swapped"), n_alleles_triallelic=nrow(tri),
  n_alleles_ref_differs=sum(ac$status=="ref_differs"),
  n_wild_assigned_accessions=nrow(assigned), n_wild_mono_all=mono_w_all, n_wild_mono_kept=mono_w_kept,
  n_wild_mono_removed=mono_w_all-mono_w_kept, n_elite_lines_shown=nrow(elite_cfg),
  n_elite_mono_all=mono_l_all, n_elite_mono_kept=mono_l_kept, n_elite_mono_removed=mono_l_all-mono_l_kept)
fwrite(ov, file.path(OUT,"Table_site_overlap.tsv"), sep="\t")
# --- gene windows + haplotype groups (per trait) -----------------------------
gr <- fread(file.path(OUT,"mgmin3_gene_results.tsv"))
gwt <- rbindlist(lapply(c("fiber","starch"), function(tr){ g <- gr[trait==tr & epsilon==EPS]; lo <- LOCI[[tr]]
  data.table(gene_id=GENE, short_name=SHORT, trait=tr, chr=CHR, gene_start=GS, gene_end=GE, strand="-",
    win_start=WS, win_end=WE, locus_id=lo$locus, lead_SNP=lo$lead, lead_pos=lo$lead_pos,
    lead_neg_log10p=lo$p, dist_to_lead_bp=max(0,max(lo$lead_pos-GE, GS-lo$lead_pos)),
    n_snps_window=g$n_snps_window, n_ind=g$n_assigned, n_unassigned=g$n_unassigned,
    n_groups=g$haplotype_groups, group_sizes=g$group_sizes, kw_p_raw=g$kw_p_raw,
    fdr_q=NA_real_, eta_squared=g$eta_squared, delta_top_bottom_sd=g$delta_top_bottom_sd,
    description="Secreted GDSL esterase/lipase (SGNH hydrolase)", window_bp=WIN,
    crosshap_run="03_01_7H_branch_MGmin3_eps0.9", crosshap_epsilon=EPS, crosshap_MGmin=MG)}))
fwrite(gwt, file.path(OUT,"Table_gene_windows.tsv"), sep="\t")
hg <- fread(file.path(OUT,"mgmin3_haplotype_groups.tsv"))[epsilon==EPS]
fwrite(hg[, .(gene_id=GENE, short_name=SHORT, trait, hap, n, mean_pheno=mean, median_pheno=median, sd_pheno=sd)],
       file.path(OUT,"Table_haplotype_groups.tsv"), sep="\t")
# --- elite tables ------------------------------------------------------------
fwrite(elite_cfg, file.path(OUT,"Table_elite_lines.tsv"), sep="\t")
em <- L$gt[match(shared_keys, L$key), , drop=FALSE]
scr <- data.table(line_name=elite_cfg$line_name, sample_id=elite_cfg$sample_id,
  n_shared_sites=length(shared_keys),
  n_called=vapply(elite_cfg$sample_id, function(s) sum(gt_state(em[,s])!="Missing"), integer(1)))
scr[, call_rate_pct := round(100*n_called/n_shared_sites,1)][, passes_85 := call_rate_pct>=85]
fwrite(scr, file.path(OUT,"Table_elite_line_screen.tsv"), sep="\t")
wide <- data.table(chr=CHR, pos=as.integer(sub("^[^:]+:([0-9]+):.*","\\1",shared_keys)),
                   site_key=shared_keys, is_gwas_signal_snp=as.integer(sub("^[^:]+:([0-9]+):.*","\\1",shared_keys)) %in% SIG)
for (i in seq_len(nrow(elite_cfg))) wide[[elite_cfg$line_name[i]]] <- gt_state(em[, elite_cfg$sample_id[i]])
setorder(wide, pos); fwrite(wide, file.path(OUT, paste0("Table_elite_genotypes_wide__",GENE,".tsv")), sep="\t")
fwrite(data.table(gene_id=GENE, chr=CHR, gene_start=GS, gene_end=GE, win_start=WS, win_end=WE,
  requested="chr7H:573601051-573608728", pad_bp=2000L, n_returned=nrow(L$fx),
  n_in_window=length(Lsnv), n_samples=ncol(L$gt), source="DivBrowse barley_pangenome_v2 /vcf_export",
  fetched_utc=format(Sys.time(), tz="UTC", usetz=TRUE)),
  file.path(OUT,"elite_vcf_provenance.tsv"), sep="\t")
cat(sprintf("04 OK | shared %d | wild-only %d | elite-only %d (indel %d) | triallelic %d\n",
  length(shared_keys), length(wild_only), length(elite_only), sum(!L$snv), nrow(tri)))
print(summ[, .(file, removed_from, n_sites, n_polymorphic, n_monomorphic, n_undetermined)])
print(scr[, .(line_name, n_called, n_shared_sites, call_rate_pct, passes_85)])
