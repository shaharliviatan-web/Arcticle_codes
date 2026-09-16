# Build the LD matrix directly (no haplotyping), then sweep epsilon/MGmin.
Sys.setenv(TMPDIR="/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({library(yaml);library(data.table);library(crosshap);library(dplyr);library(tibble)})
sd_ <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/04_USED_haplotype_analysis_crosshap/01_scripts"
source(file.path(sd_,"R","utils.R")); source(file.path(sd_,"R","run_crosshap.R"))
cfg <- yaml::read_yaml(file.path(dirname(sd_),"00_config","config.yaml"))
gw <- as.data.table(read_gene_windows(file.path(dirname(sd_),"00_config","gene_windows.tsv")))[trait=="protein"]
MG <- c(2L,3L); EPS <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,3.0,5.0)
tmp <- file.path(cfg$output_root,"04_runs","_tmp_protein_eps"); ensure_dir(tmp)
rp <- fread(pheno_path(cfg,"protein"),header=FALSE); colnames(rp) <- c("FID","IID","trait")
pheno <- rp %>% mutate(Ind=as.character(IID)) %>% select(Ind,Pheno=trait) %>% distinct(Ind,.keep_all=TRUE)
out <- list()
for (i in seq_len(nrow(gw))) {
  gf <- gw$gene_file[i]; gn <- make_gene_name(gf)
  raw <- raw_vcf_path(cfg,"protein",gf); imp <- imputed_vcf_path(cfg,"protein",gf)
  vr <- tryCatch(read_vcf_robust(raw), error=function(e) NULL)
  vi <- tryCatch(read_vcf_robust(imp), error=function(e) NULL)
  nraw <- if (is.null(vr)) 0 else nrow(vr)
  if (is.null(vr)||is.null(vi)||nrow(vr)==0||nrow(vi)==0) { cat(sprintf("%-28s %2d SNPs  -> no variants\n",gn,nraw)); next }
  ir <- paste0(vr[["#CHROM"]],":",vr[["POS"]]); ii <- paste0(vi[["#CHROM"]],":",vi[["POS"]])
  common <- intersect(ir,ii)
  if (length(common)<2) { cat(sprintf("%-28s %2d SNPs  -> only %d common variant(s), cannot cluster\n",gn,nraw,length(common))); next }
  vr <- vr[ir %in% common,]; vi <- vi[ii %in% common,]
  vr$ID <- make.unique(paste0(vr[["#CHROM"]],":",vr[["POS"]]))
  pre <- file.path(tmp,paste0("ld_",gn,"_"))
  hdr <- paste0(pre,"h.vcf"); body <- paste0(pre,"b.vcf"); rdy <- paste0(pre,"r.vcf")
  writeLines(system(paste0("zgrep \"^#\" ",shQuote(imp)),intern=TRUE), hdr)
  fwrite(vi,body,sep="\t",col.names=FALSE,quote=FALSE)
  system(paste("cat",shQuote(hdr),shQuote(body),">",shQuote(rdy)))
  system(paste(shQuote(cfg$plink_bin),"--vcf",shQuote(rdy),"--r2 square --keep-allele-order",
               "--allow-extra-chr --double-id --silent --out",shQuote(paste0(pre,"ld"))))
  ldf <- paste0(pre,"ld.ld"); if (!file.exists(ldf)) { cat(sprintf("%-28s PLINK LD failed\n",gn)); next }
  LD <- read_LD(ldf, vcf=vr)
  got <- 0
  for (M in MG) for (E in EPS) {
    r <- tryCatch(suppressWarnings(suppressMessages(run_haplotyping(vcf=vr,LD=LD,pheno=pheno,epsilon=E,MGmin=M,
         minHap=cfg$minHap,hetmiss_as=cfg$hetmiss_as,keep_outliers=cfg$keep_outliers))),error=function(e) NULL)
    if (is.null(r)) next
    H <- r[[paste0("Haplotypes_MGmin",M,"_E",E)]]; if (is.null(H)||is.null(H$Indfile)||!nrow(H$Indfile)) next
    hp <- as.character(H$Indfile$hap); ph <- suppressWarnings(as.numeric(H$Indfile$Pheno))
    k <- hp!="0" & !is.na(ph); ng <- length(unique(hp[k])); if (ng<2) next
    kt <- tryCatch(kruskal.test(ph[k]~hp[k]),error=function(e) NULL); if (is.null(kt)) next
    n <- sum(k); got <- got+1
    out[[length(out)+1]] <- data.table(gene=gn,n_snps=length(common),MGmin=M,eps=E,n_groups=ng,
      n_assigned=sum(hp!="0"), kw_p=kt$p.value, eta2=max(0,min(1,(unname(kt$statistic)-ng+1)/(n-ng))))
  }
  cat(sprintf("%-28s %2d SNPs  -> %d of %d configs gave a grouping\n",gn,length(common),got,length(MG)*length(EPS)))
}
d <- if (length(out)) rbindlist(out) else data.table()
saveRDS(d,"/mnt/data/shahar/.tmp/protein_eps.rds")
if (!nrow(d)) cat("\n*** NO protein gene yields a valid haplotype grouping at ANY of the 32 configs tested.\n") else {
  cat("\n=== configs where a protein gene becomes testable ===\n"); print(d[order(kw_p)])
  cat("\nNOTE: raw KW p only. Protein's BH family would be these genes alone.\n")
}
