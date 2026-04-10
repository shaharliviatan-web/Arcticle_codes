suppressPackageStartupMessages({
  library(crosshap)
  library(data.table)
})

run_crosshap_intersection_wild_elite <- function(target, cfg, tmp_dir) {
  target <- normalize_target(target)

  merged_path <- make_merged_vcf_path(cfg$merged_intersection_vcf_dir, target$merged_vcf_file)
  ld_path <- make_imputed_vcf_path(cfg$ld_vcf_root, target$trait, target$gene_file)
  pheno_file <- make_pheno_file(cfg$pheno_root, target$trait)

  stopifnot(file.exists(merged_path))
  stopifnot(file.exists(ld_path))
  stopifnot(file.exists(pheno_file))

  pheno <- load_pheno_wild_elite(pheno_file)

  vcf_merged <- read_vcf(merged_path)
  if (nrow(vcf_merged) == 0) stop("Merged intersection VCF has 0 variants: ", merged_path)

  vcf_ld <- read_vcf(ld_path)
  if (nrow(vcf_ld) == 0) stop("LD VCF has 0 variants: ", ld_path)

  ord_merged <- order(vcf_merged[["#CHROM"]], vcf_merged[["POS"]], vcf_merged[["REF"]], vcf_merged[["ALT"]])
  ord_ld <- order(vcf_ld[["#CHROM"]], vcf_ld[["POS"]], vcf_ld[["REF"]], vcf_ld[["ALT"]])

  vcf_merged <- vcf_merged[ord_merged, ]
  vcf_ld <- vcf_ld[ord_ld, ]

  merged_pos <- paste0(vcf_merged[["#CHROM"]], ":", vcf_merged[["POS"]])
  ld_pos <- paste0(vcf_ld[["#CHROM"]], ":", vcf_ld[["POS"]])
  common_pos <- intersect(merged_pos, ld_pos)

  if (length(common_pos) < target$MGmin) {
    stop("Too few common variants between merged intersection VCF and LD VCF: ", length(common_pos))
  }

  keep_merged <- merged_pos %in% common_pos
  keep_ld <- ld_pos %in% common_pos

  vcf_merged <- vcf_merged[keep_merged, ]
  vcf_ld <- vcf_ld[keep_ld, ]

  merged_pos <- paste0(vcf_merged[["#CHROM"]], ":", vcf_merged[["POS"]])
  ld_pos <- paste0(vcf_ld[["#CHROM"]], ":", vcf_ld[["POS"]])

  vcf_merged$ID <- make.unique(merged_pos)
  vcf_ld$ID <- make.unique(ld_pos)

  ensure_dir(tmp_dir)

  tmp_prefix <- file.path(tmp_dir, paste0("tmp_plink_ld_", gsub("[^A-Za-z0-9_]+", "_", target$gene_file), "_"))
  on.exit({
    tmp_files <- list.files(tmp_dir, pattern = paste0("^", basename(tmp_prefix)), full.names = TRUE)
    if (length(tmp_files) > 0) suppressWarnings(file.remove(tmp_files))
  }, add = TRUE)

  tmp_header <- paste0(tmp_prefix, "header.vcf")
  tmp_body <- paste0(tmp_prefix, "body.vcf")
  tmp_ready <- paste0(tmp_prefix, "ready_for_plink.vcf")
  tmp_ld_out <- paste0(tmp_prefix, "ld")

  header_lines <- system(paste0("zgrep \"^#\" ", shQuote(ld_path)), intern = TRUE)
  writeLines(header_lines, tmp_header)
  fwrite(vcf_ld, tmp_body, sep = "\t", col.names = FALSE, quote = FALSE)
  system(paste("cat", shQuote(tmp_header), shQuote(tmp_body), ">", shQuote(tmp_ready)))

  system(paste(
    shQuote(cfg$plink_bin),
    "--vcf", shQuote(tmp_ready),
    "--r2 square --out", shQuote(tmp_ld_out),
    "--allow-extra-chr --double-id --silent"
  ))

  ld_file <- paste0(tmp_ld_out, ".ld")
  if (!file.exists(ld_file)) stop("PLINK LD file not created: ", ld_file)

  LD <- read_LD(ld_file, vcf = vcf_merged)

  HapObject <- run_haplotyping(
    vcf = vcf_merged,
    LD = LD,
    pheno = pheno,
    epsilon = target$epsilon,
    MGmin = target$MGmin,
    minHap = cfg$minHap,
    hetmiss_as = cfg$hetmiss_as,
    keep_outliers = cfg$keep_outliers
  )

  list(
    target = target,
    merged_path = merged_path,
    ld_path = ld_path,
    pheno_file = pheno_file,
    common_pos = common_pos,
    vcf_merged_ids = vcf_merged$ID,
    HapObject = HapObject,
    pheno = pheno
  )
}
