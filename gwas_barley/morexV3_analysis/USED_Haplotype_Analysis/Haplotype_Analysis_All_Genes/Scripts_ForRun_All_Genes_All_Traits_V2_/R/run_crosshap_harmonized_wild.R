suppressPackageStartupMessages({
  library(crosshap)
  library(data.table)
  library(dplyr)
})

run_crosshap_harmonized_wild <- function(
  trait,
  gene_file,
  epsilon_vector,
  MGmin,
  minHap,
  hetmiss_as,
  keep_outliers,
  vcf_root,
  pheno_root,
  plink_bin,
  tmp_dir
) {
  base_name <- make_base_name(trait)

  raw_dir <- make_raw_vcf_dir(vcf_root, trait)
  imp_dir <- make_imputed_vcf_dir(vcf_root, trait)

  raw_path <- file.path(raw_dir, gene_file)
  imp_path <- file.path(imp_dir, gene_file)
  pheno_file <- make_pheno_file(pheno_root, trait)

  stopifnot(file.exists(raw_path))
  stopifnot(file.exists(imp_path))
  stopifnot(file.exists(pheno_file))

  raw_pheno <- fread(pheno_file, header = FALSE)
  colnames(raw_pheno) <- c("FID", "IID", "trait")

  pheno <- raw_pheno %>%
    mutate(Ind = paste(FID, IID, sep = "_")) %>%
    select(Ind, Pheno = trait) %>%
    distinct(Ind, .keep_all = TRUE)

  vcf_raw <- read_vcf(raw_path)
  if (nrow(vcf_raw) == 0) stop("vcf_raw has 0 variants: ", raw_path)

  vcf_imp <- read_vcf(imp_path)
  if (nrow(vcf_imp) == 0) stop("vcf_imp has 0 variants: ", imp_path)

  ord_raw <- order(vcf_raw[["#CHROM"]], vcf_raw[["POS"]], vcf_raw[["REF"]], vcf_raw[["ALT"]])
  ord_imp <- order(vcf_imp[["#CHROM"]], vcf_imp[["POS"]], vcf_imp[["REF"]], vcf_imp[["ALT"]])

  vcf_raw <- vcf_raw[ord_raw, ]
  vcf_imp <- vcf_imp[ord_imp, ]

  raw_id_base <- paste0(vcf_raw[["#CHROM"]], ":", vcf_raw[["POS"]])
  imp_id_base <- paste0(vcf_imp[["#CHROM"]], ":", vcf_imp[["POS"]])

  common_ids <- intersect(raw_id_base, imp_id_base)
  if (length(common_ids) < MGmin) {
    stop("Too few common variants after raw/imputed intersection: ", length(common_ids))
  }

  keep_raw <- raw_id_base %in% common_ids
  keep_imp <- imp_id_base %in% common_ids

  vcf_raw <- vcf_raw[keep_raw, ]
  vcf_imp <- vcf_imp[keep_imp, ]

  raw_id_base <- paste0(vcf_raw[["#CHROM"]], ":", vcf_raw[["POS"]])
  imp_id_base <- paste0(vcf_imp[["#CHROM"]], ":", vcf_imp[["POS"]])

  # Keep the old tool's behavior exactly here because it is already proven to work
  # for this dataset and this test gene.
  vcf_raw$ID <- make.unique(raw_id_base)
  vcf_imp$ID <- make.unique(imp_id_base)

  ensure_dir(tmp_dir)

  tmp_prefix <- file.path(
    tmp_dir,
    paste0("tmp_plink_ld_", gsub("[^A-Za-z0-9_]+", "_", gene_file), "_")
  )

  on.exit({
    tmp_files <- list.files(
      tmp_dir,
      pattern = paste0("^", basename(tmp_prefix)),
      full.names = TRUE
    )
    if (length(tmp_files) > 0) suppressWarnings(file.remove(tmp_files))
  }, add = TRUE)

  tmp_header <- paste0(tmp_prefix, "header.vcf")
  tmp_body <- paste0(tmp_prefix, "body.vcf")
  tmp_ready <- paste0(tmp_prefix, "ready_for_plink.vcf")
  tmp_ld_out <- paste0(tmp_prefix, "ld")

  # Use the old working tool's exact header extraction and PLINK preparation path.
  header_lines <- system(paste0("zgrep \"^#\" ", shQuote(imp_path)), intern = TRUE)
  writeLines(header_lines, tmp_header)
  data.table::fwrite(vcf_imp, tmp_body, sep = "\t", col.names = FALSE, quote = FALSE)
  system(paste("cat", shQuote(tmp_header), shQuote(tmp_body), ">", shQuote(tmp_ready)))

  system(paste(
    shQuote(plink_bin),
    "--vcf", shQuote(tmp_ready),
    "--r2 square --out", shQuote(tmp_ld_out),
    "--allow-extra-chr --double-id --silent"
  ))

  ld_file <- paste0(tmp_ld_out, ".ld")
  if (!file.exists(ld_file)) {
    stop("PLINK LD file not created: ", ld_file)
  }

  LD <- read_LD(ld_file, vcf = vcf_raw)

  HapObject <- run_haplotyping(
    vcf = vcf_raw,
    LD = LD,
    pheno = pheno,
    epsilon = epsilon_vector,
    MGmin = MGmin,
    minHap = minHap,
    hetmiss_as = hetmiss_as,
    keep_outliers = keep_outliers
  )

  list(
    trait = trait,
    base_name = base_name,
    gene_file = gene_file,
    raw_path = raw_path,
    imp_path = imp_path,
    pheno_file = pheno_file,
    common_ids = common_ids,
    vcf_raw_ids = vcf_raw$ID,
    HapObject = HapObject
  )
}