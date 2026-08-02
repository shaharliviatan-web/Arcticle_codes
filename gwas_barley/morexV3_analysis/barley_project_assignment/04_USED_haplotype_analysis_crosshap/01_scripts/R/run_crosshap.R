suppressPackageStartupMessages({
  library(crosshap)
  library(data.table)
  library(dplyr)
  library(tibble)
})

# ---------------------------------------------------------------------------
# run_crosshap.R -- run crosshap for one gene x MGmin (all epsilon at once).
#
# Design (kept from the proven tool):
#   * RAW per-gene VCF  -> haplotyping + genotype display
#   * IMPUTED per-gene VCF -> PLINK LD (r2 is flip-invariant)
#   * variants intersected by CHROM:POS; IDs = make.unique("chr:pos")
#
# Cleaned for step 04:
#   * inputs are bcftools-cut, single-IID, gene+/-window VCFs (no harmonization)
#   * phenotype Ind = IID (single) to match VCF sample names; 290-join asserted
#   * PLINK = absolute path + --keep-allele-order
# ---------------------------------------------------------------------------

VCF_META_COLS <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

# Robust VCF reader: reproduces crosshap::read_vcf's transforms (collapse genotype
# to first field, "/"->"|", numeric POS, "-"->"." in colnames) but strips the "##"
# meta lines first so data.table::fread never mis-detects the separator from a
# comma-containing header line (e.g. ##FORMAT=<ID=AD,Number=R,...>). The "#CHROM"
# line becomes the header. Works for plain or bgzipped VCFs.
read_vcf_robust <- function(path) {
  cmd <- paste("zcat -f", shQuote(path), "| grep -v '^##'")
  vcf <- data.table::fread(cmd = cmd, header = TRUE, sep = "\t", nThread = 4) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                ~ base::gsub(":.*", "", base::gsub("/", "|", .)))) %>%
    dplyr::mutate(POS = as.numeric(.data$POS))
  colnames(vcf) <- gsub("-", ".", colnames(vcf))
  vcf
}

run_crosshap <- function(cfg, trait, gene_file, MGmin, tmp_dir) {
  raw_path   <- raw_vcf_path(cfg, trait, gene_file)
  imp_path   <- imputed_vcf_path(cfg, trait, gene_file)
  pheno_file <- pheno_path(cfg, trait)

  stopifnot(file.exists(raw_path), file.exists(imp_path), file.exists(pheno_file))

  ## ---- Phenotype: Ind = IID (single), matching VCF sample names ----
  raw_pheno <- fread(pheno_file, header = FALSE)
  colnames(raw_pheno) <- c("FID", "IID", "trait")
  pheno <- raw_pheno %>%
    mutate(Ind = as.character(IID)) %>%
    select(Ind, Pheno = trait) %>%
    distinct(Ind, .keep_all = TRUE)

  ## ---- Read VCFs ----
  vcf_raw <- read_vcf_robust(raw_path)
  if (nrow(vcf_raw) == 0) stop("vcf_raw has 0 variants: ", raw_path)
  vcf_imp <- read_vcf_robust(imp_path)
  if (nrow(vcf_imp) == 0) stop("vcf_imp has 0 variants: ", imp_path)

  ## ---- Assert sample/phenotype join (catches ID-convention bugs early) ----
  vcf_samples <- setdiff(colnames(vcf_raw), VCF_META_COLS)
  n_shared <- length(intersect(vcf_samples, pheno$Ind))
  if (n_shared != length(vcf_samples)) {
    stop(sprintf("Sample/pheno join mismatch: only %d of %d VCF samples have a phenotype row (Ind=IID convention?).",
                 n_shared, length(vcf_samples)))
  }

  ## ---- Intersect raw/imputed by CHROM:POS (proven logic, unchanged) ----
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

  vcf_raw <- vcf_raw[raw_id_base %in% common_ids, ]
  vcf_imp <- vcf_imp[imp_id_base %in% common_ids, ]
  raw_id_base <- paste0(vcf_raw[["#CHROM"]], ":", vcf_raw[["POS"]])
  imp_id_base <- paste0(vcf_imp[["#CHROM"]], ":", vcf_imp[["POS"]])

  vcf_raw$ID <- make.unique(raw_id_base)
  vcf_imp$ID <- make.unique(imp_id_base)

  ## ---- PLINK LD from the imputed VCF ----
  ensure_dir(tmp_dir)
  tmp_prefix <- file.path(tmp_dir, paste0("tmp_plink_ld_", gsub("[^A-Za-z0-9_]+", "_", gene_file),
                                          "_MG", MGmin, "_"))
  on.exit({
    tf <- list.files(tmp_dir, pattern = paste0("^", basename(tmp_prefix)), full.names = TRUE)
    if (length(tf) > 0) suppressWarnings(file.remove(tf))
  }, add = TRUE)

  tmp_header <- paste0(tmp_prefix, "header.vcf")
  tmp_body   <- paste0(tmp_prefix, "body.vcf")
  tmp_ready  <- paste0(tmp_prefix, "ready_for_plink.vcf")
  tmp_ld_out <- paste0(tmp_prefix, "ld")

  header_lines <- system(paste0("zgrep \"^#\" ", shQuote(imp_path)), intern = TRUE)
  writeLines(header_lines, tmp_header)
  data.table::fwrite(vcf_imp, tmp_body, sep = "\t", col.names = FALSE, quote = FALSE)
  system(paste("cat", shQuote(tmp_header), shQuote(tmp_body), ">", shQuote(tmp_ready)))

  rc <- system(paste(
    shQuote(cfg$plink_bin),
    "--vcf", shQuote(tmp_ready),
    "--r2 square",
    "--keep-allele-order",
    "--allow-extra-chr --double-id --silent",
    "--out", shQuote(tmp_ld_out)
  ))

  ld_file <- paste0(tmp_ld_out, ".ld")
  if (!file.exists(ld_file)) stop("PLINK LD file not created: ", ld_file, " (plink rc=", rc, ")")

  LD <- read_LD(ld_file, vcf = vcf_raw)

  ## ---- Haplotyping (raw VCF), all epsilon at once ----
  HapObject <- run_haplotyping(
    vcf = vcf_raw, LD = LD, pheno = pheno,
    epsilon = cfg$epsilon_vector, MGmin = MGmin, minHap = cfg$minHap,
    hetmiss_as = cfg$hetmiss_as, keep_outliers = cfg$keep_outliers
  )

  list(trait = trait, gene_file = gene_file, raw_path = raw_path, imp_path = imp_path,
       pheno_file = pheno_file, common_ids = common_ids, n_common = length(common_ids),
       vcf_raw_ids = vcf_raw$ID, HapObject = HapObject)
}
