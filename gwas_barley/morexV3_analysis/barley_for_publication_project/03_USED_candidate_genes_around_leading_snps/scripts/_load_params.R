#!/usr/bin/env Rscript
# _load_params.R -- make config/params.sh readable from R.
# Sourced at the top of every R script here. Parses the `export KEY=VALUE` lines
# so bash and R can never drift apart on a parameter value.

load_params <- function(path = NULL) {
  if (is.null(path)) {
    path <- file.path(
      "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project",
      "03_USED_candidate_genes_around_leading_snps", "config", "params.sh")
  }
  stopifnot(file.exists(path))
  ln <- readLines(path, warn = FALSE)
  ln <- grep("^export [A-Za-z_][A-Za-z0-9_]*=", ln, value = TRUE)
  key <- sub("^export ([A-Za-z_][A-Za-z0-9_]*)=.*$", "\\1", ln)
  val <- sub("^export [A-Za-z_][A-Za-z0-9_]*=", "", ln)
  val <- gsub('^"|"$', "", val)                       # strip surrounding quotes
  p <- as.list(setNames(val, key))
  # Expand $VAR / ${VAR} references in declaration order.
  for (i in seq_along(p)) {
    for (k in names(p)[seq_len(i - 1)]) {
      p[[i]] <- gsub(paste0("\\$\\{", k, "\\}"), p[[k]], p[[i]])
      p[[i]] <- gsub(paste0("\\$", k, "\\b"),    p[[k]], p[[i]])
    }
  }
  p$FLANK_BP <- as.integer(p$FLANK_BP)
  p$FLANK_SENSITIVITY_SET <- as.integer(strsplit(trimws(p$FLANK_SENSITIVITY_SET), "\\s+")[[1]])
  Sys.setenv(TMPDIR = p$TMPDIR)
  p
}
