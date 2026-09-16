# =============================================================================
# _load_params.R -- make config/params.sh readable from R
# =============================================================================
# Bash sources config/params.sh directly. R parses the SAME file here, so a
# parameter can never drift between the two. Copied in spirit from step 03's
# scripts/_load_params.R.
#
# Usage:  P <- load_params()            # named list of scalars
#         P$WINDOW_BP ; P$TARGET_GENES  # arrays come back as character vectors
# =============================================================================

load_params <- function(params_file = NULL) {
  if (is.null(params_file)) {
    args <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args, value = TRUE)
    here <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
    params_file <- file.path(here, "..", "config", "params.sh")
  }
  if (!file.exists(params_file)) stop("params.sh not found: ", params_file)
  params_file <- normalizePath(params_file)

  lines <- readLines(params_file, warn = FALSE)
  P <- list()

  # --- plain scalars:  NAME="value"  /  NAME=value  ---------------------------
  scal <- grep('^[A-Za-z_][A-Za-z0-9_]*=', lines)
  for (i in scal) {
    l <- lines[i]
    nm <- sub('=.*$', '', l)
    v  <- sub('^[^=]*=', '', l)
    v  <- sub('\\s+#.*$', '', v)            # strip trailing comment
    v  <- gsub('^"|"$', '', trimws(v))
    P[[nm]] <- v
  }

  # --- bash arrays:  NAME=( ... )  spanning to the closing ')' -----------------
  arr_start <- grep('^[A-Za-z_][A-Za-z0-9_]*=\\(', lines)
  for (i in arr_start) {
    nm <- sub('=\\(.*$', '', lines[i])
    j <- i
    while (j <= length(lines) && !grepl('\\)', lines[j])) j <- j + 1
    body <- paste(lines[i:j], collapse = " ")
    body <- sub('^[^(]*\\(', '', body)
    body <- sub('\\).*$', '', body)
    vals <- regmatches(body, gregexpr('"[^"]*"', body))[[1]]
    P[[nm]] <- gsub('^"|"$', '', vals)
  }

  # --- expand ${VAR} references, repeatedly until stable ----------------------
  expand_once <- function(x, env) {
    for (k in names(env)) {
      if (length(env[[k]]) != 1) next
      x <- gsub(paste0("\\$\\{", k, "\\}"), env[[k]], x, fixed = FALSE)
    }
    x
  }
  for (round in 1:6) {
    P <- lapply(P, function(v) expand_once(v, P))
  }

  # --- typed conveniences -----------------------------------------------------
  num <- c("WINDOW_BP", "DIVBROWSE_PAD_BP", "DIVBROWSE_PAD_MAX_BP",
           "DIVBROWSE_RETRY_FACTOR", "DIVBROWSE_TIMEOUT", "CROSSHAP_MGMIN",
           "BIOSAMPLES_PARALLEL", "FIG_WIDTH_IN", "FIG_HEIGHT_IN", "FIG_DPI",
           "CALL_RATE_MIN", "TILE_HEIGHT")
  for (n in num) if (!is.null(P[[n]])) P[[n]] <- as.numeric(P[[n]])

  P
}

# Split "gene_id:short_name" entries of TARGET_GENES into a data.frame.
target_genes_df <- function(P) {
  parts <- strsplit(P$TARGET_GENES, ":", fixed = TRUE)
  data.frame(
    gene_id    = vapply(parts, `[`, character(1), 1),
    short_name = vapply(parts, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
}

# Read config/elite_lines.tsv, dropping the '#' commentary header block.
read_elite_lines <- function(P) {
  ln <- readLines(P$ELITE_LINES_TSV, warn = FALSE)
  ln <- ln[!grepl('^#', ln)]
  ln <- ln[nzchar(trimws(ln))]
  utils::read.delim(text = paste(ln, collapse = "\n"),
                    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}
