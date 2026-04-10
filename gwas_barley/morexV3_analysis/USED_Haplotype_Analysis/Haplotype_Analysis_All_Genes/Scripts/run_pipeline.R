suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(dplyr)
})

argv <- commandArgs(trailingOnly = FALSE)
file_arg <- argv[grepl("^--file=", argv)]
if (length(file_arg) == 0) stop("Run with: Rscript [run_pipeline.R](http://_vscodecontentref_/7) (not source()).")

script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
script_dir <- dirname(script_path)

source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "run_crosshap_harmonized_wild.R"))
source(file.path(script_dir, "R", "plot_combined_pdf.R"))
source(file.path(script_dir, "R", "plot_heatmaps.R"))

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path(script_dir, "config.yaml")
targets_path <- if (length(args) >= 2) args[2] else file.path(script_dir, "targets.yaml")

cfg <- yaml::read_yaml(config_path)

cfg$traits <- cfg$traits %||% c("betaglucan", "fiber", "starch")
cfg$mgmin_values <- cfg$mgmin_values %||% c(2, 3)
cfg$epsilon_vector <- cfg$epsilon_vector %||% c(0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 0.85)
cfg$minHap <- cfg$minHap %||% 9
cfg$hetmiss_as <- cfg$hetmiss_as %||% "allele"
cfg$keep_outliers <- cfg$keep_outliers %||% FALSE
cfg$plink_bin <- cfg$plink_bin %||% "plink"
cfg$overwrite_cache <- cfg$overwrite_cache %||% FALSE
cfg$continue_on_error <- cfg$continue_on_error %||% TRUE
cfg$fail_fast <- cfg$fail_fast %||% FALSE
cfg$stats <- cfg$stats %||% list()
cfg$stats$alpha <- cfg$stats$alpha %||% 0.05
cfg$stats$within_gene_method <- cfg$stats$within_gene_method %||% "holm"
cfg$stats$across_gene_fdr_method <- cfg$stats$across_gene_fdr_method %||% "BH"
cfg$stats$across_gene_bonferroni <- cfg$stats$across_gene_bonferroni %||% TRUE

run_root <- file.path(cfg$output_root, "Runs", cfg$run_id)
cache_root <- file.path(run_root, "Cache")
combined_root <- file.path(run_root, "CombinedPDF")
heatmaps_root <- file.path(run_root, "Heatmaps")
logs_root <- file.path(run_root, "Logs")
stats_root <- file.path(run_root, "Stats")
tmp_root <- file.path(run_root, "tmp")

ensure_dir(run_root)
ensure_dir(cache_root)
ensure_dir(combined_root)
ensure_dir(heatmaps_root)
ensure_dir(logs_root)
ensure_dir(stats_root)
ensure_dir(tmp_root)

run_ts <- file_ts()
event_log_path <- file.path(
  logs_root,
  paste0("event_log_", cfg$run_id, "_", run_ts, ".csv")
)
init_event_log(event_log_path)

status_from_invalid_reason <- function(reason) {
  if (is.na(reason) || !nzchar(reason)) return("no_valid_kw_pvalue")
  if (grepl("Not enough groups|No rows after filtering|empty", reason, ignore.case = TRUE)) {
    return("no_valid_haplotype_grouping")
  }
  "no_valid_kw_pvalue"
}

empty_per_test_stats <- function() {
  data.table(
    trait = character(),
    gene_file = character(),
    gene_name = character(),
    title = character(),
    MGmin = integer(),
    epsilon = character(),
    test_signature = character(),
    n_ind = integer(),
    n_groups = integer(),
    group_sizes = character(),
    kw_p_raw = numeric(),
    kw_p_holm_within_gene = numeric(),
    is_best_within_gene = logical(),
    combined_pdf_path = character(),
    heatmap_pdf_path = character(),
    notes = character()
  )
}

empty_gene_summary <- function() {
  data.table(
    trait = character(),
    gene_file = character(),
    gene_name = character(),
    title = character(),
    n_total_valid_tests_before_collapse = integer(),
    n_unique_valid_tests = integer(),
    best_MGmin = integer(),
    best_epsilon = character(),
    best_n_ind = integer(),
    best_n_groups = integer(),
    best_group_sizes = character(),
    best_kw_p_raw = numeric(),
    best_kw_p_holm_within_gene = numeric(),
    trait_fdr_p = numeric(),
    trait_bonferroni_p = numeric(),
    significant_fdr = logical(),
    significant_bonferroni = logical()
  )
}

empty_successful_runs <- function() {
  data.table(
    trait = character(),
    gene_file = character(),
    gene_name = character(),
    title = character(),
    MGmin = integer(),
    cache_rds = character(),
    combined_pdf_path = character()
  )
}

log_msg(
  "Run started.",
  level = "INFO",
  stage = "run",
  event = "run_start",
  event_log_path = event_log_path
)

log_msg(
  paste0("Run ID: ", cfg$run_id),
  level = "INFO",
  stage = "run",
  event = "run_config",
  event_log_path = event_log_path
)

log_msg(
  paste0("Event log path: ", event_log_path),
  level = "INFO",
  stage = "run",
  event = "event_log_created",
  event_log_path = event_log_path
)

log_msg(
  paste0("Traits: ", paste(cfg$traits, collapse = ", ")),
  level = "INFO",
  stage = "run",
  event = "run_config",
  event_log_path = event_log_path
)

log_msg(
  paste0("MGmin values: ", paste(cfg$mgmin_values, collapse = ", ")),
  level = "INFO",
  stage = "run",
  event = "run_config",
  event_log_path = event_log_path
)

log_msg(
  paste0("Epsilon vector: ", paste(cfg$epsilon_vector, collapse = ", ")),
  level = "INFO",
  stage = "run",
  event = "run_config",
  event_log_path = event_log_path
)

targets_list <- lapply(cfg$traits, function(trait) {
  raw_dir <- make_raw_vcf_dir(cfg$vcf_root, trait)
  if (!dir.exists(raw_dir)) {
    stop("Raw VCF directory does not exist for trait ", trait, ": ", raw_dir)
  }
  list_trait_targets(cfg$vcf_root, trait)
})

targets_dt <- rbindlist(targets_list, fill = TRUE)
targets_dt <- unique(targets_dt[, .(trait, gene_file, title)])
setorder(targets_dt, trait, gene_file)

if (nrow(targets_dt) == 0) {
  stop("No .vcf.gz targets discovered in raw harmonized VCF directories.")
}


write_targets_yaml(as.data.frame(targets_dt), targets_path)

log_msg(
  paste0("Auto-generated targets file: ", targets_path),
  level = "INFO",
  stage = "discovery",
  event = "targets_written",
  event_log_path = event_log_path
)

for (i in seq_len(nrow(targets_dt))) {
  log_msg(
    paste0("Discovered target ", i, "/", nrow(targets_dt)),
    level = "INFO",
    trait = targets_dt$trait[i],
    gene_file = targets_dt$gene_file[i],
    stage = "discovery",
    event = "target_discovered",
    event_log_path = event_log_path
  )
}

failures_dt <- init_failures_df()
per_test_dt <- empty_per_test_stats()
gene_summary_dt <- empty_gene_summary()
successful_runs_dt <- empty_successful_runs()

for (i in seq_len(nrow(targets_dt))) {
  trait <- targets_dt$trait[i]
  gene_file <- targets_dt$gene_file[i]
  title <- targets_dt$title[i]
  gene_name <- make_gene_name(gene_file)

  log_msg(
    paste0("Processing gene ", i, "/", nrow(targets_dt)),
    level = "INFO",
    trait = trait,
    gene_file = gene_file,
    stage = "gene",
    event = "gene_start",
    event_log_path = event_log_path
  )

  gene_valid_dt <- empty_per_test_stats()[0]

  for (MGmin in cfg$mgmin_values) {
    log_msg(
      "Starting MGmin analysis.",
      level = "INFO",
      trait = trait,
      gene_file = gene_file,
      MGmin = MGmin,
      stage = "mgmin",
      event = "mgmin_start",
      event_log_path = event_log_path
    )

    cache_dir <- file.path(cache_root, trait, gene_name, paste0("MGmin_", MGmin))
    ensure_dir(cache_dir)

    cache_rds <- file.path(cache_dir, "HapObject.rds")
    res <- NULL

    if (file.exists(cache_rds) && !isTRUE(cfg$overwrite_cache)) {
      log_msg(
        paste0("Using cached HapObject: ", cache_rds),
        level = "INFO",
        trait = trait,
        gene_file = gene_file,
        MGmin = MGmin,
        stage = "cache",
        event = "cache_used",
        event_log_path = event_log_path
      )

      res <- tryCatch(
        readRDS(cache_rds),
        error = function(e) {
          msg <- paste0("Failed to read cache RDS: ", conditionMessage(e))
          failures_dt <<- add_failure(
            failures_dt,
            trait = trait,
            gene_file = gene_file,
            MGmin = MGmin,
            epsilon = NA_character_,
            stage = "cache",
            status = "cache_read_failed",
            message = msg
          )
          log_msg(
            msg,
            level = "ERROR",
            trait = trait,
            gene_file = gene_file,
            MGmin = MGmin,
            stage = "cache",
            event = "failure_encountered",
            event_log_path = event_log_path
          )
          NULL
        }
      )
    } else {
      log_msg(
        "Running CrossHap.",
        level = "INFO",
        trait = trait,
        gene_file = gene_file,
        MGmin = MGmin,
        stage = "crosshap",
        event = "crosshap_start",
        event_log_path = event_log_path
      )

      res <- tryCatch(
        run_crosshap_harmonized_wild(
          trait = trait,
          gene_file = gene_file,
          epsilon_vector = cfg$epsilon_vector,
          MGmin = MGmin,
          minHap = cfg$minHap,
          hetmiss_as = cfg$hetmiss_as,
          keep_outliers = cfg$keep_outliers,
          vcf_root = cfg$vcf_root,
          pheno_root = cfg$pheno_root,
          plink_bin = cfg$plink_bin,
          tmp_dir = tmp_root
        ),
        error = function(e) {
          msg <- paste0("CrossHap failed: ", conditionMessage(e))
          failures_dt <<- add_failure(
            failures_dt,
            trait = trait,
            gene_file = gene_file,
            MGmin = MGmin,
            epsilon = NA_character_,
            stage = "crosshap",
            status = "crosshap_failed",
            message = msg
          )
          log_msg(
            msg,
            level = "ERROR",
            trait = trait,
            gene_file = gene_file,
            MGmin = MGmin,
            stage = "crosshap",
            event = "failure_encountered",
            event_log_path = event_log_path
          )
          NULL
        }
      )

      if (!is.null(res)) {
        saveRDS(res, cache_rds)
        log_msg(
          paste0("Saved cache: ", cache_rds),
          level = "INFO",
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          stage = "cache",
          event = "cache_written",
          event_log_path = event_log_path
        )
      }
    }

    if (is.null(res)) {
      if (isTRUE(cfg$fail_fast) && !isTRUE(cfg$continue_on_error)) {
        stop("Stopping because fail_fast is TRUE and continue_on_error is FALSE.")
      }
      next
    }

    comb_dir <- file.path(combined_root, trait, gene_name, paste0("MGmin_", MGmin))
    ensure_dir(comb_dir)

    combined_pdf_path <- file.path(
      comb_dir,
      paste0(
        "CrosshapTree+Violin__",
        gene_name,
        "__MGmin", MGmin,
        "__EpsVec", length(cfg$epsilon_vector),
        ".pdf"
      )
    )

    successful_runs_dt <- rbind(
      successful_runs_dt,
      data.table(
        trait = trait,
        gene_file = gene_file,
        gene_name = gene_name,
        title = title,
        MGmin = as.integer(MGmin),
        cache_rds = cache_rds,
        combined_pdf_path = combined_pdf_path
      ),
      fill = TRUE
    )

    tryCatch(
      write_combined_pdf(
        HapObject = res$HapObject,
        out_pdf = combined_pdf_path,
        title = title,
        trait = trait,
        gene_file = gene_file,
        gene_name = gene_name,
        MGmin = MGmin,
        epsilon_vector = cfg$epsilon_vector,
        mgmin_test_stats = NULL,
        gene_summary_row = NULL
      ),
      error = function(e) {
        msg <- paste0("Combined PDF failed: ", conditionMessage(e))
        failures_dt <<- add_failure(
          failures_dt,
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = NA_character_,
          stage = "combined_pdf",
          status = "combined_pdf_failed",
          message = msg
        )
        log_msg(
          msg,
          level = "ERROR",
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          stage = "combined_pdf",
          event = "failure_encountered",
          event_log_path = event_log_path
        )
      }
    )

    if (file.exists(combined_pdf_path)) {
      log_msg(
        paste0("Combined PDF written: ", combined_pdf_path),
        level = "INFO",
        trait = trait,
        gene_file = gene_file,
        MGmin = MGmin,
        stage = "combined_pdf",
        event = "pdf_generated",
        event_log_path = event_log_path
      )
    }

    heatmap_dir <- file.path(heatmaps_root, trait, gene_name, paste0("MGmin_", MGmin))
    ensure_dir(heatmap_dir)

    for (eps in cfg$epsilon_vector) {
      label <- paste0("Haplotypes_MGmin", MGmin, "_E", eps)
      heatmap_pdf_path <- file.path(
        heatmap_dir,
        paste0(
          "Heatmaps__",
          gene_name,
          "__MGmin", MGmin,
          "__Eps", eps_tag(eps),
          ".pdf"
        )
      )

      if (is.null(res$HapObject[[label]])) {
        msg <- "No HapObject entry for this epsilon; excluded from stats."
        failures_dt <- add_failure(
          failures_dt,
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = as.character(eps),
          stage = "epsilon_validation",
          status = "no_valid_haplotype_grouping",
          message = msg
        )
        log_msg(
          msg,
          level = "WARN",
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = as.character(eps),
          stage = "epsilon_validation",
          event = "epsilon_processed",
          event_log_path = event_log_path
        )
        next
      }

      tryCatch(
        write_heatmaps_for_eps(
          HapObject = res$HapObject,
          label = label,
          gene_file = gene_file,
          title = title,
          trait = trait,
          eps = eps,
          MGmin = MGmin,
          raw_path = res$raw_path,
          common_ids = res$common_ids,
          vcf_raw_ids = res$vcf_raw_ids,
          out_pdf = heatmap_pdf_path
        ),
        error = function(e) {
          msg <- paste0("Heatmap failed: ", conditionMessage(e))
          failures_dt <<- add_failure(
            failures_dt,
            trait = trait,
            gene_file = gene_file,
            MGmin = MGmin,
            epsilon = as.character(eps),
            stage = "heatmap",
            status = "heatmap_failed",
            message = msg
          )
          log_msg(
            msg,
            level = "ERROR",
            trait = trait,
            gene_file = gene_file,
            MGmin = MGmin,
            epsilon = as.character(eps),
            stage = "heatmap",
            event = "failure_encountered",
            event_log_path = event_log_path
          )
        }
      )

      if (file.exists(heatmap_pdf_path)) {
        log_msg(
          paste0("Heatmap PDF written: ", heatmap_pdf_path),
          level = "INFO",
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = as.character(eps),
          stage = "heatmap",
          event = "heatmap_generated",
          event_log_path = event_log_path
        )
      }

      kw_info <- compute_valid_kw_test(res$HapObject, label)

      if (!isTRUE(kw_info$valid)) {
        status <- status_from_invalid_reason(kw_info$reason)
        failures_dt <- add_failure(
          failures_dt,
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = as.character(eps),
          stage = "kw_validation",
          status = status,
          message = kw_info$reason
        )
        log_msg(
          paste0("Invalid omnibus test excluded from stats: ", kw_info$reason),
          level = "WARN",
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = as.character(eps),
          stage = "kw_validation",
          event = "epsilon_processed",
          event_log_path = event_log_path
        )
        next
      }

      test_signature <- paste(
        formatC(kw_info$kw_p_raw, digits = 15, format = "fg", flag = "#"),
        kw_info$n_groups,
        kw_info$sorted_group_sizes,
        sep = "__"
      )

      gene_valid_dt <- rbind(
        gene_valid_dt,
        data.table(
          trait = trait,
          gene_file = gene_file,
          gene_name = gene_name,
          title = title,
          MGmin = as.integer(MGmin),
          epsilon = as.character(eps),
          test_signature = test_signature,
          n_ind = as.integer(kw_info$n_ind),
          n_groups = as.integer(kw_info$n_groups),
          group_sizes = kw_info$group_sizes,
          kw_p_raw = as.numeric(kw_info$kw_p_raw),
          kw_p_holm_within_gene = NA_real_,
          is_best_within_gene = FALSE,
          combined_pdf_path = combined_pdf_path,
          heatmap_pdf_path = heatmap_pdf_path,
          notes = ""
        ),
        fill = TRUE
      )

      log_msg(
        paste0(
          "Valid KW test found. p=",
          format_p_plain(kw_info$kw_p_raw),
          " | groups=",
          kw_info$n_groups,
          " | sizes=",
          kw_info$group_sizes
        ),
        level = "INFO",
        trait = trait,
        gene_file = gene_file,
        MGmin = MGmin,
        epsilon = as.character(eps),
        stage = "kw_validation",
        event = "valid_kw_test_found",
        event_log_path = event_log_path
      )
    }
  }

  n_total_valid_before_collapse <- nrow(gene_valid_dt)

  if (n_total_valid_before_collapse == 0) {
    msg <- "Gene produced zero valid omnibus tests across all MGmin and epsilon values."
    failures_dt <- add_failure(
      failures_dt,
      trait = trait,
      gene_file = gene_file,
      MGmin = NA_integer_,
      epsilon = NA_character_,
      stage = "gene_stats",
      status = "no_valid_tests_for_gene",
      message = msg
    )
    log_msg(
      msg,
      level = "WARN",
      trait = trait,
      gene_file = gene_file,
      stage = "gene_stats",
      event = "failure_encountered",
      event_log_path = event_log_path
    )
    next
  }

  gene_valid_dt[, epsilon_num := as_numeric_epsilon(epsilon)]
  setorder(gene_valid_dt, epsilon_num, -MGmin, test_signature)

  unique_idx <- !duplicated(gene_valid_dt$test_signature)
  gene_unique_dt <- copy(gene_valid_dt[unique_idx])

  gene_valid_dt[, epsilon_num := NULL]
  gene_unique_dt[, epsilon_num := NULL]

  dup_counts <- gene_valid_dt[, .N, by = test_signature]
  gene_unique_dt <- merge(gene_unique_dt, dup_counts, by = "test_signature", all.x = TRUE, sort = FALSE)
  setnames(gene_unique_dt, "N", "duplicate_count")

  for (k in seq_len(nrow(gene_unique_dt))) {
    sig <- gene_unique_dt$test_signature[k]
    dup_rows <- gene_valid_dt[test_signature == sig]

    if (nrow(dup_rows) > 1) {
      dup_desc <- paste0("MGmin", dup_rows$MGmin, ":eps", dup_rows$epsilon, collapse = ";")
      gene_unique_dt$notes[k] <- paste0(
        "Duplicate signature retained by smallest epsilon, then larger MGmin. Collapsed from ",
        dup_desc
      )

      log_msg(
        paste0(
          "Collapsed duplicate tests using explicit rule: smallest epsilon, then larger MGmin. Signature=",
          sig,
          " | rows=",
          dup_desc
        ),
        level = "INFO",
        trait = trait,
        gene_file = gene_file,
        MGmin = gene_unique_dt$MGmin[k],
        epsilon = gene_unique_dt$epsilon[k],
        stage = "collapse",
        event = "duplicate_collapsed",
        event_log_path = event_log_path
      )
    }
  }

  gene_unique_dt$kw_p_holm_within_gene <- safe_p_adjust(
    gene_unique_dt$kw_p_raw,
    method = cfg$stats$within_gene_method
  )

  log_msg(
    paste0(
      "Within-gene Holm correction complete. Valid before collapse=",
      n_total_valid_before_collapse,
      " | unique after collapse=",
      nrow(gene_unique_dt)
    ),
    level = "INFO",
    trait = trait,
    gene_file = gene_file,
    stage = "within_gene_correction",
    event = "holm_correction_complete",
    event_log_path = event_log_path
  )

  best_holm_p <- min(gene_unique_dt$kw_p_holm_within_gene, na.rm = TRUE)
  best_idx <- which(gene_unique_dt$kw_p_holm_within_gene == best_holm_p)

  if (length(best_idx) != 1) {
    tied_dt <- copy(gene_unique_dt[best_idx])

    tied_dt[, epsilon_num_tie := as_numeric_epsilon(epsilon)]
    setorder(tied_dt, kw_p_raw, epsilon_num_tie, -MGmin)

    chosen_row <- tied_dt[1]

    tie_msg <- paste0(
      "AMBIGUOUS_BEST_TEST_RESOLVED: trait=",
      trait,
      " | gene_file=",
      gene_file,
      " | matching_holm_p=",
      best_holm_p,
      " | chosen_MGmin=",
      chosen_row$MGmin,
      " | chosen_epsilon=",
      chosen_row$epsilon,
      " | chosen_raw_p=",
      chosen_row$kw_p_raw,
      " | tie_break_rule=smallest_raw_p_then_smallest_epsilon_then_larger_MGmin"
    )

    log_msg(
      tie_msg,
      level = "WARN",
      trait = trait,
      gene_file = gene_file,
      MGmin = chosen_row$MGmin,
      epsilon = chosen_row$epsilon,
      stage = "within_gene_summary",
      event = "ambiguous_best_test_resolved",
      event_log_path = event_log_path
    )

    best_idx <- which(
      gene_unique_dt$MGmin == chosen_row$MGmin &
        gene_unique_dt$epsilon == chosen_row$epsilon &
        gene_unique_dt$kw_p_raw == chosen_row$kw_p_raw &
        gene_unique_dt$kw_p_holm_within_gene == chosen_row$kw_p_holm_within_gene
    )[1]
  }

  gene_unique_dt$is_best_within_gene[best_idx] <- TRUE
  gene_unique_dt$duplicate_count <- NULL

  per_test_dt <- rbind(per_test_dt, gene_unique_dt, fill = TRUE)

  best_row <- gene_unique_dt[best_idx]

  gene_summary_dt <- rbind(
    gene_summary_dt,
    data.table(
      trait = trait,
      gene_file = gene_file,
      gene_name = gene_name,
      title = title,
      n_total_valid_tests_before_collapse = as.integer(n_total_valid_before_collapse),
      n_unique_valid_tests = as.integer(nrow(gene_unique_dt)),
      best_MGmin = as.integer(best_row$MGmin),
      best_epsilon = as.character(best_row$epsilon),
      best_n_ind = as.integer(best_row$n_ind),
      best_n_groups = as.integer(best_row$n_groups),
      best_group_sizes = as.character(best_row$group_sizes),
      best_kw_p_raw = as.numeric(best_row$kw_p_raw),
      best_kw_p_holm_within_gene = as.numeric(best_row$kw_p_holm_within_gene),
      trait_fdr_p = NA_real_,
      trait_bonferroni_p = NA_real_,
      significant_fdr = NA,
      significant_bonferroni = NA
    ),
    fill = TRUE
  )

  log_msg(
    paste0(
      "Selected representative test: MGmin=",
      best_row$MGmin,
      " | epsilon=",
      best_row$epsilon,
      " | raw p=",
      format_p_plain(best_row$kw_p_raw),
      " | Holm p=",
      format_p_plain(best_row$kw_p_holm_within_gene)
    ),
    level = "INFO",
    trait = trait,
    gene_file = gene_file,
    MGmin = best_row$MGmin,
    epsilon = best_row$epsilon,
    stage = "within_gene_summary",
    event = "best_test_selected",
    event_log_path = event_log_path
  )
}

if (nrow(gene_summary_dt) > 0) {
  gene_summary_dt$trait_fdr_p <- NA_real_
  gene_summary_dt$trait_bonferroni_p <- NA_real_
  gene_summary_dt$significant_fdr <- FALSE
  gene_summary_dt$significant_bonferroni <- FALSE

  for (trait_name in unique(gene_summary_dt$trait)) {
    idx <- which(gene_summary_dt$trait == trait_name)

    fdr_vals <- safe_p_adjust(
      gene_summary_dt$best_kw_p_holm_within_gene[idx],
      method = cfg$stats$across_gene_fdr_method
    )

    bonf_vals <- safe_p_adjust(
      gene_summary_dt$best_kw_p_holm_within_gene[idx],
      method = "bonferroni"
    )

    gene_summary_dt$trait_fdr_p[idx] <- fdr_vals
    gene_summary_dt$trait_bonferroni_p[idx] <- bonf_vals
    gene_summary_dt$significant_fdr[idx] <- fdr_vals <= cfg$stats$alpha
    gene_summary_dt$significant_bonferroni[idx] <- bonf_vals <= cfg$stats$alpha

    log_msg(
      paste0(
        "Trait-level correction complete. Genes=",
        length(idx),
        " | FDR method=",
        cfg$stats$across_gene_fdr_method,
        " | alpha=",
        cfg$stats$alpha
      ),
      level = "INFO",
      trait = trait_name,
      stage = "trait_correction",
      event = "trait_level_correction_complete",
      event_log_path = event_log_path
    )
  }
}

if (nrow(successful_runs_dt) > 0) {
  successful_runs_dt <- unique(successful_runs_dt)
  setorder(successful_runs_dt, trait, gene_file, MGmin)

  for (i in seq_len(nrow(successful_runs_dt))) {
    trait <- successful_runs_dt$trait[i]
    gene_file <- successful_runs_dt$gene_file[i]
    gene_name <- successful_runs_dt$gene_name[i]
    title <- successful_runs_dt$title[i]
    MGmin <- successful_runs_dt$MGmin[i]
    cache_rds <- successful_runs_dt$cache_rds[i]
    combined_pdf_path <- successful_runs_dt$combined_pdf_path[i]

    res <- tryCatch(
      readRDS(cache_rds),
      error = function(e) {
        msg <- paste0("Failed to read cache during corrected PDF regeneration: ", conditionMessage(e))
        failures_dt <<- add_failure(
          failures_dt,
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = NA_character_,
          stage = "combined_pdf_postprocess",
          status = "cache_read_failed",
          message = msg
        )
        log_msg(
          msg,
          level = "ERROR",
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          stage = "combined_pdf_postprocess",
          event = "failure_encountered",
          event_log_path = event_log_path
        )
        NULL
      }
    )

    if (is.null(res)) next

    current_trait <- trait
    current_gene_file <- gene_file
    current_MGmin <- MGmin

    mgmin_test_stats <- as.data.frame(
      per_test_dt[
        trait == current_trait &
          gene_file == current_gene_file &
          MGmin == current_MGmin
      ],
      stringsAsFactors = FALSE
    )

    gene_summary_row <- as.data.frame(
      gene_summary_dt[
        trait == current_trait &
          gene_file == current_gene_file
      ],
      stringsAsFactors = FALSE
    )

    tryCatch(
      write_combined_pdf(
        HapObject = res$HapObject,
        out_pdf = combined_pdf_path,
        title = title,
        trait = trait,
        gene_file = gene_file,
        gene_name = gene_name,
        MGmin = MGmin,
        epsilon_vector = cfg$epsilon_vector,
        mgmin_test_stats = mgmin_test_stats,
        gene_summary_row = gene_summary_row
      ),
      error = function(e) {
        msg <- paste0("Corrected PDF regeneration failed: ", conditionMessage(e))
        failures_dt <<- add_failure(
          failures_dt,
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          epsilon = NA_character_,
          stage = "combined_pdf_postprocess",
          status = "combined_pdf_postprocess_failed",
          message = msg
        )
        log_msg(
          msg,
          level = "ERROR",
          trait = trait,
          gene_file = gene_file,
          MGmin = MGmin,
          stage = "combined_pdf_postprocess",
          event = "failure_encountered",
          event_log_path = event_log_path
        )
      }
    )

    if (file.exists(combined_pdf_path)) {
      log_msg(
        paste0("Combined PDF regenerated with corrected values: ", combined_pdf_path),
        level = "INFO",
        trait = trait,
        gene_file = gene_file,
        MGmin = MGmin,
        stage = "combined_pdf_postprocess",
        event = "pdf_generated",
        event_log_path = event_log_path
      )
    }
  }
}

setorder(per_test_dt, trait, gene_file, MGmin, epsilon)
setorder(gene_summary_dt, trait, best_kw_p_holm_within_gene, best_kw_p_raw)
setorder(failures_dt, trait, gene_file, MGmin, epsilon, stage)

per_test_path <- file.path(stats_root, "per_test_stats.csv")
gene_summary_path <- file.path(stats_root, "gene_summary.csv")
failures_path <- file.path(stats_root, "failures_notes.csv")

fwrite(per_test_dt, per_test_path)
fwrite(gene_summary_dt, gene_summary_path)
fwrite(failures_dt, failures_path)

log_msg(
  paste0("Wrote stats table: ", per_test_path),
  level = "INFO",
  stage = "output",
  event = "stats_written",
  event_log_path = event_log_path
)

log_msg(
  paste0("Wrote stats table: ", gene_summary_path),
  level = "INFO",
  stage = "output",
  event = "stats_written",
  event_log_path = event_log_path
)

log_msg(
  paste0("Wrote stats table: ", failures_path),
  level = "INFO",
  stage = "output",
  event = "stats_written",
  event_log_path = event_log_path
)

log_msg(
  paste0(
    "Run finished. Targets=",
    nrow(targets_dt),
    " | per_test_rows=",
    nrow(per_test_dt),
    " | gene_summary_rows=",
    nrow(gene_summary_dt),
    " | failures_rows=",
    nrow(failures_dt)
  ),
  level = "INFO",
  stage = "run",
  event = "run_end",
  event_log_path = event_log_path
)