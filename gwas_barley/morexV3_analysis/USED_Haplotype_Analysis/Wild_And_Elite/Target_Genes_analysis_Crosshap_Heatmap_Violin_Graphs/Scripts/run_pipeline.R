suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

argv <- commandArgs(trailingOnly = FALSE)
file_arg <- argv[grepl("^--file=", argv)]
if (length(file_arg) == 0) stop("Run with: Rscript run_pipeline.R")

script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
script_dir <- dirname(script_path)

source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "run_crosshap_intersection_wild_elite.R"))
source(file.path(script_dir, "R", "plot_combined_pdf.R"))
source(file.path(script_dir, "R", "plot_heatmaps.R"))

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path(script_dir, "config.yaml")
targets_path <- if (length(args) >= 2) args[2] else file.path(script_dir, "targets.yaml")

cfg <- yaml::read_yaml(config_path)
cfg$plink_bin <- cfg$plink_bin %||% "plink"
cfg$minHap <- cfg$minHap %||% 9
cfg$hetmiss_as <- cfg$hetmiss_as %||% "allele"
cfg$keep_outliers <- cfg$keep_outliers %||% FALSE
cfg$continue_on_error <- cfg$continue_on_error %||% TRUE
cfg$overwrite_cache <- cfg$overwrite_cache %||% FALSE
cfg$fail_fast <- cfg$fail_fast %||% FALSE
cfg$stats <- cfg$stats %||% list()
cfg$stats$pairwise_method <- cfg$stats$pairwise_method %||% "holm"

targets_yaml <- yaml::read_yaml(targets_path)
targets_list <- lapply(targets_yaml$targets, normalize_target)
if (length(targets_list) == 0) stop("No targets found in targets.yaml")

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

event_log_path <- file.path(logs_root, paste0("event_log_", cfg$run_id, "_", file_ts(), ".csv"))
init_event_log(event_log_path)

log_msg("Run started.", stage = "run", event = "run_start", event_log_path = event_log_path)
log_msg(paste0("Run ID: ", cfg$run_id), stage = "run", event = "run_config", event_log_path = event_log_path)
log_msg(paste0("Targets: ", length(targets_list)), stage = "run", event = "run_config", event_log_path = event_log_path)

failures_dt <- init_failures_df()
summary_dt <- data.table(
  trait = character(),
  gene_file = character(),
  gene_name = character(),
  title = character(),
  epsilon = numeric(),
  MGmin = integer(),
  kw_p_raw = numeric(),
  n_ind = integer(),
  n_groups = integer(),
  group_sizes = character(),
  combined_pdf_path = character(),
  heatmap_pdf_path = character(),
  elite_assignments_csv = character(),
  cache_rds = character(),
  status = character(),
  note = character()
)

for (i in seq_along(targets_list)) {
  target <- targets_list[[i]]
  gene_name <- make_gene_name(target$gene_file)

  log_msg(
    paste0("Processing target ", i, "/", length(targets_list)),
    trait = target$trait,
    gene_file = target$gene_file,
    stage = "gene",
    event = "gene_start",
    event_log_path = event_log_path
  )

  cache_dir <- file.path(cache_root, target$trait, gene_name)
  ensure_dir(cache_dir)
  cache_rds <- file.path(cache_dir, "HapObject.rds")
  context_rds <- file.path(cache_dir, "RunContext.rds")

  combined_pdf_path <- file.path(combined_root, target$trait, paste0(gene_name, "__MGmin", target$MGmin, "__E", target$epsilon, ".pdf"))
  heatmap_pdf_path <- file.path(heatmaps_root, target$trait, paste0(gene_name, "__MGmin", target$MGmin, "__E", target$epsilon, "__heatmaps.pdf"))
  elite_csv_path <- file.path(combined_root, target$trait, paste0(gene_name, "__MGmin", target$MGmin, "__E", target$epsilon, "__elite_assignments.csv"))

  target_status <- "success"
  target_note <- ""

  tryCatch({
    if (file.exists(cache_rds) && file.exists(context_rds) && !isTRUE(cfg$overwrite_cache)) {
      HapObject <- readRDS(cache_rds)
      run_context <- readRDS(context_rds)
      run_context$HapObject <- HapObject
      log_msg("Using cached HapObject.", trait = target$trait, gene_file = target$gene_file, stage = "cache", event = "cache_hit", event_log_path = event_log_path)
    } else {
      run_context <- run_crosshap_intersection_wild_elite(target, cfg, tmp_root)
      HapObject <- run_context$HapObject
      saveRDS(HapObject, cache_rds)
      saveRDS(run_context[setdiff(names(run_context), "HapObject")], context_rds)
      log_msg("Created and cached HapObject.", trait = target$trait, gene_file = target$gene_file, stage = "cache", event = "cache_write", event_log_path = event_log_path)
    }

    pdf_info <- write_combined_pdf(
      HapObject = run_context$HapObject,
      target = target,
      pheno = run_context$pheno,
      out_pdf = combined_pdf_path,
      pairwise_method = cfg$stats$pairwise_method
    )

    ensure_dir(dirname(elite_csv_path))
    utils::write.csv(pdf_info$elite_assignments, elite_csv_path, row.names = FALSE)

    write_heatmaps_for_target(
      HapObject = run_context$HapObject,
      target = target,
      merged_path = run_context$merged_path,
      common_pos = run_context$common_pos,
      vcf_merged_ids = run_context$vcf_merged_ids,
      out_pdf = heatmap_pdf_path
    )

    summary_dt <- rbind(
      summary_dt,
      data.table(
        trait = target$trait,
        gene_file = target$gene_file,
        gene_name = gene_name,
        title = target$title,
        epsilon = target$epsilon,
        MGmin = target$MGmin,
        kw_p_raw = pdf_info$kw_p,
        n_ind = pdf_info$n_ind,
        n_groups = pdf_info$n_groups,
        group_sizes = pdf_info$group_sizes,
        combined_pdf_path = combined_pdf_path,
        heatmap_pdf_path = heatmap_pdf_path,
        elite_assignments_csv = elite_csv_path,
        cache_rds = cache_rds,
        status = target_status,
        note = target_note
      ),
      fill = TRUE
    )

    log_msg("Target completed.", trait = target$trait, gene_file = target$gene_file, stage = "gene", event = "gene_done", event_log_path = event_log_path)
  }, error = function(e) {
    target_status <<- "failed"
    target_note <<- conditionMessage(e)
    failures_dt <<- add_failure(failures_dt, target$trait, target$gene_file, stage = "gene", status = "failed", message = target_note)

    summary_dt <<- rbind(
      summary_dt,
      data.table(
        trait = target$trait,
        gene_file = target$gene_file,
        gene_name = gene_name,
        title = target$title,
        epsilon = target$epsilon,
        MGmin = target$MGmin,
        kw_p_raw = NA_real_,
        n_ind = NA_integer_,
        n_groups = NA_integer_,
        group_sizes = NA_character_,
        combined_pdf_path = combined_pdf_path,
        heatmap_pdf_path = heatmap_pdf_path,
        elite_assignments_csv = elite_csv_path,
        cache_rds = cache_rds,
        status = target_status,
        note = target_note
      ),
      fill = TRUE
    )

    log_msg(target_note, level = "ERROR", trait = target$trait, gene_file = target$gene_file, stage = "gene", event = "gene_failed", event_log_path = event_log_path)
    if (!isTRUE(cfg$continue_on_error) || isTRUE(cfg$fail_fast)) stop(e)
  })
}

fwrite(summary_dt, file.path(stats_root, "per_gene_summary.csv"))
fwrite(failures_dt, file.path(stats_root, "failures_notes.csv"))

log_msg("Run finished.", stage = "run", event = "run_done", event_log_path = event_log_path)