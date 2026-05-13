#!/usr/bin/env Rscript

# list_runs.R
# Summarise tcga_runs/ contents so an agent can find prior work without
# walking the directory manually.
#
# Usage:
#   Rscript tcga_toolkit/scripts/list_runs.R
#   Rscript tcga_toolkit/scripts/list_runs.R --status failed
#   Rscript tcga_toolkit/scripts/list_runs.R --task run_deg
#   Rscript tcga_toolkit/scripts/list_runs.R --since 2026-05-01
#   Rscript tcga_toolkit/scripts/list_runs.R --json --limit 50

args <- commandArgs(trailingOnly = TRUE)
status_filter <- NULL
task_filter <- NULL
since_filter <- NULL
limit <- NA_integer_
as_json <- FALSE
output_root_override <- NULL
for (i in seq_along(args)) {
  a <- args[[i]]
  if (a == "--status" && i < length(args)) status_filter <- args[[i + 1L]]
  if (a == "--task" && i < length(args)) task_filter <- args[[i + 1L]]
  if (a == "--since" && i < length(args)) since_filter <- args[[i + 1L]]
  if (a == "--limit" && i < length(args)) limit <- suppressWarnings(as.integer(args[[i + 1L]]))
  if (a == "--output-root" && i < length(args)) output_root_override <- args[[i + 1L]]
  if (a == "--json") as_json <- TRUE
}

script_flag <- grep("^--file=", commandArgs(), value = TRUE)[1]
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_flag), perl = TRUE)
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
options(tcga_toolkit.root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
if (!is.null(output_root_override)) {
  options(tcga_toolkit.output_root = normalizePath(output_root_override, mustWork = TRUE))
}
source(file.path(script_dir, "common.R"))

runs_dir <- file.path(output_root(), "tcga_runs")
if (!dir.exists(runs_dir)) {
  if (as_json) cat("[]\n") else cat(sprintf("No tcga_runs/ at %s\n", runs_dir))
  quit(status = 0)
}

dirs <- list.dirs(runs_dir, recursive = FALSE, full.names = TRUE)

since_time <- NULL
if (!is.null(since_filter)) {
  since_time <- tryCatch(as.POSIXct(since_filter, tz = "UTC"),
                         error = function(e) NULL)
  if (is.null(since_time) || is.na(since_time)) {
    stop(sprintf("Invalid --since value: %s (expected YYYY-MM-DD)", since_filter), call. = FALSE)
  }
}

rows <- list()
for (d in dirs) {
  meta_path <- file.path(d, "run_metadata.json")
  meta <- if (file.exists(meta_path)) {
    tryCatch(jsonlite::fromJSON(meta_path, simplifyVector = FALSE),
             error = function(e) NULL)
  } else NULL

  mtime <- file.info(d)$mtime
  results_dir <- file.path(d, "results")
  plots_dir <- file.path(d, "plots")

  row <- list(
    run_dir = d,
    task_id = basename(d),
    task = meta$task %||% NA,
    toolkit_version = meta$toolkit_version %||% NA,
    status = meta$status %||% NA,
    completed_at = meta$completed_at %||% NA,
    modified_at = format(mtime, "%Y-%m-%d %H:%M:%S"),
    n_result_files = if (dir.exists(results_dir)) length(list.files(results_dir)) else 0L,
    n_plot_files = if (dir.exists(plots_dir)) length(list.files(plots_dir)) else 0L,
    error = meta$error %||% NA
  )

  if (!is.null(status_filter) && !identical(row$status, status_filter)) next
  if (!is.null(task_filter) && !identical(row$task, task_filter)) next
  if (!is.null(since_time) && !is.na(mtime) && mtime < since_time) next

  rows[[length(rows) + 1L]] <- row
}

rows <- rows[order(vapply(rows, function(r) r$modified_at %||% "", character(1)), decreasing = TRUE)]
if (!is.na(limit) && limit > 0 && length(rows) > limit) {
  rows <- rows[seq_len(limit)]
}

if (as_json) {
  cat(jsonlite::toJSON(rows, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")
} else {
  if (!length(rows)) {
    cat("No matching runs.\n")
    quit(status = 0)
  }
  cat(sprintf("%-46s %-26s %-9s %-19s files\n", "task_id", "task", "status", "completed_at"))
  cat(strrep("-", 110), "\n", sep = "")
  for (r in rows) {
    cat(sprintf("%-46s %-26s %-9s %-19s r=%s p=%s\n",
                substr(r$task_id, 1, 46),
                substr(r$task %||% "?", 1, 26),
                r$status %||% "?",
                r$completed_at %||% "?",
                r$n_result_files,
                r$n_plot_files))
    if (!is.na(r$error) && !is.null(r$error)) {
      cat(sprintf("    error: %s\n", substr(as.character(r$error), 1, 140)))
    }
  }
}
