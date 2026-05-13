#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

keep <- NULL
older_than_days <- NULL
output_root <- NULL
failed_only <- FALSE
dry_run <- FALSE

for (i in seq_along(args)) {
  if (args[[i]] == "--keep" && i < length(args)) {
    keep <- as.integer(args[[i + 1L]])
  }
  if (args[[i]] == "--older-than" && i < length(args)) {
    older_than_days <- as.integer(args[[i + 1L]])
  }
  if (args[[i]] == "--output-root" && i < length(args)) {
    output_root <- args[[i + 1L]]
  }
  if (args[[i]] == "--failed-only") {
    failed_only <- TRUE
  }
  if (args[[i]] == "--dry-run") {
    dry_run <- TRUE
  }
}

script_flag <- grep("^--file=", commandArgs(), value = TRUE)[1]
if (is.na(script_flag)) {
  stop("Unable to resolve script path.", call. = FALSE)
}
script_dir <- normalizePath(dirname(sub("^--file=", "", script_flag)), mustWork = TRUE)
options(tcga_toolkit.root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
if (!is.null(output_root)) {
  options(tcga_toolkit.output_root = normalizePath(output_root, mustWork = TRUE))
}
source(file.path(script_dir, "common.R"))

runs_root <- file.path(output_root(), "tcga_runs")
if (!dir.exists(runs_root)) {
  message("No tcga_runs directory found.")
  quit(status = 0)
}

dirs <- list.dirs(runs_root, recursive = FALSE, full.names = FALSE)
if (!length(dirs)) {
  message("No run directories to clean up.")
  quit(status = 0)
}

run_infos <- lapply(dirs, function(d) {
  meta_path <- file.path(runs_root, d, "run_metadata.json")
  info <- list(
    dir = d,
    path = file.path(runs_root, d),
    mtime = file.info(file.path(runs_root, d))$mtime,
    has_meta = file.exists(meta_path)
  )
  if (info$has_meta) {
    meta <- tryCatch(jsonlite::fromJSON(meta_path, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(meta)) {
      info$task <- meta$task %||% NA_character_
      info$status <- meta$status %||% "unknown"
      info$completed_at <- meta$completed_at %||% NA_character_
    }
  }
  info
})

# Build a data frame for easier filtering
df <- data.frame(
  dir = sapply(run_infos, `[[`, "dir"),
  path = sapply(run_infos, `[[`, "path"),
  mtime = as.POSIXct(sapply(run_infos, `[[`, "mtime"), origin = "1970-01-01"),
  status = sapply(run_infos, function(x) x$status %||% "unknown"),
  stringsAsFactors = FALSE
)

to_delete <- rep(FALSE, nrow(df))
reasons <- rep(NA_character_, nrow(df))

# --failed-only
if (failed_only) {
  mask <- df$status == "failed"
  to_delete <- to_delete | mask
  reasons[mask] <- "failed run"
}

# --older-than
if (!is.null(older_than_days)) {
  cutoff <- Sys.time() - older_than_days * 86400
  mask <- df$mtime < cutoff
  to_delete <- to_delete | mask
  reasons[mask] <- ifelse(is.na(reasons[mask]), sprintf("older than %s days", older_than_days), paste(reasons[mask], sprintf("older than %s days", older_than_days), sep = "; "))
}

# --keep N (most recent by mtime)
if (!is.null(keep)) {
  sorted <- order(df$mtime, decreasing = TRUE)
  if (keep < nrow(df)) {
    drop_idx <- sorted[(keep + 1):nrow(df)]
    to_delete[drop_idx] <- TRUE
    reasons[drop_idx] <- ifelse(is.na(reasons[drop_idx]), sprintf("beyond top %s recent", keep), paste(reasons[drop_idx], sprintf("beyond top %s recent", keep), sep = "; "))
  }
}

# If no filter flags were provided, show usage and exit
if (is.null(keep) && is.null(older_than_days) && !failed_only) {
  cat("Usage: Rscript tcga_toolkit/scripts/cleanup_runs.R [OPTIONS]\n")
  cat("\nOptions:\n")
  cat("  --keep N              Keep the N most recent runs, delete the rest\n")
  cat("  --older-than DAYS     Delete runs older than DAYS days\n")
  cat("  --failed-only         Delete only runs with status 'failed'\n")
  cat("  --output-root DIR     Target a different output root (default: project root)\n")
  cat("  --dry-run             Print what would be deleted without removing anything\n")
  cat("\nExamples:\n")
  cat("  Rscript tcga_toolkit/scripts/cleanup_runs.R --keep 20\n")
  cat("  Rscript tcga_toolkit/scripts/cleanup_runs.R --older-than 30 --failed-only\n")
  cat("  Rscript tcga_toolkit/scripts/cleanup_runs.R --keep 10 --dry-run\n")
  quit(status = 0)
}

if (!any(to_delete)) {
  message("No run directories matched the cleanup criteria.")
  quit(status = 0)
}

del_df <- df[to_delete, ]
del_df$reason <- reasons[to_delete]

if (dry_run) {
  cat("[DRY RUN] The following directories would be deleted:\n")
  for (i in seq_len(nrow(del_df))) {
    cat(sprintf("  - %s (%s)\n", del_df$dir[i], del_df$reason[i]))
  }
  cat(sprintf("\nTotal: %s directories\n", nrow(del_df)))
} else {
  cat("Deleting matched run directories:\n")
  for (i in seq_len(nrow(del_df))) {
    cat(sprintf("  - %s (%s)\n", del_df$dir[i], del_df$reason[i]))
    unlink(del_df$path[i], recursive = TRUE)
  }
  cat(sprintf("\nDeleted: %s directories\n", nrow(del_df)))
}
