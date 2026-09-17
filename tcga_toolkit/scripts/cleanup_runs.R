#!/usr/bin/env Rscript

# Preview by default. Applying moves complete native run directories into a
# recoverable archive; it never unlinks run data. All supplied filters are AND.
args <- commandArgs(trailingOnly = TRUE)
keep <- NULL
older_than_days <- NULL
root_arg <- NULL
protected_file <- NULL
archive_plan_file <- NULL
failed_only <- FALSE
apply <- FALSE
dry_run <- FALSE
i <- 1L
while (i <= length(args)) {
  flag <- args[[i]]
  if (flag %in% c("--keep", "--older-than", "--output-root", "--protected-runs", "--archive-plan")) {
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) stop("Missing value for ", flag, call. = FALSE)
    value <- args[[i + 1L]]
    if (flag %in% c("--keep", "--older-than") && !grepl("^[0-9]+$", value)) stop(flag, " requires a nonnegative integer.", call. = FALSE)
    if (flag %in% c("--keep", "--older-than") && !is.finite(as.numeric(value))) stop(flag, " requires a finite integer.", call. = FALSE)
    if (flag == "--keep") keep <- as.numeric(value)
    if (flag == "--older-than") older_than_days <- as.numeric(value)
    if (flag == "--output-root") root_arg <- value
    if (flag == "--protected-runs") protected_file <- value
    if (flag == "--archive-plan") archive_plan_file <- value
    i <- i + 2L
  } else {
    if (!flag %in% c("--failed-only", "--apply", "--dry-run", "--help")) stop("Unknown option: ", flag, call. = FALSE)
    if (flag == "--failed-only") failed_only <- TRUE
    if (flag == "--apply") apply <- TRUE
    if (flag == "--dry-run") dry_run <- TRUE
    i <- i + 1L
  }
}
if (apply && dry_run) stop("Use either --apply or --dry-run, not both.", call. = FALSE)
if ("--help" %in% args || (is.null(keep) && is.null(older_than_days) && !failed_only)) {
  cat(paste0("Usage: Rscript tcga_toolkit/scripts/cleanup_runs.R [filters] [options]\n",
             "Filters combine with AND:\n",
             "  --keep N              Consider runs outside the N most recent completed runs\n",
             "  --older-than DAYS     Consider runs whose directory mtime is older than DAYS\n",
             "  --failed-only         Consider only failed runs\n",
             "Options:\n",
             "  --output-root DIR     Project containing tcga_runs/\n",
             "  --protected-runs FILE Additional protected native run IDs; # comments allowed\n",
             "  --archive-plan FILE   Fresh JSON plan from manage-research-project lifecycle helper\n",
             "  --dry-run             Preview only (the default)\n",
             "  --apply               Move allowlisted candidates into archive/tcga_runs/; requires archive-plan\n",
             "Only candidate_native_directories from a verified plan authorize moving runs.\n",
             "Protected IDs can add protection but cannot authorize an archive.\n"))
  quit(status = 0)
}
if (apply && is.null(archive_plan_file)) stop("--apply requires --archive-plan FILE from the project lifecycle helper.", call. = FALSE)

script_flag <- grep("^--file=", commandArgs(), value = TRUE)[1L]
if (is.na(script_flag)) stop("Unable to resolve script path.", call. = FALSE)
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_flag), perl = TRUE)
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
options(tcga_toolkit.root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
if (!is.null(root_arg)) options(tcga_toolkit.output_root = normalizePath(root_arg, mustWork = TRUE))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "archive_plan.R"))
verified_plan <- if (!is.null(archive_plan_file)) validate_archive_plan(archive_plan_file, output_root()) else NULL
runs_root <- file.path(output_root(), "tcga_runs")
if (!dir.exists(runs_root)) {
  message("No tcga_runs directory found.")
  quit(status = 0)
}
if (nzchar(Sys.readlink(runs_root))) stop("Refusing a symlinked tcga_runs root.", call. = FALSE)
runs_root <- normalizePath(runs_root)
protected <- character()
if (!is.null(protected_file)) {
  protected <- trimws(readLines(protected_file, warn = FALSE))
  protected <- unique(protected[nzchar(protected) & !startsWith(protected, "#")])
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", protected))) stop("Protected entries must be native run IDs, one per line.", call. = FALSE)
}
if (!is.null(verified_plan)) {
  protected_paths <- verified_plan$protected[grepl("^tcga_runs/[^/]+$", verified_plan$protected)]
  protected <- union(protected, sub("^tcga_runs/", "", protected_paths))
}
dirs <- list.dirs(runs_root, recursive = FALSE, full.names = FALSE)
if (!length(dirs)) {
  message("No run directories to archive.")
  quit(status = 0)
}
metadata <- lapply(dirs, function(id) {
  path <- file.path(runs_root, id)
  if (nzchar(Sys.readlink(path))) return(NULL)
  tryCatch(jsonlite::read_json(file.path(path, "run_metadata.json"), simplifyVector = FALSE), error = function(e) NULL)
})
names(metadata) <- dirs
# Expand known dependencies of protected pipeline/report runs. Missing legacy
# manifests cannot establish complete ancestry; callers must supply those IDs.
repeat {
  prior <- protected
  for (id in protected) {
    meta <- metadata[[id]]
    if (is.null(meta)) next
    protected <- union(protected, unlist(meta$dependency_run_ids %||% character(), use.names = FALSE))
    if (!is.null(meta$parent_run_id)) protected <- union(protected, meta$parent_run_id)
    inputs_path <- file.path(runs_root, id, "run_inputs.tsv")
    if (file.exists(inputs_path)) {
      inputs <- tryCatch(utils::read.delim(inputs_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (!is.null(inputs) && "path" %in% names(inputs)) {
        prefix <- paste0(runs_root, "/")
        values <- inputs$path[startsWith(inputs$path, prefix)]
        protected <- union(protected, sub("/.*$", "", substring(values, nchar(prefix) + 1L)))
      }
    }
  }
  if (setequal(prior, protected)) break
}
df <- data.frame(id = dirs, path = file.path(runs_root, dirs),
                 status = vapply(metadata, function(meta) meta$status %||% "unknown", character(1)),
                 mtime = file.info(file.path(runs_root, dirs))$mtime,
                 stringsAsFactors = FALSE)
eligible <- df$status %in% c("success", "failed") & !nzchar(Sys.readlink(df$path))
selected <- eligible
if (failed_only) selected <- selected & df$status == "failed"
if (!is.null(older_than_days)) selected <- selected & df$mtime < Sys.time() - older_than_days * 86400
if (!is.null(keep)) {
  recent <- which(eligible)[order(df$mtime[eligible], decreasing = TRUE)]
  kept <- head(recent, keep)
  selected[kept] <- FALSE
}
selected[is.na(selected)] <- FALSE
selected <- selected & !df$id %in% protected
if (!is.null(verified_plan)) selected <- selected & paste0("tcga_runs/", df$id) %in% verified_plan$candidates
targets <- df[selected, , drop = FALSE]
if (!nrow(targets)) {
  message("No unprotected completed run directories matched all criteria.")
  quit(status = 0)
}
cat(if (apply) "Archiving selected native runs:\n" else "[DRY RUN] Native runs that would be archived:\n")
for (id in targets$id) cat("  - ", id, "\n", sep = "")
if (!apply) {
  if (is.null(verified_plan)) cat("Lifecycle plan not supplied; this preview cannot establish archive eligibility.\n")
  quit(status = 0)
}
archive_root <- file.path(output_root(), "archive", "tcga_runs")
for (path in c(file.path(output_root(), "archive"), archive_root)) {
  link <- Sys.readlink(path)
  if (!is.na(link) && nzchar(link)) stop("Refusing a symlinked archive root.", call. = FALSE)
}
ensure_dir(archive_root)
batch <- tempfile(paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-"), tmpdir = archive_root)
if (!dir.create(batch)) stop("Cannot create archive batch.", call. = FALSE)
manifest <- data.frame(run_id = targets$id, source_path = targets$path,
                       archive_path = file.path(batch, targets$id), status = targets$status,
                       archive_status = "pending", archived_at = "", stringsAsFactors = FALSE)
write_manifest <- function() utils::write.table(manifest, file.path(batch, "ARCHIVE_MANIFEST.tsv"),
                                                sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(sort(protected), file.path(batch, "PROTECTED_RUNS.txt"))
if (!file.copy(archive_plan_file, file.path(batch, "archive_plan.json"), overwrite = FALSE)) stop("Cannot preserve archive plan.", call. = FALSE)
writeLines(c("Restore an archived run with file.rename(archive_path, source_path).",
             "First verify source_path does not exist. Never overwrite a current run.",
             "ARCHIVE_MANIFEST.tsv records each move; no run files were deleted."), file.path(batch, "RESTORE.txt"))
write_manifest()
for (i in seq_len(nrow(manifest))) {
  if (!file.rename(manifest$source_path[[i]], manifest$archive_path[[i]])) {
    stop("Archive move failed; earlier moves are recorded and recoverable in ", batch, call. = FALSE)
  }
  manifest$archive_status[[i]] <- "archived"
  manifest$archived_at[[i]] <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  write_manifest()
}
cat("Archive manifest: ", file.path(batch, "ARCHIVE_MANIFEST.tsv"), "\n", sep = "")
