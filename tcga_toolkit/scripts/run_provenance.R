# Native TCGA run contract. Source after common.R. These helpers do not copy
# cohort data and do not claim to observe every file read by third-party code.

provenance_md5 <- function(path) {
  value <- unname(tools::md5sum(path))
  if (length(value) != 1L || is.na(value)) fail("Cannot hash provenance file: %s", path)
  value
}

provenance_rows <- function() {
  data.frame(path = character(), role = character(), md5 = character(),
             bytes = numeric(), status = character(), stringsAsFactors = FALSE)
}

write_provenance_rows <- function(rows, path) {
  if (any(grepl("[\t\r\n]", rows$path))) fail("Provenance paths cannot contain tabs or newlines.")
  # File sizes are contract integers, independent of R's scientific-notation options.
  rows$bytes <- ifelse(is.na(rows$bytes), "", sprintf("%.0f", rows$bytes))
  utils::write.table(rows, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

new_provenance_context <- function(config, config_path) {
  if (isTRUE(getOption("tcga_toolkit.overwrite", FALSE))) {
    fail("Native runs are immutable; --overwrite is no longer supported. Start a new run.")
  }
  task <- config$task %||% fail("Config must include a task field.")
  task_id <- normalize_task_id(task, config$task_id %||% config$project %||% config$label)
  root <- file.path(output_root(), "tcga_runs")
  ensure_dir(root)
  run_dir <- file.path(root, task_id)
  # Atomic directory creation also protects two invocations within one second.
  if (!dir.create(run_dir, recursive = FALSE, showWarnings = FALSE)) {
    fail("Run directory already exists or cannot be created: %s. Use a new task_id.", run_dir)
  }
  for (part in c("results", "plots", "objects")) ensure_dir(file.path(run_dir, part))
  list(task = task, task_id = task_id, run_dir = normalizePath(run_dir),
       results_dir = file.path(run_dir, "results"), plots_dir = file.path(run_dir, "plots"),
       objects_dir = file.path(run_dir, "objects"),
       config_path = normalizePath(config_path, mustWork = TRUE), toolkit_version = toolkit_version())
}

provenance_backend_identity <- function() {
  git <- Sys.which("git")
  if (!nzchar(git)) return(list(backend_revision = "unavailable", backend_dirty = NULL))
  git_value <- function(args) {
    value <- suppressWarnings(system2(git, c("-C", shQuote(project_root()), args),
                                     stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(value, "status")) && attr(value, "status") != 0L) return(NULL)
    value
  }
  revision <- git_value(c("rev-parse", "HEAD"))
  dirty <- git_value(c("status", "--porcelain", "--untracked-files=normal"))
  list(backend_revision = if (length(revision)) revision[[1L]] else "unavailable",
       backend_dirty = if (is.null(dirty)) NULL else length(dirty) > 0L)
}

collect_provenance_inputs <- function(config) {
  rows <- provenance_rows()
  visited_directories <- character()
  add <- function(path, role) {
    if (!is.character(path) || length(path) != 1L || !nzchar(path)) return(invisible(NULL))
    if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://|\\{\\{", path)) {
      rows[nrow(rows) + 1L, ] <<- list(path, role, "", NA_real_, "unresolved")
      return(invisible(NULL))
    }
    path <- path.expand(path)
    if (!grepl("^/", path)) path <- file.path(getwd(), path)
    path <- normalizePath(path, mustWork = FALSE)
    if (!file.exists(path)) {
      rows[nrow(rows) + 1L, ] <<- list(path, role, "", NA_real_, "missing")
    } else if (dir.exists(path)) {
      # Report inputs and other explicitly declared directories are expanded so
      # copied figures and tables retain their concrete dependency identities.
      key <- paste(role, path, sep = "\t")
      if (key %in% visited_directories) return(invisible(NULL))
      visited_directories <<- c(visited_directories, key)
      files <- list.files(path, recursive = FALSE, full.names = TRUE, all.files = TRUE,
                          no.. = TRUE, include.dirs = TRUE)
      if (!length(files)) rows[nrow(rows) + 1L, ] <<- list(path, role, "", NA_real_, "directory_empty")
      for (file in sort(files)) add(file, role)
    } else {
      rows[nrow(rows) + 1L, ] <<- list(normalizePath(path), role, provenance_md5(path),
                                      unname(file.info(path)$size), "recorded")
    }
  }
  walk <- function(obj, prefix = "") {
    if (!is.list(obj)) return(invisible(NULL))
    for (i in seq_along(obj)) {
      key <- if (!is.null(names(obj)) && nzchar(names(obj)[[i]])) names(obj)[[i]] else as.character(i)
      label <- if (nzchar(prefix)) paste(prefix, key, sep = ".") else key
      value <- obj[[i]]
      is_path <- grepl("(_files?|_paths?|_dirs?)$|^(input_deg|immune_gmt|provenance_inputs)$", key)
      is_output <- grepl("^(output|results|plots|objects|archive|cache)_", key)
      if (is_path && !is_output) {
        for (item in unlist(value, use.names = FALSE)) add(item, paste0("declared:", label))
      } else if (is.list(value)) walk(value, label)
    }
  }
  walk(config)
  # These are candidates, not a claim that every file was consumed. Tasks can
  # choose prepared data, cached objects or third-party references internally.
  projects <- unique(unlist(c(config$project, config$projects), use.names = FALSE))
  for (project in projects) {
    candidates <- project_files(project)
    for (kind in setdiff(names(candidates), "raw_dir")) {
      if (file.exists(candidates[[kind]])) add(candidates[[kind]], paste0("cohort_candidate:", kind))
    }
  }
  if (nrow(rows)) rows <- unique(rows[order(rows$path, rows$role), , drop = FALSE])
  rows
}

collect_provenance_artifacts <- function(run_dir) {
  root <- normalizePath(run_dir, mustWork = TRUE)
  walk <- function(directory) {
    children <- list.files(directory, recursive = FALSE, full.names = TRUE,
                            all.files = TRUE, no.. = TRUE, include.dirs = TRUE)
    if (any(nzchar(Sys.readlink(children)))) fail("Native run artifacts cannot be symbolic links.")
    directories <- children[file.info(children)$isdir]
    c(children, unlist(lapply(directories, walk), use.names = FALSE))
  }
  entries <- walk(root)
  paths <- entries[!file.info(entries)$isdir]
  rel <- substring(paths, nchar(root) + 2L)
  keep <- !rel %in% c("run_metadata.json", "run_inputs.tsv", "run_artifacts.tsv")
  paths <- paths[keep]
  rel <- rel[keep]
  rows <- provenance_rows()
  for (i in order(rel)) {
    path <- rel[[i]]
    role <- if (path %in% c("config_used.json", "sessionInfo.txt")) "provenance" else
      if (grepl("^plots/|\\.(pdf|png|svg|tiff?|jpe?g)$", path, ignore.case = TRUE)) "figure" else
      if (grepl("^report|\\.(html|qmd|md)$", path, ignore.case = TRUE)) "report" else
      if (grepl("^objects/|\\.(rds|rda|rdata)$", path, ignore.case = TRUE)) "object" else
      if (grepl("\\.(csv|tsv|xlsx?)$", path, ignore.case = TRUE)) "table" else "other"
    rows[nrow(rows) + 1L, ] <- list(path, role, provenance_md5(paths[[i]]),
                                   unname(file.info(paths[[i]])$size), "recorded")
  }
  rows
}

execute_provenance_task <- function(config, ctx, handler) {
  metadata <- c(list(manifest_schema_version = 1L, backend = "tcga-toolkit",
                     task = ctx$task, task_id = ctx$task_id, toolkit_version = ctx$toolkit_version,
                     config_path = ctx$config_path, run_dir = ctx$run_dir,
                     started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                     status = "running", provenance_status = "pending",
                     config_file = "config_used.json", environment_file = "sessionInfo.txt",
                     input_manifest_file = "run_inputs.tsv", artifact_manifest_file = "run_artifacts.tsv",
                     provenance_scope = "declared-config-paths-and-known-local-cohort-candidates; not a complete runtime IO audit"),
                provenance_backend_identity())
  if (!is.null(ctx$parent_run_id)) metadata$parent_run_id <- ctx$parent_run_id
  metadata_path <- file.path(ctx$run_dir, "run_metadata.json")
  if (file.exists(metadata_path) || file.exists(file.path(ctx$run_dir, "config_used.json"))) {
    fail("Refusing to overwrite an executed native run: %s", ctx$run_dir)
  }
  result <- NULL
  error_message <- NULL
  inputs <- NULL
  prepared <- FALSE
  tryCatch({
    write_json(config, file.path(ctx$run_dir, metadata$config_file))
    metadata$config_md5 <- provenance_md5(file.path(ctx$run_dir, metadata$config_file))
    inputs <- collect_provenance_inputs(config)
    write_provenance_rows(inputs, file.path(ctx$run_dir, metadata$input_manifest_file))
    metadata$input_manifest_md5 <- provenance_md5(file.path(ctx$run_dir, metadata$input_manifest_file))
    write_json(metadata, metadata_path)
    prepared <- TRUE
    result <- handler(config, ctx)
  }, error = function(e) error_message <<- conditionMessage(e))
  provenance_error <- NULL
  tryCatch({
    writeLines(capture.output(utils::sessionInfo()), file.path(ctx$run_dir, metadata$environment_file))
    if (!prepared) fail("Input provenance preparation did not complete.")
    # A completed task must still refer to the exact config snapshot it used.
    if (!identical(metadata$config_md5, provenance_md5(file.path(ctx$run_dir, metadata$config_file)))) {
      fail("Config snapshot was changed during execution.")
    }
    if (!identical(metadata$input_manifest_md5, provenance_md5(file.path(ctx$run_dir, metadata$input_manifest_file)))) {
      fail("Input manifest was changed during execution.")
    }
    declared <- inputs[inputs$status == "recorded" & grepl("^declared:", inputs$role), , drop = FALSE]
    for (i in seq_len(nrow(declared))) {
      if (!file.exists(declared$path[[i]]) || !identical(provenance_md5(declared$path[[i]]), declared$md5[[i]])) {
        fail("Declared input changed during execution: %s", declared$path[[i]])
      }
    }
    if (is.list(result) && !is.null(result$step_contexts)) {
      metadata$dependency_run_ids <- unname(as.list(vapply(result$step_contexts, `[[`, character(1), "task_id")))
    }
    artifacts <- collect_provenance_artifacts(ctx$run_dir)
    write_provenance_rows(artifacts, file.path(ctx$run_dir, metadata$artifact_manifest_file))
    metadata$artifact_manifest_md5 <- provenance_md5(file.path(ctx$run_dir, metadata$artifact_manifest_file))
  }, error = function(e) provenance_error <<- conditionMessage(e))
  errors <- unique(c(error_message, provenance_error))
  metadata$status <- if (length(errors)) "failed" else "success"
  metadata$provenance_status <- if (is.null(provenance_error)) "complete" else "failed"
  metadata$error <- if (length(errors)) paste(errors, collapse = "; ") else NULL
  metadata$completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  write_json(metadata, metadata_path)
  list(result = result, status = metadata$status, error = metadata$error, metadata = metadata)
}
