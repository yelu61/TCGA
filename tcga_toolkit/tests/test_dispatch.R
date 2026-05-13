#!/usr/bin/env Rscript

# Smoke test: run_task.R must source successfully and the dispatch table must
# contain every advertised task. No real analysis runs.

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), mustWork = TRUE)
toolkit_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
scripts_dir <- file.path(toolkit_root, "scripts")

options(tcga_toolkit.root = toolkit_root)

`%||%` <- function(x, y) if (is.null(x)) y else x

passes <- 0L
failures <- 0L
log <- function(name, ok, detail = "") {
  if (ok) {
    passes <<- passes + 1L
    cat(sprintf("  PASS  %s\n", name))
  } else {
    failures <<- failures + 1L
    cat(sprintf("  FAIL  %s\n", name))
    if (nzchar(detail)) cat(sprintf("        %s\n", detail))
  }
}

# Load common.R and every task file (the actual run_task.R parses argv so we
# replicate just the file-loading steps).
source(file.path(scripts_dir, "common.R"))
task_files <- list.files(scripts_dir, pattern = "^task_.*\\.R$", full.names = TRUE)
for (tf in task_files) {
  err <- tryCatch({ source(tf); NULL }, error = function(e) conditionMessage(e))
  log(sprintf("source %s", basename(tf)), is.null(err), detail = err %||% "")
}

# Match advertised tasks against actual functions in the environment.
advertised_out <- suppressWarnings(system2(
  "Rscript",
  c(file.path(scripts_dir, "list_tasks.R"), "--json"),
  stdout = TRUE, stderr = TRUE
))
advertised <- tryCatch({
  parsed <- jsonlite::fromJSON(paste(advertised_out, collapse = "\n"), simplifyVector = FALSE)
  names(parsed)
}, error = function(e) character())

missing_fn <- c()
for (task in advertised) {
  fn <- sprintf("task_%s", task)
  if (!exists(fn, mode = "function")) {
    missing_fn <- c(missing_fn, task)
  }
}
log("every advertised task has a task_<name> function defined",
    length(missing_fn) == 0L,
    detail = if (length(missing_fn)) paste("missing:", paste(missing_fn, collapse = ", ")) else "")

cat(sprintf("\nDispatch: %s passed, %s failed.\n", passes, failures))
if (failures > 0L) quit(status = 1)
