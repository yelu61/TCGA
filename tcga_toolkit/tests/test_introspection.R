#!/usr/bin/env Rscript

# Smoke test: introspection scripts must run without errors and return the
# expected shape of output. None of these touch real cohort data.

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), mustWork = TRUE)
toolkit_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
scripts_dir <- file.path(toolkit_root, "scripts")

run_cmd <- function(args) {
  cmd <- paste("Rscript", paste(shQuote(args), collapse = " "))
  out <- suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = FALSE))
  status <- attr(out, "status") %||% 0L
  list(status = status, output = out)
}

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

# 1. list_tasks
r <- run_cmd(c(file.path(scripts_dir, "list_tasks.R"), "--json"))
ok <- r$status == 0 && any(grepl("\"audit_data\"", r$output, fixed = TRUE))
log("list_tasks.R --json contains audit_data", ok,
    detail = paste(utils::head(r$output, 3), collapse = " | "))

r <- run_cmd(c(file.path(scripts_dir, "list_tasks.R"), "--task", "prognostic_model"))
ok <- r$status == 0 && any(grepl("Lasso-Cox", r$output, fixed = TRUE))
log("list_tasks.R --task prognostic_model mentions Lasso-Cox", ok)

# 2. check_deps
r <- run_cmd(c(file.path(scripts_dir, "check_deps.R"), "--task", "cohort_qc"))
ok <- r$status == 0
log("check_deps.R --task cohort_qc runs", ok,
    detail = paste(utils::head(r$output, 3), collapse = " | "))

# 3. inspect_project --list
r <- run_cmd(c(file.path(scripts_dir, "inspect_project.R"), "--list"))
ok <- r$status == 0 && length(r$output) > 0
log("inspect_project.R --list runs", ok)

# 4. list_runs
r <- run_cmd(c(file.path(scripts_dir, "list_runs.R"), "--limit", "3"))
ok <- r$status == 0
log("list_runs.R --limit 3 runs", ok)

cat(sprintf("\nIntrospection: %s passed, %s failed.\n", passes, failures))
if (failures > 0L) quit(status = 1)
