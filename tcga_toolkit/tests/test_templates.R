#!/usr/bin/env Rscript

# Smoke test: every JSON in tcga_toolkit/templates/ must pass validate_config
# (schema only — referenced data files may or may not exist locally).

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), mustWork = TRUE)
toolkit_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
scripts_dir <- file.path(toolkit_root, "scripts")
templates_dir <- file.path(toolkit_root, "templates")

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

templates <- list.files(templates_dir, pattern = "\\.json$", full.names = TRUE)
data_absence_pattern <- "(does not exist|not found|no expression rda|entry not found)"
for (tmpl in templates) {
  cmd <- paste("Rscript", shQuote(file.path(scripts_dir, "validate_config.R")), "--config", shQuote(tmpl))
  out <- suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = FALSE))
  status <- attr(out, "status") %||% 0L
  # A template passes if the schema itself is OK. Missing local data files
  # (CNV / methylation / external cohort / project rda) are tolerated — they
  # are user-provided, not part of the toolkit distribution.
  schema_problems <- out[grepl("^  - ", out)]
  schema_problems <- schema_problems[!grepl(data_absence_pattern, schema_problems)]
  ok <- status == 0 || length(schema_problems) == 0L
  log(basename(tmpl), ok, detail = paste(utils::tail(out, 4), collapse = " | "))
}

cat(sprintf("\nTemplate validation: %s passed, %s failed.\n", passes, failures))
if (failures > 0L) quit(status = 1)
