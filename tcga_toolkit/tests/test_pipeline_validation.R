#!/usr/bin/env Rscript
# Pipeline report paths are future outputs at preflight, with bounded references.
flag <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script <- gsub("~\\+~", " ", sub("^--file=", "", flag), perl = TRUE)
validator <- normalizePath(file.path(dirname(script), "../scripts/validate_config.R"))
tmp <- tempfile("pipeline-preflight-")
dir.create(tmp)
config_path <- file.path(tmp, "config.json")
check <- function(config, expected_ok, pattern) {
  jsonlite::write_json(config, config_path, auto_unbox = TRUE)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(validator, "--config", config_path, "--json")), stdout = TRUE, stderr = TRUE))
  result <- jsonlite::fromJSON(paste(output, collapse = "\n"))
  stopifnot(identical(result$ok, expected_ok), grepl(pattern, paste(output, collapse = "\n")))
}
pipeline <- list(task = "pipeline", steps = list(list(task = "audit_data"),
  list(task = "render_report", run_dirs = list("{{step0.run_dir}}"))))
check(pipeline, TRUE, "deferred")
for (bad in c("{{step1.run_dir}}", "{{step99.run_dir}}", "{{step0.typo}}")) {
  pipeline$steps[[2]]$run_dirs <- list(bad)
  check(pipeline, FALSE, "prior pipeline step")
}
check(list(task = "render_report", run_dirs = list("{{step0.run_dir}}")), FALSE, "prior pipeline step")
check(list(task = "render_report", run_dirs = list(file.path(tmp, "missing"))), FALSE, "entry not found")
unlink(tmp, recursive = TRUE)
cat("PASS: 6 pipeline report preflight checks; valid prior references deferred, future/unknown/standalone references rejected.\n")
