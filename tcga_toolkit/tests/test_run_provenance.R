#!/usr/bin/env Rscript
# Isolated native-run, pipeline and reversible cleanup fixtures. No cohort analysis.
flag <- grep("^--file=", commandArgs(), value = TRUE)[1L]
script <- gsub("~\\+~", " ", sub("^--file=", "", flag), perl = TRUE)
scripts <- normalizePath(file.path(dirname(script), "..", "scripts"))
options(tcga_toolkit.root = normalizePath(file.path(scripts, "..")), tcga_toolkit.overwrite = FALSE)
source(file.path(scripts, "common.R"))
source(file.path(scripts, "run_provenance.R"))
source(file.path(scripts, "task_pipeline.R"))
tmp <- tempfile("tcga-contract-")
dir.create(tmp)
options(tcga_toolkit.output_root = tmp)
config_file <- file.path(tmp, "config.json")
write_json(list(task = "audit_data"), config_file)
large_bytes <- provenance_rows()
large_bytes[1L, ] <- list("large.tsv", "table", "fixture", 1e12, "recorded")
write_provenance_rows(large_bytes, file.path(tmp, "large-bytes.tsv"))
stopifnot(grepl("\t1000000000000\t", readLines(file.path(tmp, "large-bytes.tsv"))[[2L]]))
input <- file.path(tmp, "external expression.tsv")
writeLines("gene\tsample\nG1\t5", input)
input_hash <- provenance_md5(input)
empty_dir <- file.path(tmp, "empty")
dir.create(empty_dir)
config <- list(task = "fixture", task_id = "success", expression_file = input,
               clinical_file = file.path(tmp, "missing.tsv"), signature_file = "https://example.invalid/genes.tsv",
               training_dir = empty_dir)
ctx <- new_provenance_context(config, config_file)
output <- execute_provenance_task(config, ctx, function(config, ctx) {
  writeLines("x,y\n1,2", file.path(ctx$results_dir, "table.csv"))
  writeLines("figure", file.path(ctx$plots_dir, "figure.svg"))
  list()
})
stopifnot(output$status == "success", output$metadata$manifest_schema_version == 1L,
          output$metadata$provenance_status == "complete")
inputs <- read.delim(file.path(ctx$run_dir, "run_inputs.tsv"), stringsAsFactors = FALSE)
artifacts <- read.delim(file.path(ctx$run_dir, "run_artifacts.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(names(inputs), c("path", "role", "md5", "bytes", "status")),
          setequal(inputs$status, c("recorded", "missing", "unresolved", "directory_empty")),
          inputs$md5[inputs$status == "recorded"] == input_hash,
          all(!startsWith(artifacts$path, "/")),
          setequal(artifacts$role, c("table", "figure", "provenance")),
          !any(artifacts$path %in% c("run_metadata.json", "run_inputs.tsv", "run_artifacts.tsv")))
for (field in c("config", "input_manifest", "artifact_manifest")) {
  stopifnot(output$metadata[[paste(field, "md5", sep = "_")]] ==
              provenance_md5(file.path(ctx$run_dir, output$metadata[[paste(field, "file", sep = "_")]])))
}
throws <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  !is.null(error) && grepl(pattern, error)
}
stopifnot(throws(execute_provenance_task(config, ctx, function(...) NULL), "overwrite"))
saved_id <- normalize_task_id
normalize_task_id <- function(...) "collision"
collision <- new_provenance_context(config, config_file)
writeLines("keep", file.path(collision$run_dir, "sentinel"))
stopifnot(throws(new_provenance_context(config, config_file), "already exists"),
          readLines(file.path(collision$run_dir, "sentinel")) == "keep")
normalize_task_id <- saved_id

config$task_id <- "handler-failure"
failed_ctx <- new_provenance_context(config, config_file)
failed <- execute_provenance_task(config, failed_ctx, function(config, ctx) {
  writeLines("partial", file.path(ctx$results_dir, "partial.txt"))
  stop("intentional handler failure")
})
stopifnot(failed$status == "failed", failed$metadata$provenance_status == "complete",
          grepl("intentional", failed$error), file.exists(file.path(failed_ctx$run_dir, "run_artifacts.tsv")))

config$task_id <- "changed-input"
changed_ctx <- new_provenance_context(config, config_file)
changed <- execute_provenance_task(config, changed_ctx, function(...) {
  writeLines("changed", input)
  list()
})
stopifnot(changed$status == "failed", changed$metadata$provenance_status == "failed",
          grepl("input changed", changed$error))
writeLines("gene\tsample\nG1\t5", input)

config$task_id <- "linked-artifact"
linked_ctx <- new_provenance_context(config, config_file)
linked <- execute_provenance_task(config, linked_ctx, function(config, ctx) {
  stopifnot(file.symlink(input, file.path(ctx$results_dir, "link.tsv")))
  list()
})
stopifnot(linked$status == "failed", grepl("symbolic links", linked$error))

dispatch <- list(fixture = function(config, ctx) {
  writeLines("gene\tvalue\nG1\t2", file.path(ctx$results_dir, "output.tsv"))
  list()
})
pipeline_config <- list(task = "pipeline", task_id = "pipeline", steps = list(
  list(task = "fixture"), list(task = "fixture", input_deg = "{{step0.results_dir}}/output.tsv")))
pipeline_ctx <- new_provenance_context(pipeline_config, config_file)
pipeline <- execute_provenance_task(pipeline_config, pipeline_ctx, task_pipeline)
stopifnot(pipeline$status == "success", length(pipeline$metadata$dependency_run_ids) == 2L)
for (child in pipeline$result$step_contexts) {
  meta <- jsonlite::read_json(file.path(child$run_dir, "run_metadata.json"))
  stopifnot(meta$manifest_schema_version == 1L, meta$parent_run_id == pipeline_ctx$task_id,
            meta$provenance_status == "complete", meta$status == "success")
}
cat("PASS native provenance: snapshots, digests, statuses, collisions, partial failure, mutation and pipeline\n")

run_cli <- function(script, args) {
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                   c(shQuote(file.path(scripts, script)), vapply(args, shQuote, character(1))),
                                   stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status") %||% 0L, output = paste(out, collapse = "\n"))
}
cli <- run_cli("run_task.R", c("--config", config_file, "--output-root", tmp))
stopifnot(cli$status == 0L)
audit <- list.dirs(file.path(tmp, "tcga_runs"), recursive = FALSE, full.names = TRUE)
audit <- audit[grepl("-audit_data$", audit)]
stopifnot(length(audit) == 1L, jsonlite::read_json(file.path(audit, "run_metadata.json"))$manifest_schema_version == 1L)
overwrite <- run_cli("run_task.R", c("--config", config_file, "--output-root", tmp, "--overwrite"))
stopifnot(overwrite$status != 0L, grepl("immutable", overwrite$output))
cat("PASS actual CLI fixture and --overwrite rejection\n")

cleanup_root <- file.path(tmp, "cleanup")
dir.create(file.path(cleanup_root, "tcga_runs"), recursive = TRUE)
make_run <- function(id, status, days = 40, extra = list()) {
  path <- file.path(cleanup_root, "tcga_runs", id)
  dir.create(path)
  writeLines(id, file.path(path, "result.txt"))
  write_json(c(list(task_id = id, status = status), extra), file.path(path, "run_metadata.json"))
  Sys.setFileTime(path, Sys.time() - days * 86400)
  path
}
old_failed <- make_run("old-failed", "failed")
old_success <- make_run("old-success", "success")
recent_failed <- make_run("recent-failed", "failed", 0)
protected_run <- make_run("protected-failed", "failed")
dependency <- make_run("dependency-failed", "failed")
report <- make_run("protected-report", "success", extra = list(dependency_run_ids = "dependency-failed"))
unknown <- make_run("unknown", "running")
unregistered <- make_run("unregistered-failed", "failed")
protected_file <- file.path(cleanup_root, "protected.txt")
writeLines(c("protected-failed", "protected-report"), protected_file)
base_args <- c("--output-root", cleanup_root, "--older-than", "30", "--failed-only", "--protected-runs", protected_file)
preview <- run_cli("cleanup_runs.R", base_args)
stopifnot(preview$status == 0L, grepl("DRY RUN", preview$output), grepl("old-failed", preview$output),
          !grepl("old-success|recent-failed|protected-failed|dependency-failed|unknown", preview$output),
          all(dir.exists(c(old_failed, old_success, recent_failed, protected_run, dependency, unknown))))
no_protection <- run_cli("cleanup_runs.R", c("--output-root", cleanup_root, "--failed-only", "--apply"))
stopifnot(no_protection$status != 0L, dir.exists(old_failed))
invalid <- run_cli("cleanup_runs.R", c("--output-root", cleanup_root, "--keep", "-1"))
stopifnot(invalid$status != 0L)
dir.create(file.path(cleanup_root, "provenance"))
dir.create(file.path(cleanup_root, "results"))
writeLines("fixture registry", file.path(cleanup_root, "provenance", "RUNS.tsv"))
writeLines("fixture artifacts", file.path(cleanup_root, "provenance", "ARTIFACTS.tsv"))
writeLines("manual note", file.path(cleanup_root, "results", "note.txt"))
fresh_files <- c("provenance/RUNS.tsv", "provenance/ARTIFACTS.tsv", "results/note.txt")
freshness <- list(sha256 = setNames(lapply(fresh_files, function(path) {
  digest::digest(file = file.path(cleanup_root, path), algo = "sha256", serialize = FALSE)
}), fresh_files), scope_files = list(results = list("results/note.txt")))
plan_file <- file.path(cleanup_root, "archive-plan.json")
write_json(list(project_root = normalizePath(cleanup_root),
                candidate_native_directories = list("tcga_runs/old-failed", "tcga_runs/dependency-failed"),
                protected_native_directories = list("tcga_runs/protected-report", "tcga_runs/protected-failed"),
                freshness = freshness), plan_file)
plan_args <- c(base_args, "--archive-plan", plan_file)
writeLines("changed manual note", file.path(cleanup_root, "results", "note.txt"))
stale <- run_cli("cleanup_runs.R", c(plan_args, "--apply"))
stopifnot(stale$status != 0L, grepl("Stale archive plan", stale$output), dir.exists(old_failed))
writeLines("manual note", file.path(cleanup_root, "results", "note.txt"))
writeLines("new manual work", file.path(cleanup_root, "results", "new.txt"))
stale <- run_cli("cleanup_runs.R", c(plan_args, "--apply"))
stopifnot(stale$status != 0L, grepl("Stale archive plan", stale$output), dir.exists(old_failed))
unlink(file.path(cleanup_root, "results", "new.txt"))
old_hash <- provenance_md5(file.path(old_failed, "result.txt"))
archived <- run_cli("cleanup_runs.R", c(plan_args, "--apply"))
stopifnot(archived$status == 0L, !dir.exists(old_failed),
          all(dir.exists(c(old_success, recent_failed, protected_run, dependency, unknown, unregistered))))
archive_manifest <- list.files(file.path(cleanup_root, "archive"), pattern = "ARCHIVE_MANIFEST.tsv", recursive = TRUE, full.names = TRUE)
stopifnot(length(archive_manifest) == 1L)
record <- read.delim(archive_manifest, stringsAsFactors = FALSE)
stopifnot(nrow(record) == 1L, record$run_id == "old-failed", record$archive_status == "archived",
          provenance_md5(file.path(record$archive_path, "result.txt")) == old_hash)
stopifnot(file.rename(record$archive_path, record$source_path), dir.exists(old_failed))
cat("PASS cleanup: dry-run default, AND filters, plan allowlist/freshness, protections, dependency closure, archive and restore\n")
cat("Fixture root: ", tmp, "\n", sep = "")
