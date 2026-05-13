#!/usr/bin/env Rscript

# validate_config.R
# Dry-run validation of a task config without executing analysis.
# Verifies: task is known, required fields present, referenced files exist,
# project data exists for project-bound tasks, pipeline steps are well-formed.
#
# Usage:
#   Rscript tcga_toolkit/scripts/validate_config.R --config <path.json>
#   Rscript tcga_toolkit/scripts/validate_config.R --config <path.json> --json

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
as_json <- FALSE
for (i in seq_along(args)) {
  if (args[[i]] == "--config" && i < length(args)) {
    config_path <- args[[i + 1L]]
  }
  if (args[[i]] == "--json") {
    as_json <- TRUE
  }
}
if (is.null(config_path)) {
  stop("Usage: Rscript tcga_toolkit/scripts/validate_config.R --config <path.json> [--json]", call. = FALSE)
}

script_flag <- grep("^--file=", commandArgs(), value = TRUE)[1]
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_flag), perl = TRUE)
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
options(tcga_toolkit.root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
source(file.path(script_dir, "common.R"))

known_tasks <- c(
  "audit_data", "prepare_bulk_rna", "run_deg", "run_enrichment",
  "run_clinical_assoc", "run_tme", "run_gtex_compare",
  "validate_sc_signatures", "pan_cancer_expression", "stage_analysis",
  "survival_map", "subtype_analysis", "gene_correlation_heatmap",
  "pipeline", "render_report", "maf_summary", "mutation_survival",
  "prognostic_model", "ssgsea_score", "tmb_analysis", "cohort_qc",
  "cnv_summary", "methylation_diff", "external_validate",
  "immune_phenotype", "wgcna_modules", "drug_response"
)

issues <- list()
notes <- list()
add_issue <- function(msg) issues[[length(issues) + 1L]] <<- msg
add_note <- function(msg) notes[[length(notes) + 1L]] <<- msg

check_file <- function(path, label) {
  if (is.null(path) || (is.character(path) && !any(nzchar(path)))) {
    return(invisible(NULL))
  }
  # Skip pipeline variable placeholders (resolved at runtime)
  if (grepl("\\{\\{", path)) {
    return(invisible(NULL))
  }
  # Resolve relative to project root when cwd differs from repo root
  root <- project_root()
  resolved <- path
  if (!is.null(root) && !file.exists(path)) {
    alt <- file.path(root, path)
    if (file.exists(alt)) resolved <- alt
  }
  if (!file.exists(resolved)) {
    add_issue(sprintf("%s does not exist: %s", label, path))
  } else {
    add_note(sprintf("%s present: %s", label, path))
  }
}

check_project_data <- function(project) {
  if (is.null(project)) return(invisible(NULL))
  if (is.list(project)) project <- unlist(project)
  if (!is.character(project) || !any(nzchar(project))) return(invisible(NULL))
  for (p in project[nzchar(project)]) {
    files <- project_files(p)
    if (!file.exists(files$mrna) && !file.exists(files$deg_input)) {
      add_issue(sprintf("Project %s: no expression rda found at %s or %s", p, files$mrna, files$deg_input))
    } else {
      add_note(sprintf("Project %s: expression rda present", p))
    }
  }
}

validate_one <- function(cfg, prefix = "") {
  task <- cfg$task
  if (is.null(task) || !nzchar(task)) {
    add_issue(sprintf("%stask field is missing", prefix))
    return(invisible(NULL))
  }
  if (!(task %in% known_tasks)) {
    add_issue(sprintf("%sunknown task: %s (run list_tasks.R to see valid tasks)", prefix, task))
    return(invisible(NULL))
  }

  tryCatch(validate_config(task, cfg), error = function(e) {
    add_issue(sprintf("%s%s", prefix, conditionMessage(e)))
  })

  if (task == "pipeline") {
    steps <- cfg$steps %||% list()
    if (!length(steps)) {
      add_issue(sprintf("%spipeline.steps is empty", prefix))
    }
    for (i in seq_along(steps)) {
      validate_one(steps[[i]], prefix = sprintf("%sstep[%s].", prefix, i - 1))
    }
    return(invisible(NULL))
  }

  proj <- cfg[["project"]]
  if (!is.null(proj) && is.character(proj) && any(nzchar(proj))) {
    check_project_data(proj)
  }
  if (!is.null(cfg$projects)) {
    check_project_data(cfg$projects)
  }
  check_file(cfg$signature_file, "signature_file")
  check_file(cfg$input_deg, "input_deg")
  check_file(cfg$gene_list_file, "gene_list_file")
  if (!is.null(cfg$gmt_paths)) {
    for (p in unlist(cfg$gmt_paths)) check_file(p, "gmt_paths entry")
  }
  if (!is.null(cfg$run_dirs)) {
    for (p in unlist(cfg$run_dirs)) {
      if (!dir.exists(p)) {
        add_issue(sprintf("%srun_dirs entry not found: %s", prefix, p))
      } else {
        add_note(sprintf("%srun_dir present: %s", prefix, p))
      }
    }
  }
}

raw <- tryCatch(read_config(config_path), error = function(e) {
  add_issue(sprintf("Failed to parse config: %s", conditionMessage(e)))
  NULL
})

if (!is.null(raw)) {
  validate_one(raw)
}

ok <- !length(issues)
result <- list(
  config_path = normalizePath(config_path, mustWork = FALSE),
  ok = ok,
  issues = unlist(issues),
  notes = unlist(notes)
)

if (as_json) {
  cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")
} else {
  cat(sprintf("Config: %s\n", result$config_path))
  cat(sprintf("Status: %s\n", ifelse(ok, "OK", "FAILED")))
  if (length(result$issues)) {
    cat("\nIssues:\n")
    for (msg in result$issues) cat(sprintf("  - %s\n", msg))
  }
  if (length(result$notes)) {
    cat("\nNotes:\n")
    for (msg in result$notes) cat(sprintf("  - %s\n", msg))
  }
}

if (!ok) {
  quit(status = 1)
}
