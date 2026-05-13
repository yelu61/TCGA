#!/usr/bin/env Rscript

# check_deps.R
# Report which R / Bioconductor packages are required by each task and which
# ones are missing on this machine. Emits install commands you can paste.
#
# Usage:
#   Rscript tcga_toolkit/scripts/check_deps.R
#   Rscript tcga_toolkit/scripts/check_deps.R --task run_deg
#   Rscript tcga_toolkit/scripts/check_deps.R --json

args <- commandArgs(trailingOnly = TRUE)
task_filter <- NULL
as_json <- FALSE
for (i in seq_along(args)) {
  if (args[[i]] == "--task" && i < length(args)) task_filter <- args[[i + 1L]]
  if (args[[i]] == "--json") as_json <- TRUE
}

# CRAN  -> install.packages
# Bioc  -> BiocManager::install
# GH    -> remotes::install_github(...) note string
dep_specs <- list(
  core = list(
    cran = c("jsonlite", "ggplot2"),
    bioc = c("SummarizedExperiment")
  ),
  audit_data = list(),
  prepare_bulk_rna = list(),
  run_deg = list(
    bioc = "DESeq2"
  ),
  run_enrichment = list(
    bioc = "clusterProfiler"
  ),
  run_clinical_assoc = list(
    cran = c("survival", "survminer")
  ),
  run_tme = list(
    gh = "IOBR/IOBR"
  ),
  run_gtex_compare = list(),
  validate_sc_signatures = list(
    cran = c("survival", "survminer")
  ),
  pan_cancer_expression = list(),
  stage_analysis = list(),
  survival_map = list(
    cran = c("survival", "survminer")
  ),
  subtype_analysis = list(),
  gene_correlation_heatmap = list(
    cran = c("readxl", "pheatmap")
  ),
  pipeline = list(),
  render_report = list(),
  maf_summary = list(
    bioc = "maftools"
  ),
  mutation_survival = list(
    cran = "survival",
    bioc = "maftools"
  ),
  prognostic_model = list(
    cran = c("glmnet", "survival", "survminer", "timeROC", "randomForestSRC")
  ),
  ssgsea_score = list(
    bioc = c("GSVA", "GSEABase")
  ),
  tmb_analysis = list(
    cran = "survival",
    bioc = "maftools"
  ),
  cohort_qc = list(
    cran = c("ggplot2", "matrixStats"),
    bioc = "sva"
  ),
  cnv_summary = list(
    cran = "ggplot2"
  ),
  methylation_diff = list(
    cran = "ggplot2"
  ),
  external_validate = list(
    cran = c("survival", "survminer")
  ),
  immune_phenotype = list(
    bioc = c("GSVA", "GSEABase"),
    gh = "IOBR/IOBR"
  ),
  wgcna_modules = list(
    cran = "WGCNA",
    bioc = c("impute", "preprocessCore")
  ),
  drug_response = list(
    cran = "oncoPredict",
    bioc = "sva"
  )
)

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

resolve_deps <- function(task) {
  spec <- dep_specs[[task]] %||% list()
  list(
    cran = unique(c(dep_specs$core$cran, spec$cran)),
    bioc = unique(c(dep_specs$core$bioc, spec$bioc)),
    gh = unique(spec$gh)
  )
}

tasks <- names(dep_specs)
tasks <- setdiff(tasks, "core")
if (!is.null(task_filter)) {
  if (!(task_filter %in% tasks)) {
    stop(sprintf("Unknown task: %s", task_filter), call. = FALSE)
  }
  tasks <- task_filter
}

per_task <- list()
missing_cran <- character()
missing_bioc <- character()
missing_gh <- character()

for (task in tasks) {
  deps <- resolve_deps(task)
  cran_status <- vapply(deps$cran, is_installed, logical(1))
  bioc_status <- vapply(deps$bioc, is_installed, logical(1))
  gh_pkgs <- vapply(deps$gh, function(spec) sub("^.+/", "", spec), character(1))
  gh_status <- vapply(gh_pkgs, is_installed, logical(1))

  missing_cran <- unique(c(missing_cran, deps$cran[!cran_status]))
  missing_bioc <- unique(c(missing_bioc, deps$bioc[!bioc_status]))
  missing_gh <- unique(c(missing_gh, deps$gh[!gh_status]))

  per_task[[task]] <- list(
    cran = lapply(seq_along(deps$cran), function(i) list(pkg = deps$cran[[i]], installed = unname(cran_status[[i]]))),
    bioc = lapply(seq_along(deps$bioc), function(i) list(pkg = deps$bioc[[i]], installed = unname(bioc_status[[i]]))),
    gh = lapply(seq_along(deps$gh), function(i) list(spec = deps$gh[[i]], installed = unname(gh_status[[i]])))
  )
}

install_cmds <- list()
if (length(missing_cran)) {
  install_cmds$cran <- sprintf('install.packages(c(%s))',
                               paste(sprintf('"%s"', missing_cran), collapse = ", "))
}
if (length(missing_bioc)) {
  install_cmds$bioc <- c(
    'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")',
    sprintf('BiocManager::install(c(%s))', paste(sprintf('"%s"', missing_bioc), collapse = ", "))
  )
}
if (length(missing_gh)) {
  install_cmds$gh <- c(
    'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")',
    sprintf('remotes::install_github(c(%s))', paste(sprintf('"%s"', missing_gh), collapse = ", "))
  )
}

result <- list(
  R_version = as.character(getRversion()),
  per_task = per_task,
  missing = list(cran = missing_cran, bioc = missing_bioc, gh = missing_gh),
  install = install_cmds
)

if (as_json) {
  cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")
} else {
  cat(sprintf("R version: %s\n\n", result$R_version))
  for (task in names(per_task)) {
    cat(sprintf("[%s]\n", task))
    spec <- per_task[[task]]
    for (entry in spec$cran) {
      cat(sprintf("  CRAN   %-30s %s\n", entry$pkg, ifelse(entry$installed, "ok", "MISSING")))
    }
    for (entry in spec$bioc) {
      cat(sprintf("  Bioc   %-30s %s\n", entry$pkg, ifelse(entry$installed, "ok", "MISSING")))
    }
    for (entry in spec$gh) {
      cat(sprintf("  GH     %-30s %s\n", entry$spec, ifelse(entry$installed, "ok", "MISSING")))
    }
    cat("\n")
  }
  if (length(install_cmds)) {
    cat("Install commands:\n")
    if (!is.null(install_cmds$cran)) cat("  ", install_cmds$cran, "\n")
    if (!is.null(install_cmds$bioc)) for (l in install_cmds$bioc) cat("  ", l, "\n")
    if (!is.null(install_cmds$gh)) for (l in install_cmds$gh) cat("  ", l, "\n")
  } else {
    cat("All dependencies satisfied.\n")
  }
}
