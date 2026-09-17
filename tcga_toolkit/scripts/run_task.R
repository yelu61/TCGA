#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
output_root <- NULL
overwrite <- FALSE
for (i in seq_along(args)) {
  if (args[[i]] == "--config" && i < length(args)) {
    config_path <- args[[i + 1L]]
  }
  if (args[[i]] == "--output-root" && i < length(args)) {
    output_root <- args[[i + 1L]]
  }
  if (args[[i]] == "--overwrite") {
    overwrite <- TRUE
  }
}
if (is.null(config_path)) {
  stop("Usage: Rscript tcga_toolkit/scripts/run_task.R --config <config.json> [--output-root <dir>]", call. = FALSE)
}
if (overwrite) stop("Native runs are immutable; --overwrite is no longer supported. Start a new run.", call. = FALSE)

script_flag <- grep("^--file=", commandArgs(), value = TRUE)[1]
if (is.na(script_flag)) {
  stop("Unable to resolve script path.", call. = FALSE)
}
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_flag), perl = TRUE)
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
options(tcga_toolkit.root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
if (!is.null(output_root)) {
  options(tcga_toolkit.output_root = normalizePath(output_root, mustWork = TRUE))
}
options(tcga_toolkit.overwrite = FALSE)

source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "run_provenance.R"))
source(file.path(script_dir, "task_audit_data.R"))
source(file.path(script_dir, "task_prepare_bulk_rna.R"))
source(file.path(script_dir, "task_run_deg.R"))
source(file.path(script_dir, "task_run_enrichment.R"))
source(file.path(script_dir, "task_run_clinical_assoc.R"))
source(file.path(script_dir, "task_run_tme.R"))
source(file.path(script_dir, "task_run_gtex_compare.R"))
source(file.path(script_dir, "task_validate_sc_signatures.R"))
source(file.path(script_dir, "task_pan_cancer_expression.R"))
source(file.path(script_dir, "task_stage_analysis.R"))
source(file.path(script_dir, "task_survival_map.R"))
source(file.path(script_dir, "task_subtype_analysis.R"))
source(file.path(script_dir, "task_gene_correlation_heatmap.R"))
source(file.path(script_dir, "task_gene_pair_coexpression.R"))
source(file.path(script_dir, "task_pipeline.R"))
source(file.path(script_dir, "task_render_report.R"))
source(file.path(script_dir, "task_maf_summary.R"))
source(file.path(script_dir, "task_mutation_survival.R"))
source(file.path(script_dir, "task_prognostic_model.R"))
source(file.path(script_dir, "task_ssgsea_score.R"))
source(file.path(script_dir, "task_tmb_analysis.R"))
source(file.path(script_dir, "task_cohort_qc.R"))
source(file.path(script_dir, "task_cnv_summary.R"))
source(file.path(script_dir, "task_methylation_diff.R"))
source(file.path(script_dir, "task_external_validate.R"))
source(file.path(script_dir, "task_immune_phenotype.R"))
source(file.path(script_dir, "task_wgcna_modules.R"))
source(file.path(script_dir, "task_drug_response.R"))

config <- read_config(config_path)
validate_config(config$task %||% "", config)
ctx <- new_provenance_context(config, config_path)

dispatch <- list(
  audit_data = task_audit_data,
  prepare_bulk_rna = task_prepare_bulk_rna,
  run_deg = task_run_deg,
  run_enrichment = task_run_enrichment,
  run_clinical_assoc = task_run_clinical_assoc,
  run_tme = task_run_tme,
  run_gtex_compare = task_run_gtex_compare,
  validate_sc_signatures = task_validate_sc_signatures,
  pan_cancer_expression = task_pan_cancer_expression,
  stage_analysis = task_stage_analysis,
  survival_map = task_survival_map,
  subtype_analysis = task_subtype_analysis,
  gene_correlation_heatmap = task_gene_correlation_heatmap,
  gene_pair_coexpression = task_gene_pair_coexpression,
  pipeline = task_pipeline,
  render_report = task_render_report,
  maf_summary = task_maf_summary,
  mutation_survival = task_mutation_survival,
  prognostic_model = task_prognostic_model,
  ssgsea_score = task_ssgsea_score,
  tmb_analysis = task_tmb_analysis,
  cohort_qc = task_cohort_qc,
  cnv_summary = task_cnv_summary,
  methylation_diff = task_methylation_diff,
  external_validate = task_external_validate,
  immune_phenotype = task_immune_phenotype,
  wgcna_modules = task_wgcna_modules,
  drug_response = task_drug_response
)

handler <- dispatch[[ctx$task]]
if (is.null(handler)) {
  fail("Unknown task: %s", ctx$task)
}

execution <- execute_provenance_task(config, ctx, handler)
if (execution$status == "failed") stop(execution$error, call. = FALSE)

cat(ctx$run_dir, "\n")
