#!/usr/bin/env Rscript

# list_tasks.R
# Print every task understood by run_task.R together with its required and
# optional config fields. Useful for Claude/agents to know what to put in a
# config file without reading scripts.
#
# Usage:
#   Rscript tcga_toolkit/scripts/list_tasks.R
#   Rscript tcga_toolkit/scripts/list_tasks.R --json
#   Rscript tcga_toolkit/scripts/list_tasks.R --task run_deg

args <- commandArgs(trailingOnly = TRUE)
as_json <- "--json" %in% args
task_filter <- NULL
for (i in seq_along(args)) {
  if (args[[i]] == "--task" && i < length(args)) {
    task_filter <- args[[i + 1L]]
  }
}

tasks <- list(
  audit_data = list(
    summary = "Inventory local TCGA/GTEx data and flag stale projects.",
    required = character(),
    optional = c("gdc_release_version", "gdc_release_date")
  ),
  prepare_bulk_rna = list(
    summary = "Build a project-level prepared expression bundle (counts + tpm + clinical).",
    required = "project",
    optional = c("counts_assay", "tpm_assay")
  ),
  run_deg = list(
    summary = "DESeq2 differential expression between two clinical groups.",
    required = c("project", "numerator", "denominator"),
    optional = c("group_column", "padj_cutoff", "log2fc_cutoff", "min_count", "min_fraction")
  ),
  run_enrichment = list(
    summary = "ORA/GSEA enrichment against local GMT collections.",
    required = "input_deg",
    optional = c("mode", "gmt_paths", "padj_cutoff", "log2fc_cutoff")
  ),
  run_clinical_assoc = list(
    summary = "Associate a gene or signature score with clinical variables and overall survival.",
    required = c("project", "gene_or_signature_file"),
    optional = c("clinical_columns", "signature_method")
  ),
  run_tme = list(
    summary = "Tumour microenvironment deconvolution (IOBR / ESTIMATE / CIBERSORT).",
    required = "project",
    optional = c("methods", "perm")
  ),
  run_gtex_compare = list(
    summary = "Compare a gene's expression between TCGA tumour/normal and GTEx tissue.",
    required = c("project", "gene"),
    optional = "gtex_tissues"
  ),
  validate_sc_signatures = list(
    summary = "Score a single-cell-derived gene signature on TCGA bulk and link to clinical/survival.",
    required = c("project", "signature_file"),
    optional = c("signature_method", "clinical_columns")
  ),
  pan_cancer_expression = list(
    summary = "Compare gene/signature expression across many TCGA projects.",
    required = "gene_or_signature_file",
    optional = "projects"
  ),
  stage_analysis = list(
    summary = "Test gene/signature against AJCC pathologic stage in a project.",
    required = c("project", "gene_or_signature_file"),
    optional = "stage_column"
  ),
  survival_map = list(
    summary = "Pan-cancer survival map: median-split KM and Cox per project.",
    required = "gene_or_signature_file",
    optional = "projects"
  ),
  subtype_analysis = list(
    summary = "Test gene/signature across molecular subtypes within a project.",
    required = c("project", "gene_or_signature_file"),
    optional = "subtype_column"
  ),
  gene_correlation_heatmap = list(
    summary = "Pearson/Spearman correlation between a target gene and a gene list across projects.",
    required = c("target_gene", "gene_list_file"),
    optional = c("sheet", "gene_column", "projects", "method", "sample_filter")
  ),
  pipeline = list(
    summary = "Run a sequence of tasks; later steps can reference earlier outputs.",
    required = "steps",
    optional = character()
  ),
  render_report = list(
    summary = "Render a unified self-contained HTML report from one or more run directories.",
    required = "run_dirs",
    optional = character()
  ),
  maf_summary = list(
    summary = "Mutation burden, oncoplot and top-mutated gene frequencies from project MAF.",
    required = "project",
    optional = "top_n"
  ),
  mutation_survival = list(
    summary = "Survival comparison between mutant and wildtype samples for a single gene.",
    required = c("project", "gene"),
    optional = character()
  ),
  prognostic_model = list(
    summary = "Train a multivariable prognostic model (Lasso-Cox or RSF) with train/test split and time-dependent ROC.",
    required = c("project", "feature_genes_or_signature"),
    optional = c("method", "test_fraction", "seed", "time_points", "feature_genes", "signature_file", "covariates")
  ),
  ssgsea_score = list(
    summary = "Single-sample pathway scoring (ssGSEA / GSVA) against any GMT collection; outputs sample x pathway matrix.",
    required = c("project", "gmt_paths"),
    optional = c("method", "min_size", "max_size", "clinical_columns")
  ),
  tmb_analysis = list(
    summary = "Tumour mutation burden per sample and TMB-stratified survival; uses local MAF or PanCancer_TMB.Rdata.",
    required = "project",
    optional = c("tmb_source", "tmb_threshold", "exome_size_mb")
  ),
  cohort_qc = list(
    summary = "Cohort QC: library size, gene detection, PCA outliers, optional ComBat-seq batch correction diagnostic.",
    required = "project",
    optional = c("batch_column", "outlier_sd", "min_library_size")
  ),
  cnv_summary = list(
    summary = "Copy-number alteration burden, top amplifications / deletions and CNV-burden survival.",
    required = "project",
    optional = c("cnv_file", "top_n", "amp_threshold", "del_threshold", "focal_threshold")
  ),
  methylation_diff = list(
    summary = "Differential methylation (Wilcoxon) between two clinical groups on a user-provided beta matrix.",
    required = "project",
    optional = c("meth_file", "group_column", "numerator", "denominator", "delta_beta_cutoff", "padj_cutoff", "top_n")
  ),
  external_validate = list(
    summary = "Validate a signature / risk score on an external cohort (CGGA, METABRIC, GEO, etc).",
    required = c("cohort_name", "expression_file", "clinical_file"),
    optional = c("feature_genes", "signature_file", "weight_file")
  ),
  immune_phenotype = list(
    summary = "Immune phenotype scoring (IOBR TME methods + immune subtype from PanCancer_subtypes).",
    required = "project",
    optional = c("tme_methods", "perm", "immune_subtype_file", "immune_gmt")
  ),
  wgcna_modules = list(
    summary = "WGCNA co-expression module detection plus module-trait correlation.",
    required = "project",
    optional = c("n_top_genes", "min_module_size", "merge_cut_height", "soft_power", "sample_filter", "trait_columns")
  ),
  drug_response = list(
    summary = "Predict drug IC50 per sample via oncoPredict (GDSC2 / GDSC1 / CTRP2) and test by group.",
    required = "project",
    optional = c("training_dir", "drug_set")
  )
)

if (!is.null(task_filter)) {
  if (!(task_filter %in% names(tasks))) {
    stop(sprintf("Unknown task: %s. Run without --task to see all options.", task_filter), call. = FALSE)
  }
  tasks <- tasks[task_filter]
}

if (as_json) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required for --json output.", call. = FALSE)
  }
  cat(jsonlite::toJSON(tasks, auto_unbox = TRUE, pretty = TRUE), "\n")
} else {
  for (name in names(tasks)) {
    spec <- tasks[[name]]
    cat(sprintf("- %s\n", name))
    cat(sprintf("    %s\n", spec$summary))
    if (length(spec$required)) {
      cat(sprintf("    required: %s\n", paste(spec$required, collapse = ", ")))
    } else {
      cat("    required: (none)\n")
    }
    if (length(spec$optional)) {
      cat(sprintf("    optional: %s\n", paste(spec$optional, collapse = ", ")))
    }
    cat("\n")
  }
}
