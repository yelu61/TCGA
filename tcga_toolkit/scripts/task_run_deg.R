task_run_deg <- function(config, ctx) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    fail("run_deg requires the DESeq2 package.")
  }

  project <- config$project %||% fail("run_deg requires a project.")
  group_column <- config$group_column %||% "Group"
  numerator <- config$numerator %||% fail("run_deg requires numerator.")
  denominator <- config$denominator %||% fail("run_deg requires denominator.")
  min_count <- as.numeric(config$min_count %||% 10)
  min_fraction <- as.numeric(config$min_fraction %||% 0.25)

  data <- standardize_project_data(project)
  clinical <- select_group_samples(data, group_column, numerator, denominator)
  clinical <- clinical[match(colnames(data$counts), clinical$sample_id), , drop = FALSE]
  keep <- !is.na(clinical[[group_column]]) & as.character(clinical[[group_column]]) %in% c(numerator, denominator)
  clinical <- clinical[keep, , drop = FALSE]

  # Patient-level dedup: one sample per patient so replicate aliquots of the
  # same tumour are not treated as independent biological replicates.
  if (isTRUE(config$dedup_patient %||% TRUE) && "patient_short" %in% colnames(clinical)) {
    dup_patient <- duplicated(clinical$patient_short)
    if (any(dup_patient)) {
      info("run_deg: dropping %s duplicate-patient sample(s), keeping first per patient.", sum(dup_patient))
      clinical <- clinical[!dup_patient, , drop = FALSE]
    }
  }

  counts <- data$counts[, clinical$sample_id, drop = FALSE]
  # Previously non-integer input (e.g. TPM) was silently round()ed, which
  # fabricates count data; refuse instead.
  assert_integer_counts(counts, sprintf("run_deg counts matrix for %s", project))

  # Optional design covariates / explicit design formula.
  covariates <- unlist(config$covariates %||% character())
  missing_cov <- setdiff(covariates, colnames(clinical))
  if (length(missing_cov)) {
    fail("run_deg covariates not found in clinical table: %s", paste(missing_cov, collapse = ", "))
  }
  for (cv in covariates) {
    if (!is.numeric(clinical[[cv]])) clinical[[cv]] <- factor(as.character(clinical[[cv]]))
  }
  if (length(covariates)) {
    complete_cov <- stats::complete.cases(clinical[, covariates, drop = FALSE])
    if (any(!complete_cov)) {
      info("run_deg: dropping %s sample(s) with NA covariate values.", sum(!complete_cov))
      clinical <- clinical[complete_cov, , drop = FALSE]
    }
  }
  counts <- counts[, clinical$sample_id, drop = FALSE]

  gene_keep <- rowSums(counts >= min_count) >= ceiling(ncol(counts) * min_fraction)
  counts <- counts[gene_keep, , drop = FALSE]
  if (!nrow(counts)) {
    fail("No genes passed the count filter for %s", project)
  }

  clinical$condition <- factor(as.character(clinical[[group_column]]), levels = c(denominator, numerator))
  rownames(clinical) <- clinical$sample_id
  design_formula <- if (!is.null(config$design_formula)) {
    stats::as.formula(config$design_formula)
  } else if (length(covariates)) {
    stats::as.formula(paste("~", paste(covariates, collapse = " + "), "+ condition"))
  } else {
    ~ condition
  }
  missing_terms <- setdiff(all.vars(design_formula), colnames(clinical))
  if (length(missing_terms)) {
    fail("design_formula variables not found in clinical table: %s", paste(missing_terms, collapse = ", "))
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = clinical, design = design_formula)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", numerator, denominator)))
  res$feature_id <- rownames(res)
  res <- merge(res, data$gene_map[, c("feature_id", "gene_id", "gene_name", "gene_type")], by = "feature_id", all.x = TRUE)
  res <- res[order(res$padj, -abs(res$log2FoldChange), na.last = TRUE), ]

  out_csv <- file.path(ctx$results_dir, sprintf("%s_%s_vs_%s_deg.csv", project, numerator, denominator))
  write_table_safe(res, out_csv)
  saveRDS(res, file.path(ctx$objects_dir, sprintf("%s_%s_vs_%s_deg.rds", project, numerator, denominator)))

  sig <- subset(res, !is.na(padj) & padj < (config$padj_cutoff %||% 0.05))

  p_volcano <- plot_volcano(res, title = sprintf("%s: %s vs %s", project, numerator, denominator), padj_cutoff = config$padj_cutoff %||% 0.05, log2fc_cutoff = as.numeric(config$log2fc_cutoff %||% 1))
  if (!is.null(p_volcano)) {
    save_plot(p_volcano, sprintf("%s_%s_vs_%s_volcano", project, numerator, denominator), ctx$plots_dir)
  }

  write_report(
    ctx,
    sprintf("DEG: %s %s vs %s", project, numerator, denominator),
    c(
      sprintf("- Samples used: `%s`", nrow(clinical)),
      sprintf("- Design: `%s`", deparse(design_formula)),
      sprintf("- Patient-level dedup: `%s`", isTRUE(config$dedup_patient %||% TRUE)),
      sprintf("- Genes tested: `%s`", nrow(res)),
      sprintf("- Significant genes (padj < %s): `%s`", config$padj_cutoff %||% 0.05, nrow(sig)),
      "",
      "## Outputs",
      sprintf("- `results/%s_%s_vs_%s_deg.csv`", project, numerator, denominator),
      sprintf("- `objects/%s_%s_vs_%s_deg.rds`", project, numerator, denominator)
    )
  )

  res
}
