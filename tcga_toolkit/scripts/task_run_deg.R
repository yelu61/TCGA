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
  counts <- round(data$counts[, clinical$sample_id, drop = FALSE])

  gene_keep <- rowSums(counts >= min_count) >= ceiling(ncol(counts) * min_fraction)
  counts <- counts[gene_keep, , drop = FALSE]
  if (!nrow(counts)) {
    fail("No genes passed the count filter for %s", project)
  }

  clinical$condition <- factor(as.character(clinical[[group_column]]), levels = c(denominator, numerator))
  rownames(clinical) <- clinical$sample_id
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = clinical, design = ~ condition)
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
