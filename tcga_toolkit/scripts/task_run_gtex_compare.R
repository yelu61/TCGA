task_run_gtex_compare <- function(config, ctx) {
  project <- config$project %||% fail("run_gtex_compare requires a project.")
  gene <- config$gene %||% fail("run_gtex_compare requires a gene.")
  data <- standardize_project_data(project)
  tcga_values <- resolve_expression_vector(data, gene, use_tpm = TRUE, log2_transform = FALSE)

  gene_map_idx <- match(tcga_values$feature_id[[1L]], data$gene_map$feature_id)
  gene_id <- data$gene_map$gene_id[[gene_map_idx]]
  gtex_values <- extract_gtex_gene_vector(gene_id)
  if (is.null(gtex_values)) {
    fail("Gene %s (%s) was not found in GTEx matrix.", gene, gene_id)
  }

  gtex_samples <- read_gtex_samples()
  target_tissues <- resolve_gtex_tissues(project, explicit_tissues = unlist(config$gtex_tissues %||% character()))
  if (!length(target_tissues)) {
    fail("No GTEx tissues resolved for %s. Provide gtex_tissues in config.", project)
  }
  gtex_df <- data.frame(
    sample_id = names(gtex_values),
    gtex_tpm = as.numeric(gtex_values),
    stringsAsFactors = FALSE
  )
  gtex_df <- merge(gtex_df, gtex_samples[, c("sample_id", "tissue")], by = "sample_id", all.x = TRUE)
  gtex_df <- gtex_df[gtex_df$tissue %in% target_tissues, , drop = FALSE]

  tcga_df <- merge(tcga_values, data$clinical[, c("sample_id", "Group")], by = "sample_id", all.x = TRUE)
  tcga_df$dataset <- ifelse(tcga_df$Group %in% c("Normal"), "TCGA_Normal", "TCGA_Tumor")
  compare_df <- rbind(
    data.frame(dataset = tcga_df$dataset, tissue = project, value = tcga_df$value, sample_id = tcga_df$sample_id, stringsAsFactors = FALSE),
    data.frame(dataset = "GTEx", tissue = gtex_df$tissue, value = gtex_df$gtex_tpm, sample_id = gtex_df$sample_id, stringsAsFactors = FALSE)
  )

  summary_df <- do.call(rbind, lapply(split(compare_df, compare_df$dataset), function(df) {
    data.frame(
      dataset = unique(df$dataset),
      n = nrow(df),
      median = stats::median(df$value, na.rm = TRUE),
      mean = mean(df$value, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  write_table_safe(compare_df, file.path(ctx$results_dir, sprintf("%s_%s_gtex_compare_values.csv", project, gene)))
  write_table_safe(summary_df, file.path(ctx$results_dir, sprintf("%s_%s_gtex_compare_summary.csv", project, gene)))

  p_box <- plot_boxplot_by_group(compare_df, "dataset", "value", title = sprintf("%s expression in %s vs GTEx", gene, project), fill_col = "dataset")
  if (!is.null(p_box)) {
    n_groups <- length(unique(compare_df$dataset))
    save_plot(p_box, sprintf("%s_%s_gtex_boxplot", project, gene), ctx$plots_dir, width = max(8, n_groups * 1.2), height = 6)
  }

  write_report(
    ctx,
    sprintf("GTEx Compare for %s in %s", gene, project),
    c(
      sprintf("- GTEx tissues used: `%s`", paste(target_tissues, collapse = ", ")),
      sprintf("- TCGA samples summarized: `%s`", nrow(tcga_df)),
      sprintf("- GTEx samples summarized: `%s`", nrow(gtex_df)),
      "",
      "## Outputs",
      sprintf("- `results/%s_%s_gtex_compare_values.csv`", project, gene),
      sprintf("- `results/%s_%s_gtex_compare_summary.csv`", project, gene)
    )
  )

  list(values = compare_df, summary = summary_df)
}
