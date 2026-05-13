task_stage_analysis <- function(config, ctx) {
  project <- config$project %||% fail("stage_analysis requires a project.")
  stage_column <- config$stage_column %||% "ajcc_pathologic_stage"
  data <- standardize_project_data(project)
  scores <- compute_feature_scores_from_config(data, config)
  merged <- merge(scores, data$clinical, by = "sample_id", all.x = TRUE)

  if (!(stage_column %in% colnames(merged))) {
    fail("Stage column %s not found in clinical table for %s", stage_column, project)
  }

  merged$stage_collapsed <- collapse_stage_vector(merged[[stage_column]])
  merged <- merged[!is.na(merged$stage_collapsed), , drop = FALSE]
  if (!nrow(merged)) {
    fail("No stage annotations available for %s after normalization.", project)
  }

  stage_summary <- do.call(rbind, lapply(split(merged, merged$stage_collapsed), function(df) {
    data.frame(
      stage = as.character(unique(df$stage_collapsed)),
      n = nrow(df),
      median = stats::median(df$score, na.rm = TRUE),
      mean = mean(df$score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  stage_summary$stage <- factor(stage_summary$stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"), ordered = TRUE)
  stage_summary <- stage_summary[order(stage_summary$stage), , drop = FALSE]

  kr <- suppressWarnings(stats::kruskal.test(score ~ stage_collapsed, data = merged))
  sp <- suppressWarnings(stats::cor.test(merged$score, as.numeric(merged$stage_collapsed), method = "spearman"))
  stats_df <- data.frame(
    project = project,
    feature_label = unique(merged$feature_label),
    n = nrow(merged),
    kruskal_p_value = kr$p.value,
    spearman_rho = unname(sp$estimate),
    spearman_p_value = sp$p.value,
    stringsAsFactors = FALSE
  )

  feature_label <- unique(merged$feature_label)[1]
  write_table_safe(merged, file.path(ctx$results_dir, sprintf("%s_%s_stage_values.csv", project, feature_label)))
  write_table_safe(stage_summary, file.path(ctx$results_dir, sprintf("%s_%s_stage_summary.csv", project, feature_label)))
  write_table_safe(stats_df, file.path(ctx$results_dir, sprintf("%s_%s_stage_stats.csv", project, feature_label)))

  p_box <- plot_boxplot_by_group(merged, "stage_collapsed", "score", title = sprintf("%s by stage in %s", feature_label, project))
  if (!is.null(p_box)) {
    save_plot(p_box, sprintf("%s_%s_stage_boxplot", project, feature_label), ctx$plots_dir)
  }

  write_report(
    ctx,
    sprintf("Stage Analysis for %s in %s", feature_label, project),
    c(
      sprintf("- Stage column: `%s`", stage_column),
      sprintf("- Samples with normalized stage: `%s`", nrow(merged)),
      sprintf("- Kruskal p-value: `%s`", signif(stats_df$kruskal_p_value, 4)),
      sprintf("- Spearman rho: `%s`", signif(stats_df$spearman_rho, 4)),
      "",
      "## Outputs",
      sprintf("- `results/%s_%s_stage_values.csv`", project, feature_label),
      sprintf("- `results/%s_%s_stage_summary.csv`", project, feature_label),
      sprintf("- `results/%s_%s_stage_stats.csv`", project, feature_label)
    )
  )

  list(values = merged, summary = stage_summary, stats = stats_df)
}
