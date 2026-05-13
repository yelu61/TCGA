task_subtype_analysis <- function(config, ctx) {
  project <- config$project %||% fail("subtype_analysis requires a project.")
  subtype_column <- config$subtype_column %||% "Subtype_Selected"
  data <- standardize_project_data(project)
  scores <- compute_feature_scores_from_config(data, config)
  clinical_sub <- join_project_subtypes(data, subtype_column = subtype_column)
  merged <- merge(scores, clinical_sub, by = "sample_id", all.x = TRUE)
  merged <- merged[!is.na(merged$subtype) & nzchar(as.character(merged$subtype)), , drop = FALSE]
  if (!nrow(merged)) {
    fail("No subtype-matched samples available for %s", project)
  }

  subtype_summary <- do.call(rbind, lapply(split(merged, merged$subtype), function(df) {
    data.frame(
      subtype = unique(df$subtype),
      n = nrow(df),
      median = stats::median(df$score, na.rm = TRUE),
      mean = mean(df$score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  subtype_summary <- subtype_summary[order(-subtype_summary$n, subtype_summary$subtype), , drop = FALSE]

  stats_df <- data.frame(
    project = project,
    feature_label = unique(merged$feature_label),
    n = nrow(merged),
    subtype_count = nrow(subtype_summary),
    stringsAsFactors = FALSE
  )
  if (nrow(subtype_summary) >= 2) {
    kr <- suppressWarnings(stats::kruskal.test(score ~ subtype, data = merged))
    stats_df$kruskal_p_value <- kr$p.value
  } else {
    stats_df$kruskal_p_value <- NA_real_
  }

  feature_label <- unique(merged$feature_label)[1]
  write_table_safe(merged, file.path(ctx$results_dir, sprintf("%s_%s_subtype_values.csv", project, feature_label)))
  write_table_safe(subtype_summary, file.path(ctx$results_dir, sprintf("%s_%s_subtype_summary.csv", project, feature_label)))
  write_table_safe(stats_df, file.path(ctx$results_dir, sprintf("%s_%s_subtype_stats.csv", project, feature_label)))

  p_box <- plot_boxplot_by_group(merged, "subtype", "score", title = sprintf("%s by subtype in %s", feature_label, project))
  if (!is.null(p_box)) {
    save_plot(p_box, sprintf("%s_%s_subtype_boxplot", project, feature_label), ctx$plots_dir)
  }

  write_report(
    ctx,
    sprintf("Subtype Analysis for %s in %s", feature_label, project),
    c(
      sprintf("- Subtype column: `%s`", subtype_column),
      sprintf("- Samples with subtype labels: `%s`", nrow(merged)),
      sprintf("- Unique subtypes: `%s`", nrow(subtype_summary)),
      sprintf("- Kruskal p-value: `%s`", ifelse(is.na(stats_df$kruskal_p_value), "NA", signif(stats_df$kruskal_p_value, 4))),
      "",
      "## Outputs",
      sprintf("- `results/%s_%s_subtype_values.csv`", project, feature_label),
      sprintf("- `results/%s_%s_subtype_summary.csv`", project, feature_label),
      sprintf("- `results/%s_%s_subtype_stats.csv`", project, feature_label)
    )
  )

  list(values = merged, summary = subtype_summary, stats = stats_df)
}
