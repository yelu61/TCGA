task_run_clinical_assoc <- function(config, ctx) {
  project <- config$project %||% fail("run_clinical_assoc requires a project.")
  data <- standardize_project_data(project)

  score_df <- if (!is.null(config$gene)) {
    resolve_expression_vector(data, config$gene)
  } else if (!is.null(config$signature_file)) {
    markers <- read_marker_table(config$signature_file)
    scores <- compute_signature_score(data, markers, method = config$signature_method %||% "mean_z")
    data.frame(sample_id = scores$sample_id, value = scores$score, feature = basename(config$signature_file), stringsAsFactors = FALSE)
  } else {
    fail("run_clinical_assoc requires either gene or signature_file.")
  }

  merged <- merge(score_df, data$clinical, by = "sample_id", all.x = TRUE)
  assoc <- clinical_associations(
    score_df = data.frame(sample_id = merged$sample_id, score = merged$value, stringsAsFactors = FALSE),
    clinical = data.frame(sample_id = merged$sample_id, merged[, setdiff(colnames(merged), c("sample_id", "value", "feature", "feature_id")), drop = FALSE], stringsAsFactors = FALSE),
    columns = unlist(config$clinical_columns %||% c("Group", "gender", "race", "ajcc_pathologic_stage", "age_at_index"))
  )
  survival <- run_survival_models(
    score_df = data.frame(sample_id = merged$sample_id, score = merged$value, stringsAsFactors = FALSE),
    clinical = data.frame(sample_id = merged$sample_id, merged[, setdiff(colnames(merged), c("sample_id", "value", "feature", "feature_id")), drop = FALSE], stringsAsFactors = FALSE)
  )

  write_table_safe(merged, file.path(ctx$results_dir, sprintf("%s_feature_values.csv", project)))
  write_table_safe(assoc, file.path(ctx$results_dir, sprintf("%s_clinical_associations.csv", project)))
  if (!is.null(survival)) {
    write_table_safe(survival$cox, file.path(ctx$results_dir, sprintf("%s_survival_cox.csv", project)))
    write_json(
      list(
        median_cutoff = survival$km$median_cutoff,
        n = survival$km$n,
        groups = as.list(survival$km$groups)
      ),
      file.path(ctx$results_dir, sprintf("%s_survival_km.json", project))
    )
  }

  # Plotting
  group_col <- "Group"
  if (nrow(assoc) > 0 && any(assoc$test %in% c("wilcox", "kruskal"))) {
    best_cat <- assoc$variable[assoc$test %in% c("wilcox", "kruskal")][1]
    if (!is.na(best_cat) && best_cat %in% colnames(merged) && length(unique(merged[[best_cat]])) >= 2) {
      group_col <- best_cat
    }
  }
  if (group_col %in% colnames(merged) && length(unique(merged[[group_col]])) >= 2) {
    p_box <- plot_boxplot_by_group(merged, group_col, "value", title = sprintf("%s by %s", config$gene %||% basename(config$signature_file), group_col))
    if (!is.null(p_box)) {
      save_plot(p_box, sprintf("%s_clinical_boxplot", project), ctx$plots_dir)
    }
  }

  p_km <- plot_km_curve(
    score_df = data.frame(sample_id = merged$sample_id, score = merged$value, stringsAsFactors = FALSE),
    clinical = data$clinical,
    title = sprintf("%s survival", config$gene %||% basename(config$signature_file))
  )
  if (!is.null(p_km)) {
    save_plot(p_km, sprintf("%s_clinical_km", project), ctx$plots_dir, width = 7, height = 6)
  }

  write_report(
    ctx,
    sprintf("Clinical Associations for %s", project),
    c(
      sprintf("- Feature source: `%s`", config$gene %||% config$signature_file),
      sprintf("- Samples scored: `%s`", nrow(merged)),
      sprintf("- Clinical tests produced: `%s`", nrow(assoc)),
      sprintf("- Survival model available: `%s`", !is.null(survival)),
      "",
      "## Outputs",
      sprintf("- `results/%s_feature_values.csv`", project),
      sprintf("- `results/%s_clinical_associations.csv`", project)
    )
  )

  list(values = merged, associations = assoc, survival = survival)
}
