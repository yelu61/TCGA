task_validate_sc_signatures <- function(config, ctx) {
  project <- config$project %||% fail("validate_sc_signatures requires a project.")
  signature_file <- config$signature_file %||% fail("validate_sc_signatures requires signature_file.")
  marker_table <- read_marker_table(signature_file)
  data <- standardize_project_data(project)

  scores <- compute_signature_score(data, marker_table, method = config$signature_method %||% "mean_z")
  merged <- merge(scores, data$clinical, by = "sample_id", all.x = TRUE)
  assoc <- clinical_associations(
    score_df = data.frame(sample_id = scores$sample_id, score = scores$score, stringsAsFactors = FALSE),
    clinical = data$clinical,
    columns = unlist(config$clinical_columns %||% c("Group", "gender", "race", "ajcc_pathologic_stage", "age_at_index"))
  )
  survival <- run_survival_models(
    score_df = data.frame(sample_id = scores$sample_id, score = scores$score, stringsAsFactors = FALSE),
    clinical = data$clinical
  )

  write_table_safe(merged, file.path(ctx$results_dir, sprintf("%s_signature_scores.csv", project)))
  write_table_safe(assoc, file.path(ctx$results_dir, sprintf("%s_signature_clinical_assoc.csv", project)))
  if (!is.null(survival)) {
    write_table_safe(survival$cox, file.path(ctx$results_dir, sprintf("%s_signature_survival_cox.csv", project)))
  }

  if ("Group" %in% colnames(merged) && length(unique(merged$Group)) >= 2) {
    p_box <- plot_boxplot_by_group(merged, "Group", "score", title = sprintf("Signature score by Group in %s", project))
    if (!is.null(p_box)) {
      save_plot(p_box, sprintf("%s_signature_boxplot", project), ctx$plots_dir)
    }
  }
  p_km <- plot_km_curve(
    score_df = data.frame(sample_id = scores$sample_id, score = scores$score, stringsAsFactors = FALSE),
    clinical = data$clinical,
    title = sprintf("Signature survival in %s", project)
  )
  if (!is.null(p_km)) {
    save_plot(p_km, sprintf("%s_signature_km", project), ctx$plots_dir, width = 7, height = 6)
  }

  write_report(
    ctx,
    sprintf("Single-cell Signature Validation for %s", project),
    c(
      sprintf("- Signature file: `%s`", signature_file),
      sprintf("- Signature genes matched: `%s`", unique(scores$n_genes)),
      sprintf("- Samples scored: `%s`", nrow(scores)),
      sprintf("- Clinical tests produced: `%s`", nrow(assoc)),
      sprintf("- Survival model available: `%s`", !is.null(survival)),
      "",
      "## Outputs",
      sprintf("- `results/%s_signature_scores.csv`", project),
      sprintf("- `results/%s_signature_clinical_assoc.csv`", project)
    )
  )

  list(scores = merged, associations = assoc, survival = survival)
}
