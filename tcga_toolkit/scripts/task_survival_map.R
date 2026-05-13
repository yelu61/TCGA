task_survival_map <- function(config, ctx) {
  projects <- unlist(config$projects %||% available_projects(include_target = FALSE))
  rows <- list()

  for (project in projects) {
    info("Processing survival map: %s", project)
    result <- tryCatch({
      data <- standardize_project_data(project)
      scores <- compute_feature_scores_from_config(data, config)
      sv <- run_survival_models(
        score_df = data.frame(sample_id = scores$sample_id, score = scores$score, stringsAsFactors = FALSE),
        clinical = data$clinical
      )
      if (is.null(sv)) {
        return(NULL)
      }
      data.frame(
        project = project,
        feature_label = unique(scores$feature_label),
        hazard_ratio = sv$cox$hazard_ratio[[1]],
        p_value = sv$cox$p_value[[1]],
        n = sv$km$n,
        median_cutoff = sv$km$median_cutoff,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning(sprintf("Skipping %s: %s", project, e$message), call. = FALSE)
      NULL
    })
    if (!is.null(result)) {
      rows[[length(rows) + 1L]] <- result
    }
  }

  if (!length(rows)) {
    fail("No survival map results were generated.")
  }

  map_df <- do.call(rbind, rows)
  map_df$direction <- ifelse(map_df$hazard_ratio >= 1, "Risk", "Protective")
  feature_label <- unique(map_df$feature_label)[1]

  write_table_safe(map_df, file.path(ctx$results_dir, sprintf("%s_survival_map.csv", feature_label)))

  p_forest <- plot_forest(map_df, title = sprintf("Survival Map for %s", feature_label))
  if (!is.null(p_forest)) {
    save_plot(p_forest, sprintf("%s_survival_forest", feature_label), ctx$plots_dir, width = 7, height = max(6, nrow(map_df) * 0.25))
  }

  write_report(
    ctx,
    sprintf("Pan-Cancer Survival Map for %s", feature_label),
    c(
      sprintf("- Projects with valid survival models: `%s`", nrow(map_df)),
      sprintf("- Risk associations (HR >= 1): `%s`", sum(map_df$direction == "Risk")),
      sprintf("- Protective associations (HR < 1): `%s`", sum(map_df$direction == "Protective")),
      "",
      "## Outputs",
      sprintf("- `results/%s_survival_map.csv`", feature_label)
    )
  )

  map_df
}
