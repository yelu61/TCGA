task_pan_cancer_expression <- function(config, ctx) {
  projects <- unlist(config$projects %||% available_projects(include_target = FALSE))
  values_list <- list()
  summary_rows <- list()

  for (project in projects) {
    info("Processing pan-cancer expression: %s", project)
    result <- tryCatch({
      data <- standardize_project_data(project)
      scores <- compute_feature_scores_from_config(data, config)
      merged <- merge(scores, data$clinical, by = "sample_id", all.x = TRUE)
      merged$project <- project
      merged$dataset_group <- ifelse(merged$Group %in% c("Normal"), "Normal", "Tumor")
      merged
    }, error = function(e) {
      warning(sprintf("Skipping %s: %s", project, e$message), call. = FALSE)
      NULL
    })
    if (is.null(result) || !nrow(result)) {
      next
    }
    values_list[[project]] <- result

    split_groups <- split(result, result$dataset_group)
    for (group_name in names(split_groups)) {
      df <- split_groups[[group_name]]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        project = project,
        feature_label = unique(df$feature_label),
        dataset_group = group_name,
        n = nrow(df),
        median = stats::median(df$score, na.rm = TRUE),
        mean = mean(df$score, na.rm = TRUE),
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    }

    tumor <- result[result$dataset_group == "Tumor", "score", drop = TRUE]
    normal <- result[result$dataset_group == "Normal", "score", drop = TRUE]
    if (length(tumor) >= 3 && length(normal) >= 3) {
      wt <- suppressWarnings(stats::wilcox.test(tumor, normal))
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        project = project,
        feature_label = unique(result$feature_label),
        dataset_group = "Tumor_vs_Normal",
        n = length(tumor) + length(normal),
        median = stats::median(tumor, na.rm = TRUE) - stats::median(normal, na.rm = TRUE),
        mean = mean(tumor, na.rm = TRUE) - mean(normal, na.rm = TRUE),
        p_value = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(values_list)) {
    fail("No pan-cancer expression results were generated.")
  }

  values_df <- do.call(rbind, values_list)
  summary_df <- do.call(rbind, summary_rows)
  feature_label <- unique(values_df$feature_label)[1]

  write_table_safe(values_df, file.path(ctx$results_dir, sprintf("%s_pan_cancer_expression_values.csv", feature_label)))
  write_table_safe(summary_df, file.path(ctx$results_dir, sprintf("%s_pan_cancer_expression_summary.csv", feature_label)))

  plot_df <- values_df[values_df$dataset_group %in% c("Tumor", "Normal"), ]
  if (nrow(plot_df) > 0 && length(unique(plot_df$dataset_group)) >= 2) {
    p_box <- plot_boxplot_by_group(plot_df, "project", "score", xlab = "Project", ylab = "Score", title = sprintf("%s pan-cancer expression", feature_label), fill_col = "dataset_group")
    if (!is.null(p_box)) {
      n_projects <- length(unique(plot_df$project))
      save_plot(p_box, sprintf("%s_pan_cancer_boxplot", feature_label), ctx$plots_dir, width = max(10, n_projects * 0.3), height = 6)
    }
  }

  write_report(
    ctx,
    sprintf("Pan-Cancer Expression for %s", feature_label),
    c(
      sprintf("- Projects evaluated: `%s`", length(unique(values_df$project))),
      sprintf("- Samples scored: `%s`", nrow(values_df)),
      sprintf("- Feature label: `%s`", feature_label),
      "",
      "## Outputs",
      sprintf("- `results/%s_pan_cancer_expression_values.csv`", feature_label),
      sprintf("- `results/%s_pan_cancer_expression_summary.csv`", feature_label)
    )
  )

  list(values = values_df, summary = summary_df)
}
