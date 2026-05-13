task_prepare_bulk_rna <- function(config, ctx) {
  project <- config$project %||% fail("prepare_bulk_rna requires a project.")
  data <- standardize_project_data(project)

  utils::write.csv(data$counts, file.path(ctx$objects_dir, sprintf("%s_counts.csv", project)), quote = FALSE)
  utils::write.csv(data$tpm, file.path(ctx$objects_dir, sprintf("%s_tpm.csv", project)), quote = FALSE)
  write_table_safe(data$clinical, file.path(ctx$results_dir, sprintf("%s_clinical.csv", project)))
  write_table_safe(data$gene_map, file.path(ctx$results_dir, sprintf("%s_gene_map.csv", project)))
  saveRDS(data, file.path(ctx$objects_dir, sprintf("%s_prepared.rds", project)))

  group_counts <- sort(table(data$clinical$Group), decreasing = TRUE)
  group_text <- if (length(group_counts)) paste(names(group_counts), group_counts, collapse = ", ") else "None"

  write_report(
    ctx,
    sprintf("Prepared Bulk RNA for %s", project),
    c(
      sprintf("- Samples: `%s`", ncol(data$counts)),
      sprintf("- Features after deduplication: `%s`", nrow(data$counts)),
      sprintf("- Groups detected: `%s`", group_text),
      "",
      "## Outputs",
      sprintf("- `objects/%s_prepared.rds`", project),
      sprintf("- `objects/%s_counts.csv`", project),
      sprintf("- `objects/%s_tpm.csv`", project),
      sprintf("- `results/%s_clinical.csv`", project),
      sprintf("- `results/%s_gene_map.csv`", project)
    )
  )

  data
}
