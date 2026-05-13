task_maf_summary <- function(config, ctx) {
  if (!requireNamespace("maftools", quietly = TRUE)) {
    fail("task_maf_summary requires the maftools package.")
  }

  project <- config$project %||% fail("maf_summary requires a project.")
  files <- project_files(project)
  if (!file.exists(files$maf)) {
    fail("Missing MAF file for %s: %s", project, files$maf)
  }

  maf_df <- load_rdata_file(files$maf)
  if (!is.data.frame(maf_df)) {
    fail("MAF object for %s is not a data.frame.", project)
  }

  required_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")
  missing_cols <- setdiff(required_cols, colnames(maf_df))
  if (length(missing_cols)) {
    fail("MAF data for %s is missing required columns: %s", project, paste(missing_cols, collapse = ", "))
  }

  maf_obj <- maftools::read.maf(maf = as.data.frame(maf_df), isTCGA = TRUE)

  # Summary outputs
  summary_dt <- maftools::getSampleSummary(maf_obj)
  summary_df <- as.data.frame(summary_dt)
  write_table_safe(summary_df, file.path(ctx$results_dir, sprintf("%s_maf_summary.csv", project)))

  gene_summary_dt <- maftools::getGeneSummary(maf_obj)
  gene_summary_df <- as.data.frame(gene_summary_dt)
  write_table_safe(gene_summary_df, file.path(ctx$results_dir, sprintf("%s_gene_mutation_frequency.csv", project)))

  # Plots
  pdf_path <- file.path(ctx$plots_dir, sprintf("%s_maf_summary.pdf", project))
  png_path <- file.path(ctx$plots_dir, sprintf("%s_maf_summary.png", project))
  grDevices::pdf(pdf_path, width = 10, height = 8)
  maftools::plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
  grDevices::dev.off()
  grDevices::png(png_path, width = 1600, height = 1200, res = 150)
  maftools::plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
  grDevices::dev.off()

  top_n <- as.integer(config$top_n %||% 20)
  pdf_path2 <- file.path(ctx$plots_dir, sprintf("%s_oncoplot_top%s.pdf", project, top_n))
  png_path2 <- file.path(ctx$plots_dir, sprintf("%s_oncoplot_top%s.png", project, top_n))
  grDevices::pdf(pdf_path2, width = 10, height = max(6, top_n * 0.3))
  maftools::oncoplot(maf = maf_obj, top = top_n)
  grDevices::dev.off()
  grDevices::png(png_path2, width = 1600, height = max(900, top_n * 40), res = 150)
  maftools::oncoplot(maf = maf_obj, top = top_n)
  grDevices::dev.off()

  write_report(
    ctx,
    sprintf("MAF Summary for %s", project),
    c(
      sprintf("- Samples in MAF: `%s`", nrow(summary_df)),
      sprintf("- Unique mutated genes: `%s`", nrow(gene_summary_df)),
      sprintf("- Top N oncoplot: `%s`", top_n),
      "",
      "## Outputs",
      sprintf("- `results/%s_maf_summary.csv`", project),
      sprintf("- `results/%s_gene_mutation_frequency.csv`", project),
      sprintf("- `plots/%s_maf_summary.pdf/png`", project),
      sprintf("- `plots/%s_oncoplot_top%s.pdf/png`", project, top_n)
    )
  )

  list(summary = summary_df, gene_summary = gene_summary_df, maf_obj = maf_obj)
}
