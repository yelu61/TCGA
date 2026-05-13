task_mutation_survival <- function(config, ctx) {
  project <- config$project %||% fail("mutation_survival requires a project.")
  gene <- config$gene %||% fail("mutation_survival requires a gene.")

  files <- project_files(project)
  if (!file.exists(files$maf)) {
    fail("Missing MAF file for %s: %s", project, files$maf)
  }

  maf_df <- load_rdata_file(files$maf)
  if (!is.data.frame(maf_df)) {
    fail("MAF object for %s is not a data.frame.", project)
  }
  if (!("Hugo_Symbol" %in% colnames(maf_df)) || !("Tumor_Sample_Barcode" %in% colnames(maf_df))) {
    fail("MAF data for %s is missing Hugo_Symbol or Tumor_Sample_Barcode columns.", project)
  }

  data <- standardize_project_data(project)
  clinical <- data$clinical

  # Normalize barcodes: keep first 12 or 16 characters depending on format
  maf_barcodes <- as.character(maf_df$Tumor_Sample_Barcode)
  clinical_barcodes <- as.character(clinical$sample_id)

  # Try matching by truncating MAF barcodes to clinical barcode length
  maf_short <- unique(substr(maf_barcodes, 1, nchar(clinical_barcodes[1])))
  gene_maf <- maf_df[maf_df$Hugo_Symbol == gene, ]
  mutated_short <- unique(substr(as.character(gene_maf$Tumor_Sample_Barcode), 1, nchar(clinical_barcodes[1])))

  clinical$mutation_status <- ifelse(substr(clinical_barcodes, 1, nchar(clinical_barcodes[1])) %in% mutated_short, "Mutant", "WT")
  clinical$mutation_status <- factor(clinical$mutation_status, levels = c("WT", "Mutant"))

  n_mutant <- sum(clinical$mutation_status == "Mutant")
  n_wt <- sum(clinical$mutation_status == "WT")

  write_table_safe(clinical[, c("sample_id", "mutation_status")], file.path(ctx$results_dir, sprintf("%s_%s_mutation_status.csv", project, gene)))

  # Survival analysis using existing helper
  score_df <- data.frame(
    sample_id = clinical$sample_id,
    score = as.numeric(clinical$mutation_status == "Mutant"),
    stringsAsFactors = FALSE
  )
  survival <- run_survival_models(score_df, clinical)

  if (!is.null(survival)) {
    write_table_safe(survival$cox, file.path(ctx$results_dir, sprintf("%s_%s_mutation_survival_cox.csv", project, gene)))
    write_json(
      list(
        median_cutoff = survival$km$median_cutoff,
        n = survival$km$n,
        groups = as.list(survival$km$groups)
      ),
      file.path(ctx$results_dir, sprintf("%s_%s_mutation_survival_km.json", project, gene))
    )

    p_km <- plot_km_curve(score_df, clinical, title = sprintf("%s %s mutation survival", project, gene))
    if (!is.null(p_km)) {
      save_plot(p_km, sprintf("%s_%s_mutation_km", project, gene), ctx$plots_dir, width = 7, height = 6)
    }
  }

  p_box <- plot_boxplot_by_group(clinical, "mutation_status", "score", title = sprintf("%s %s mutation status", project, gene))
  # score is 0/1, boxplot not very informative; skip if not meaningful

  write_report(
    ctx,
    sprintf("Mutation Survival for %s in %s", gene, project),
    c(
      sprintf("- Mutant samples: `%s`", n_mutant),
      sprintf("- WT samples: `%s`", n_wt),
      sprintf("- Survival model available: `%s`", !is.null(survival)),
      "",
      "## Outputs",
      sprintf("- `results/%s_%s_mutation_status.csv`", project, gene),
      sprintf("- `results/%s_%s_mutation_survival_cox.csv`", project, gene)
    )
  )

  list(clinical = clinical, survival = survival)
}
