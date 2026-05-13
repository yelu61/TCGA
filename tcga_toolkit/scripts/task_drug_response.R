task_drug_response <- function(config, ctx) {
  if (!requireNamespace("oncoPredict", quietly = TRUE)) {
    fail(paste0("drug_response requires the oncoPredict CRAN package.\n",
                "install.packages('oncoPredict')"))
  }
  project <- config$project %||% fail("drug_response requires a project.")
  training_dir <- config$training_dir
  if (is.null(training_dir) || !nzchar(training_dir)) {
    candidates <- c(
      file.path(project_root(), "1-Input", "oncoPredict_training_data"),
      file.path(project_root(), "1-Input", "DataFiles", "Training Data"),
      file.path(project_root(), "0-Data", "oncoPredict_training_data")
    )
    training_dir <- candidates[dir.exists(candidates)][1]
  }
  if (is.null(training_dir) || is.na(training_dir) || !dir.exists(training_dir)) {
    fail(paste0(
      "drug_response needs GDSC training data for oncoPredict.\n",
      "Download from https://osf.io/c6tfx/ (GDSC1 / GDSC2 / CTRP2) and unpack to\n",
      "  1-Input/oncoPredict_training_data/\n",
      "Or pass training_dir in config. Each training set is a folder with\n",
      "GDSC2_Expr.rds and GDSC2_Res.rds (or CTRP2_Expr / CTRP2_Res)."
    ))
  }
  drug_set <- config$drug_set %||% "GDSC2"
  expr_path <- file.path(training_dir, sprintf("%s_Expr.rds", drug_set))
  res_path <- file.path(training_dir, sprintf("%s_Res.rds", drug_set))
  if (!file.exists(expr_path) || !file.exists(res_path)) {
    fail("Training files not found:\n  %s\n  %s", expr_path, res_path)
  }

  data <- standardize_project_data(project)
  clinical <- data$clinical
  expr <- 2^(log2(data$tpm + 1)) - 1
  rownames(expr) <- data$gene_map$gene_name %||% rownames(expr)
  expr <- expr[!is.na(rownames(expr)) & nzchar(rownames(expr)), , drop = FALSE]
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  expr <- as.matrix(expr)

  train_expr <- readRDS(expr_path)
  train_res <- readRDS(res_path)

  info("Running oncoPredict::calcPhenotype for %s drugs in %s on %s samples (this can take 10-60 min)…",
       ncol(train_res), drug_set, ncol(expr))

  prev_files <- list.files(ctx$results_dir, full.names = TRUE)
  oncoPredict::calcPhenotype(
    trainingExprData = train_expr,
    trainingPtype = train_res,
    testExprData = expr,
    batchCorrect = "eb",
    powerTransformPhenotype = TRUE,
    removeLowVaryingGenes = 0.2,
    minNumSamples = 10,
    printOutput = FALSE,
    folder = ctx$results_dir
  )
  new_files <- setdiff(list.files(ctx$results_dir, full.names = TRUE), prev_files)
  pred_file <- new_files[grepl("DrugPredictions", basename(new_files), ignore.case = TRUE)][1]
  if (is.na(pred_file)) {
    pred_file <- new_files[grepl("\\.csv$", new_files)][1]
  }
  if (is.na(pred_file) || !file.exists(pred_file)) {
    fail("oncoPredict ran but did not produce a predictions CSV.")
  }

  predictions <- utils::read.csv(pred_file, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  predictions$sample_id <- rownames(predictions)
  predictions <- predictions[, c("sample_id", setdiff(colnames(predictions), "sample_id"))]
  out_predictions <- file.path(ctx$results_dir, sprintf("%s_%s_drug_predictions.csv", project, drug_set))
  write_table_safe(predictions, out_predictions)

  drug_cols <- setdiff(colnames(predictions), "sample_id")
  summary_df <- data.frame(
    drug = drug_cols,
    median_ic50 = vapply(predictions[, drug_cols, drop = FALSE], stats::median, FUN.VALUE = 0, na.rm = TRUE),
    mean_ic50 = vapply(predictions[, drug_cols, drop = FALSE], mean, FUN.VALUE = 0, na.rm = TRUE),
    sd_ic50 = vapply(predictions[, drug_cols, drop = FALSE], stats::sd, FUN.VALUE = 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  group_diff <- NULL
  if ("Group" %in% colnames(clinical)) {
    merged <- merge(predictions, clinical[, c("sample_id", "Group")], by = "sample_id")
    groups <- unique(merged$Group)
    if (length(groups) >= 2L) {
      rows <- list()
      for (d in drug_cols) {
        v <- merged[[d]]
        if (all(is.na(v))) next
        test <- tryCatch(suppressWarnings(stats::kruskal.test(v ~ as.factor(merged$Group))),
                         error = function(e) NULL)
        if (is.null(test)) next
        rows[[length(rows) + 1L]] <- data.frame(
          drug = d,
          test = "kruskal",
          p_value = test$p.value,
          n = sum(!is.na(v)),
          stringsAsFactors = FALSE
        )
      }
      if (length(rows)) {
        group_diff <- do.call(rbind, rows)
        group_diff$padj <- stats::p.adjust(group_diff$p_value, method = "BH")
        group_diff <- group_diff[order(group_diff$padj), , drop = FALSE]
        write_table_safe(group_diff, file.path(ctx$results_dir, sprintf("%s_%s_drug_group_diff.csv", project, drug_set)))
      }
    }
  }

  write_table_safe(summary_df, file.path(ctx$results_dir, sprintf("%s_%s_drug_summary.csv", project, drug_set)))

  if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(group_diff)) {
    top_drugs <- utils::head(group_diff$drug, 12)
    plot_df <- merge(predictions, clinical[, c("sample_id", "Group")], by = "sample_id")
    long <- do.call(rbind, lapply(top_drugs, function(d) {
      data.frame(drug = d, Group = plot_df$Group, IC50 = plot_df[[d]], stringsAsFactors = FALSE)
    }))
    long <- long[!is.na(long$IC50), , drop = FALSE]
    if (nrow(long)) {
      p <- ggplot2::ggplot(long, ggplot2::aes(x = Group, y = IC50, fill = Group)) +
        ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
        ggplot2::facet_wrap(~ drug, scales = "free_y") +
        ggplot2::labs(title = sprintf("%s top differential drugs (%s)", project, drug_set), y = "Predicted IC50") +
        theme_tcga()
      save_plot(p, sprintf("%s_%s_drug_top_boxplots", project, drug_set), ctx$plots_dir, width = 12, height = 8)
    }
  }

  write_report(
    ctx,
    sprintf("Drug response prediction for %s (%s)", project, drug_set),
    c(
      sprintf("- Training set: `%s`", drug_set),
      sprintf("- Samples scored: `%s`", nrow(predictions)),
      sprintf("- Drugs scored: `%s`", length(drug_cols)),
      sprintf("- Group differential tests: `%s`", if (!is.null(group_diff)) nrow(group_diff) else 0),
      "",
      "## Outputs",
      sprintf("- `results/%s_%s_drug_predictions.csv`", project, drug_set),
      sprintf("- `results/%s_%s_drug_summary.csv`", project, drug_set)
    )
  )

  list(predictions = predictions, summary = summary_df, group_diff = group_diff,
       drug_set = drug_set)
}
