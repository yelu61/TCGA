task_cohort_qc <- function(config, ctx) {
  project <- config$project %||% fail("cohort_qc requires a project.")
  outlier_sd <- as.numeric(config$outlier_sd %||% 3)
  min_library_size <- as.numeric(config$min_library_size %||% 1e6)
  batch_column <- config$batch_column

  data <- standardize_project_data(project)
  counts <- round(data$counts)
  clinical <- data$clinical

  library_size <- colSums(counts, na.rm = TRUE)
  detection_rate <- colMeans(counts >= 1, na.rm = TRUE)

  per_sample <- data.frame(
    sample_id = colnames(counts),
    library_size = library_size,
    log10_library_size = log10(library_size + 1),
    detected_genes = colSums(counts >= 1, na.rm = TRUE),
    detection_rate = detection_rate,
    stringsAsFactors = FALSE
  )
  per_sample <- merge(per_sample,
                      clinical[, c("sample_id", "Group", "patient_short"), drop = FALSE],
                      by = "sample_id", all.x = TRUE)
  per_sample$low_library_flag <- per_sample$library_size < min_library_size

  expr_log <- log2(data$tpm + 1)
  variances <- apply(expr_log, 1, stats::var, na.rm = TRUE)
  top_idx <- order(variances, decreasing = TRUE)[seq_len(min(2000L, length(variances)))]
  expr_top <- expr_log[top_idx, , drop = FALSE]
  expr_top <- expr_top[stats::complete.cases(expr_top), , drop = FALSE]

  pca_df <- NULL
  outlier_ids <- character()
  if (nrow(expr_top) >= 100 && ncol(expr_top) >= 5) {
    pca <- tryCatch(stats::prcomp(t(expr_top), center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (!is.null(pca)) {
      pcs <- pca$x[, seq_len(min(4L, ncol(pca$x))), drop = FALSE]
      pca_df <- data.frame(sample_id = rownames(pcs), pcs, stringsAsFactors = FALSE)
      pca_df <- merge(pca_df, per_sample[, c("sample_id", "Group", "library_size", "low_library_flag")], by = "sample_id", all.x = TRUE)

      pc_z <- scale(pcs)
      flag <- apply(abs(pc_z), 1, max, na.rm = TRUE) >= outlier_sd
      pca_df$pca_outlier <- flag
      outlier_ids <- pca_df$sample_id[flag]

      var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
      ve_df <- data.frame(
        PC = sprintf("PC%s", seq_along(var_explained)),
        var_explained = var_explained,
        cumulative = cumsum(var_explained),
        stringsAsFactors = FALSE
      )
      ve_df <- ve_df[seq_len(min(10L, nrow(ve_df))), , drop = FALSE]
      write_table_safe(ve_df, file.path(ctx$results_dir, sprintf("%s_qc_pca_variance.csv", project)))
    }
  }
  per_sample$pca_outlier <- per_sample$sample_id %in% outlier_ids

  batch_summary <- NULL
  combat_status <- "not_run"
  if (!is.null(batch_column) && batch_column %in% colnames(clinical)) {
    batch_summary <- as.data.frame(table(clinical[[batch_column]], useNA = "ifany"), stringsAsFactors = FALSE)
    names(batch_summary) <- c("batch", "n_samples")
    write_table_safe(batch_summary, file.path(ctx$results_dir, sprintf("%s_qc_batch_counts.csv", project)))

    if (requireNamespace("sva", quietly = TRUE) && isTRUE(config$run_combat_seq %||% FALSE)) {
      batches <- clinical[[batch_column]][match(colnames(counts), clinical$sample_id)]
      group_col <- if ("Group" %in% colnames(clinical)) clinical$Group[match(colnames(counts), clinical$sample_id)] else NULL
      combat_status <- tryCatch({
        keep <- !is.na(batches) & nzchar(as.character(batches))
        if (sum(keep) < 10 || length(unique(batches[keep])) < 2) {
          "skipped_insufficient_batches"
        } else {
          adj_counts <- sva::ComBat_seq(
            counts = counts[, keep, drop = FALSE],
            batch = as.character(batches[keep]),
            group = if (!is.null(group_col)) as.character(group_col[keep]) else NULL
          )
          saveRDS(adj_counts, file.path(ctx$objects_dir, sprintf("%s_qc_combatseq_adjusted_counts.rds", project)))
          "completed"
        }
      }, error = function(e) sprintf("error: %s", conditionMessage(e)))
    }
  }

  write_table_safe(per_sample, file.path(ctx$results_dir, sprintf("%s_qc_per_sample.csv", project)))
  if (!is.null(pca_df)) {
    write_table_safe(pca_df, file.path(ctx$results_dir, sprintf("%s_qc_pca.csv", project)))
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_lib <- ggplot2::ggplot(per_sample, ggplot2::aes(x = log10_library_size)) +
      ggplot2::geom_histogram(bins = 40, fill = "#55A868", colour = "white") +
      ggplot2::geom_vline(xintercept = log10(min_library_size + 1), linetype = "dashed", colour = "red") +
      ggplot2::labs(title = sprintf("%s library size", project),
                    x = "log10(library size)", y = "Samples") +
      theme_tcga()
    save_plot(p_lib, sprintf("%s_qc_library_size", project), ctx$plots_dir, width = 8, height = 5)

    if (!is.null(pca_df)) {
      p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, colour = Group, shape = pca_outlier)) +
        ggplot2::geom_point(size = 2.2, alpha = 0.85) +
        ggplot2::labs(title = sprintf("%s PCA on top-variance genes", project)) +
        theme_tcga()
      save_plot(p_pca, sprintf("%s_qc_pca", project), ctx$plots_dir, width = 8, height = 6)
    }
  }

  summary_lines <- c(
    sprintf("- Samples: `%s`", ncol(counts)),
    sprintf("- Genes: `%s`", nrow(counts)),
    sprintf("- Low-library (<%s reads): `%s`", min_library_size, sum(per_sample$low_library_flag, na.rm = TRUE)),
    sprintf("- PCA outliers (>%s SD): `%s`", outlier_sd, length(outlier_ids)),
    sprintf("- Outlier sample IDs: `%s`", paste(utils::head(outlier_ids, 10), collapse = ", ")),
    sprintf("- Batch column: `%s`", batch_column %||% "(not provided)"),
    sprintf("- ComBat-seq status: `%s`", combat_status),
    "",
    "## Outputs",
    sprintf("- `results/%s_qc_per_sample.csv`", project),
    sprintf("- `results/%s_qc_pca.csv` (if PCA succeeded)", project)
  )

  write_report(ctx, sprintf("Cohort QC for %s", project), summary_lines)

  list(per_sample = per_sample, pca = pca_df, batch = batch_summary,
       combat_status = combat_status, outliers = outlier_ids)
}
