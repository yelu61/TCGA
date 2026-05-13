task_methylation_diff <- function(config, ctx) {
  project <- config$project %||% fail("methylation_diff requires a project.")
  meth_file <- config$meth_file
  if (is.null(meth_file) || !nzchar(meth_file)) {
    candidates <- c(
      file.path(project_root(), "0-Data", sprintf("%s_methylation.rda", project)),
      file.path(project_root(), "0-Data", sprintf("%s_methylation.tsv", project)),
      file.path(project_root(), "0-Data", sprintf("%s_methylation.csv", project)),
      file.path(project_root(), "0-Data", sprintf("%s_methyl.rda", project))
    )
    meth_file <- candidates[file.exists(candidates)][1]
  }

  if (is.na(meth_file) || is.null(meth_file) || !file.exists(meth_file)) {
    fail(paste0(
      "methylation_diff requires a methylation beta-value matrix.\n",
      "Place a CpG x sample matrix at one of:\n",
      "  0-Data/%s_methylation.rda  (matrix or data.frame, rows = CpG/IlmnID, cols = sample_id)\n",
      "  0-Data/%s_methylation.tsv / .csv\n",
      "Or pass an explicit `meth_file`.\n",
      "Download via TCGAbiolinks::GDCquery(project='%s', data.category='DNA Methylation', data.type='Methylation Beta Value', platform='Illumina Human Methylation 450')."
    ), project, project, project)
  }

  group_column <- config$group_column %||% "Group"
  numerator <- config$numerator %||% "Tumor"
  denominator <- config$denominator %||% "Normal"
  delta_beta_cutoff <- as.numeric(config$delta_beta_cutoff %||% 0.1)
  padj_cutoff <- as.numeric(config$padj_cutoff %||% 0.05)
  top_n <- as.integer(config$top_n %||% 50L)

  ext <- tolower(tools::file_ext(meth_file))
  if (ext %in% c("rda", "rdata")) {
    beta <- load_rdata_file(meth_file)
  } else if (ext == "tsv") {
    beta <- utils::read.delim(meth_file, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  } else if (ext == "csv") {
    beta <- utils::read.csv(meth_file, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  } else {
    fail("Unsupported methylation file format: %s", ext)
  }
  if (is.data.frame(beta)) beta <- as.matrix(beta)
  beta <- apply(beta, 2, as.numeric)
  if (max(beta, na.rm = TRUE) > 1.2) {
    info("Beta matrix appears to be on a 0-100 scale; rescaling to 0-1.")
    beta <- beta / 100
  }
  beta[beta < 0] <- 0
  beta[beta > 1] <- 1

  data <- standardize_project_data(project)
  clinical <- data$clinical
  if (!(group_column %in% colnames(clinical))) {
    fail("Group column %s not found in clinical table.", group_column)
  }

  ids_short <- substr(colnames(beta), 1, 15)
  clinical$sample_id_short <- substr(clinical$sample_id, 1, 15)
  meta <- clinical[match(ids_short, clinical$sample_id_short), , drop = FALSE]
  group <- as.character(meta[[group_column]])

  keep <- !is.na(group) & group %in% c(numerator, denominator) & !is.na(colSums(beta, na.rm = TRUE))
  if (sum(keep) < 6L) {
    fail("Need at least 3 samples per group for methylation_diff; got %s usable samples.", sum(keep))
  }
  beta <- beta[, keep, drop = FALSE]
  group <- group[keep]

  num_idx <- group == numerator
  den_idx <- group == denominator
  mean_num <- rowMeans(beta[, num_idx, drop = FALSE], na.rm = TRUE)
  mean_den <- rowMeans(beta[, den_idx, drop = FALSE], na.rm = TRUE)
  delta_beta <- mean_num - mean_den

  info("Running per-CpG Wilcoxon tests on %s CpGs (this can be slow on full 450k panels)…", nrow(beta))
  test_chunk <- function(rows) {
    vapply(rows, function(i) {
      x <- beta[i, num_idx]
      y <- beta[i, den_idx]
      if (length(stats::na.omit(x)) < 3 || length(stats::na.omit(y)) < 3) return(NA_real_)
      suppressWarnings(stats::wilcox.test(x, y)$p.value)
    }, numeric(1))
  }
  p_values <- test_chunk(seq_len(nrow(beta)))
  padj <- stats::p.adjust(p_values, method = "BH")

  res_df <- data.frame(
    cpg_id = rownames(beta),
    mean_numerator = mean_num,
    mean_denominator = mean_den,
    delta_beta = delta_beta,
    p_value = p_values,
    padj = padj,
    stringsAsFactors = FALSE
  )
  res_df <- res_df[order(res_df$padj, -abs(res_df$delta_beta), na.last = TRUE), , drop = FALSE]

  sig_df <- subset(res_df, !is.na(padj) & padj < padj_cutoff & abs(delta_beta) >= delta_beta_cutoff)
  hyper <- subset(sig_df, delta_beta > 0)
  hypo <- subset(sig_df, delta_beta < 0)

  write_table_safe(res_df, file.path(ctx$results_dir, sprintf("%s_methylation_diff.csv", project)))
  write_table_safe(utils::head(hyper, top_n),
                   file.path(ctx$results_dir, sprintf("%s_methylation_top_hyper.csv", project)))
  write_table_safe(utils::head(hypo, top_n),
                   file.path(ctx$results_dir, sprintf("%s_methylation_top_hypo.csv", project)))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_df <- res_df
    plot_df$category <- "NS"
    plot_df$category[plot_df$delta_beta >= delta_beta_cutoff & !is.na(plot_df$padj) & plot_df$padj < padj_cutoff] <- "Hyper"
    plot_df$category[plot_df$delta_beta <= -delta_beta_cutoff & !is.na(plot_df$padj) & plot_df$padj < padj_cutoff] <- "Hypo"
    plot_df$category <- factor(plot_df$category, levels = c("NS", "Hyper", "Hypo"))
    colours <- c("NS" = "grey70", "Hyper" = "#C44E52", "Hypo" = "#4C72B0")
    p_volcano <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_beta, y = -log10(padj), colour = category)) +
      ggplot2::geom_point(alpha = 0.5, size = 0.8) +
      ggplot2::scale_colour_manual(values = colours) +
      ggplot2::labs(title = sprintf("%s methylation %s vs %s", project, numerator, denominator),
                    x = "Δβ", y = "-log10 padj") +
      theme_tcga()
    save_plot(p_volcano, sprintf("%s_methylation_volcano", project), ctx$plots_dir, width = 9, height = 6)
  }

  write_report(
    ctx,
    sprintf("Methylation differential for %s (%s vs %s)", project, numerator, denominator),
    c(
      sprintf("- Methylation matrix: `%s`", meth_file),
      sprintf("- CpGs tested: `%s`", nrow(beta)),
      sprintf("- Samples used: `%s` (%s vs %s)", sum(keep), sum(num_idx), sum(den_idx)),
      sprintf("- Hyper-methylated (Δβ >= %s, padj < %s): `%s`", delta_beta_cutoff, padj_cutoff, nrow(hyper)),
      sprintf("- Hypo-methylated: `%s`", nrow(hypo)),
      "",
      "## Outputs",
      sprintf("- `results/%s_methylation_diff.csv`", project),
      sprintf("- `results/%s_methylation_top_hyper.csv`", project),
      sprintf("- `results/%s_methylation_top_hypo.csv`", project)
    )
  )

  list(results = res_df, hyper = hyper, hypo = hypo)
}
