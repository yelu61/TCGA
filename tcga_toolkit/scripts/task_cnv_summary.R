task_cnv_summary <- function(config, ctx) {
  project <- config$project %||% fail("cnv_summary requires a project.")
  top_n <- as.integer(config$top_n %||% 25L)
  amp_threshold <- as.numeric(config$amp_threshold %||% 0.3)
  del_threshold <- as.numeric(config$del_threshold %||% -0.3)
  focal_threshold <- as.integer(config$focal_threshold %||% 1L)

  cnv_file <- config$cnv_file
  if (is.null(cnv_file) || !nzchar(cnv_file)) {
    candidates <- c(
      file.path(project_root(), "0-Data", sprintf("%s_cnv.rda", project)),
      file.path(project_root(), "0-Data", sprintf("%s_gistic.rda", project)),
      file.path(project_root(), "0-Data", sprintf("%s_cnv.tsv", project)),
      file.path(project_root(), "0-Data", sprintf("%s_cnv.csv", project))
    )
    cnv_file <- candidates[file.exists(candidates)][1]
  }

  if (is.na(cnv_file) || is.null(cnv_file) || !file.exists(cnv_file)) {
    fail(paste0(
      "cnv_summary requires a CNV matrix.\n",
      "Place a gene-level CNV file at one of:\n",
      "  0-Data/%s_cnv.rda  (a matrix or data.frame with genes x samples)\n",
      "  0-Data/%s_gistic.rda\n",
      "  0-Data/%s_cnv.tsv / .csv\n",
      "Or pass an explicit `cnv_file` in the config.\n",
      "Download from GDC via TCGAbiolinks::GDCquery(project='%s', data.category='Copy Number Variation', data.type='Gene Level Copy Number')."
    ), project, project, project, project)
  }

  cnv <- NULL
  ext <- tolower(tools::file_ext(cnv_file))
  if (ext == "rda" || ext == "rdata") {
    cnv <- load_rdata_file(cnv_file)
  } else if (ext == "tsv") {
    cnv <- utils::read.delim(cnv_file, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  } else if (ext == "csv") {
    cnv <- utils::read.csv(cnv_file, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  } else {
    fail("Unsupported CNV file format: %s", ext)
  }

  if (is.data.frame(cnv)) cnv <- as.matrix(cnv)
  if (!is.matrix(cnv)) fail("CNV file must yield a matrix (genes x samples). Got class %s", class(cnv)[1])

  is_discrete <- all(stats::na.omit(as.numeric(cnv)) %in% c(-2, -1, 0, 1, 2))

  if (is_discrete) {
    amp_mat <- cnv >= focal_threshold
    del_mat <- cnv <= -focal_threshold
    high_amp_mat <- cnv >= 2L
    high_del_mat <- cnv <= -2L
  } else {
    amp_mat <- cnv >= amp_threshold
    del_mat <- cnv <= del_threshold
    high_amp_mat <- cnv >= 1
    high_del_mat <- cnv <= -1
  }

  freq_df <- data.frame(
    gene = rownames(cnv),
    n_samples = rowSums(!is.na(cnv)),
    amp_pct = rowMeans(amp_mat, na.rm = TRUE),
    del_pct = rowMeans(del_mat, na.rm = TRUE),
    high_amp_pct = rowMeans(high_amp_mat, na.rm = TRUE),
    high_del_pct = rowMeans(high_del_mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  freq_df$alt_pct <- freq_df$amp_pct + freq_df$del_pct
  freq_df <- freq_df[order(-freq_df$alt_pct), , drop = FALSE]

  top_amp <- utils::head(freq_df[order(-freq_df$amp_pct), ], top_n)
  top_del <- utils::head(freq_df[order(-freq_df$del_pct), ], top_n)

  sample_amp_counts <- colSums(amp_mat, na.rm = TRUE)
  sample_del_counts <- colSums(del_mat, na.rm = TRUE)
  per_sample <- data.frame(
    sample_id = colnames(cnv),
    amp_count = unname(sample_amp_counts),
    del_count = unname(sample_del_counts),
    cnv_burden = unname(sample_amp_counts + sample_del_counts),
    stringsAsFactors = FALSE
  )

  data <- tryCatch(standardize_project_data(project), error = function(e) NULL)
  survival_res <- NULL
  if (!is.null(data)) {
    per_sample$sample_id_short <- substr(per_sample$sample_id, 1, 15)
    clinical <- data$clinical
    clinical$sample_id_short <- substr(clinical$sample_id, 1, 15)
    merged <- merge(per_sample, clinical, by = "sample_id_short", suffixes = c("", ".clin"))
    merged$score_id <- merged$sample_id.clin %||% merged$sample_id
    if (nrow(merged) >= 30L) {
      survival_res <- tryCatch(
        run_survival_models(
          score_df = data.frame(sample_id = merged$score_id, score = merged$cnv_burden, stringsAsFactors = FALSE),
          clinical = data$clinical
        ),
        error = function(e) NULL
      )
    }
  }

  write_table_safe(freq_df, file.path(ctx$results_dir, sprintf("%s_cnv_gene_frequency.csv", project)))
  write_table_safe(top_amp, file.path(ctx$results_dir, sprintf("%s_cnv_top_amplifications.csv", project)))
  write_table_safe(top_del, file.path(ctx$results_dir, sprintf("%s_cnv_top_deletions.csv", project)))
  write_table_safe(per_sample, file.path(ctx$results_dir, sprintf("%s_cnv_per_sample.csv", project)))
  if (!is.null(survival_res)) {
    write_table_safe(survival_res$cox, file.path(ctx$results_dir, sprintf("%s_cnv_burden_cox.csv", project)))
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_top <- function(df, title, fill_col) {
      df <- df[!is.na(df[[fill_col]]) & df[[fill_col]] > 0, , drop = FALSE]
      if (!nrow(df)) return(NULL)
      df$gene <- factor(df$gene, levels = rev(df$gene))
      ggplot2::ggplot(df, ggplot2::aes(x = gene, y = .data[[fill_col]])) +
        ggplot2::geom_col(fill = if (fill_col == "amp_pct") "#C44E52" else "#4C72B0") +
        ggplot2::coord_flip() +
        ggplot2::labs(title = title, x = NULL, y = "Frequency") +
        theme_tcga()
    }
    p_amp <- plot_top(top_amp, sprintf("%s: top amplifications", project), "amp_pct")
    p_del <- plot_top(top_del, sprintf("%s: top deletions", project), "del_pct")
    if (!is.null(p_amp)) save_plot(p_amp, sprintf("%s_cnv_top_amplifications", project), ctx$plots_dir, width = 8, height = 7)
    if (!is.null(p_del)) save_plot(p_del, sprintf("%s_cnv_top_deletions", project), ctx$plots_dir, width = 8, height = 7)

    p_hist <- ggplot2::ggplot(per_sample, ggplot2::aes(x = cnv_burden)) +
      ggplot2::geom_histogram(bins = 40, fill = "#55A868", colour = "white") +
      ggplot2::labs(title = sprintf("%s CNV burden distribution", project),
                    x = "Altered genes per sample", y = "Samples") +
      theme_tcga()
    save_plot(p_hist, sprintf("%s_cnv_burden_histogram", project), ctx$plots_dir, width = 8, height = 5)
  }

  write_report(
    ctx,
    sprintf("CNV summary for %s", project),
    c(
      sprintf("- CNV file: `%s`", cnv_file),
      sprintf("- Values appear discrete (-2..+2): `%s`", is_discrete),
      sprintf("- Genes: `%s`", nrow(cnv)),
      sprintf("- Samples: `%s`", ncol(cnv)),
      sprintf("- Top amplified gene: `%s` (%s%%)", top_amp$gene[1], signif(top_amp$amp_pct[1] * 100, 3)),
      sprintf("- Top deleted gene: `%s` (%s%%)", top_del$gene[1], signif(top_del$del_pct[1] * 100, 3)),
      sprintf("- CNV-burden Cox available: `%s`", !is.null(survival_res)),
      "",
      "## Outputs",
      sprintf("- `results/%s_cnv_gene_frequency.csv`", project),
      sprintf("- `results/%s_cnv_top_amplifications.csv`", project),
      sprintf("- `results/%s_cnv_top_deletions.csv`", project),
      sprintf("- `results/%s_cnv_per_sample.csv`", project)
    )
  )

  list(frequency = freq_df, top_amp = top_amp, top_del = top_del,
       per_sample = per_sample, survival = survival_res)
}
