task_gene_correlation_heatmap <- function(config, ctx) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    fail("gene_correlation_heatmap requires ComplexHeatmap and circlize.")
  }

  target_gene <- config$target_gene %||% fail("gene_correlation_heatmap requires target_gene.")
  gene_list_file <- config$gene_list_file %||% fail("gene_correlation_heatmap requires gene_list_file.")
  gene_table <- read_gene_list_file(
    path = gene_list_file,
    sheet = config$sheet %||% NULL,
    gene_column = config$gene_column %||% NULL
  )
  genes <- unique(gene_table$gene)
  genes <- setdiff(genes, target_gene)
  projects <- unlist(config$projects %||% available_projects(include_target = FALSE))
  method <- config$method %||% "pearson"
  sample_filter <- unlist(config$sample_filter %||% c("Tumor"))
  top_n_genes <- as.integer(config$top_n_genes %||% length(genes))
  transpose_heatmap <- isTRUE(config$transpose_heatmap %||% FALSE)
  cluster_rows_flag <- isTRUE(config$cluster_rows %||% TRUE)
  cluster_columns_flag <- isTRUE(config$cluster_columns %||% TRUE)
  color_scheme <- config$color_scheme %||% "green_black_red"

  result_rows <- list()
  coverage_rows <- list()
  for (project in projects) {
    info("Processing gene correlation heatmap: %s", project)
    data <- tryCatch(standardize_project_data(project), error = function(e) NULL)
    if (is.null(data)) {
      next
    }
    cor_df <- correlate_target_with_gene_set(
      data,
      target_gene = target_gene,
      genes = genes,
      method = method,
      sample_filter = sample_filter
    )
    cor_df$project <- project
    cor_df$project_code <- tcga_project_code(project)
    result_rows[[length(result_rows) + 1L]] <- cor_df
    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      project = project,
      matched_gene_count = sum(!is.na(cor_df$correlation)),
      requested_gene_count = length(genes),
      stringsAsFactors = FALSE
    )
  }
  if (!length(result_rows)) {
    fail("No project correlations were generated.")
  }

  long_df <- do.call(rbind, result_rows)
  coverage_df <- do.call(rbind, coverage_rows)
  long_df <- long_df[order(long_df$project_code, long_df$gene), , drop = FALSE]

  gene_rank <- stats::aggregate(
    correlation ~ gene,
    data = long_df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  gene_rank$correlation[is.nan(gene_rank$correlation)] <- NA_real_
  gene_rank <- gene_rank[order(-gene_rank$correlation, gene_rank$gene, na.last = TRUE), , drop = FALSE]
  selected_genes <- head(gene_rank$gene[!is.na(gene_rank$correlation)], top_n_genes)
  if (!length(selected_genes)) {
    fail("No genes with valid correlations were found for %s.", target_gene)
  }
  long_df <- long_df[long_df$gene %in% selected_genes, , drop = FALSE]

  projects_order <- unique(long_df$project_code)
  genes_order <- selected_genes
  cor_mat <- matrix(NA_real_, nrow = length(projects_order), ncol = length(genes_order), dimnames = list(projects_order, genes_order))
  p_mat <- cor_mat
  n_mat <- cor_mat
  for (i in seq_len(nrow(long_df))) {
    row <- long_df[i, ]
    cor_mat[row$project_code, row$gene] <- row$correlation
    p_mat[row$project_code, row$gene] <- row$p_value
    n_mat[row$project_code, row$gene] <- row$n
  }

  cluster_mat <- cor_mat
  cluster_mat[is.na(cluster_mat)] <- 0
  if (transpose_heatmap) {
    display_mat <- t(cor_mat)
    cluster_mat_disp <- t(cluster_mat)
  } else {
    display_mat <- cor_mat
    cluster_mat_disp <- cluster_mat
  }

  row_cluster <- if (cluster_rows_flag && nrow(cluster_mat_disp) >= 2) stats::hclust(stats::dist(cluster_mat_disp)) else FALSE
  col_cluster <- if (cluster_columns_flag && ncol(cluster_mat_disp) >= 2) stats::hclust(stats::dist(t(cluster_mat_disp))) else FALSE

  colors <- switch(
    color_scheme,
    red_blue = circlize::colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")),
    blue_red = circlize::colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")),
    green_black_red = circlize::colorRamp2(c(-0.6, 0, 0.6), c("red", "black", "green")),
    circlize::colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))
  )
  ht <- ComplexHeatmap::Heatmap(
    display_mat,
    name = "R",
    col = colors,
    cluster_rows = row_cluster,
    cluster_columns = col_cluster,
    na_col = "#F0F0F0",
    row_names_side = "left",
    column_names_rot = 90,
    heatmap_legend_param = list(at = c(-0.5, 0, 0.5))
  )

  pdf_path <- file.path(ctx$plots_dir, sprintf("%s_gene_correlation_heatmap.pdf", target_gene))
  png_path <- file.path(ctx$plots_dir, sprintf("%s_gene_correlation_heatmap.png", target_gene))
  grDevices::pdf(pdf_path, width = max(10, ncol(display_mat) * 0.25), height = max(8, nrow(display_mat) * 0.25))
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
  grDevices::png(png_path, width = max(1600, ncol(display_mat) * 40), height = max(1200, nrow(display_mat) * 40), res = 150)
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()

  cor_df_out <- data.frame(project_code = rownames(cor_mat), as.data.frame(cor_mat, check.names = FALSE), check.names = FALSE)
  p_df_out <- data.frame(project_code = rownames(p_mat), as.data.frame(p_mat, check.names = FALSE), check.names = FALSE)
  n_df_out <- data.frame(project_code = rownames(n_mat), as.data.frame(n_mat, check.names = FALSE), check.names = FALSE)

  write_table_safe(long_df, file.path(ctx$results_dir, sprintf("%s_gene_correlation_long.csv", target_gene)))
  write_table_safe(cor_df_out, file.path(ctx$results_dir, sprintf("%s_gene_correlation_matrix.csv", target_gene)))
  write_table_safe(p_df_out, file.path(ctx$results_dir, sprintf("%s_gene_correlation_pvalues.csv", target_gene)))
  write_table_safe(n_df_out, file.path(ctx$results_dir, sprintf("%s_gene_correlation_n.csv", target_gene)))
  write_table_safe(coverage_df, file.path(ctx$results_dir, sprintf("%s_gene_correlation_coverage.csv", target_gene)))
  write_table_safe(
    data.frame(
      gene = genes,
      matched_any_project = genes %in% unique(long_df$gene[!is.na(long_df$correlation)]),
      stringsAsFactors = FALSE
    ),
    file.path(ctx$results_dir, sprintf("%s_gene_match_status.csv", target_gene))
  )

  write_report(
    ctx,
    sprintf("Gene Correlation Heatmap for %s", target_gene),
    c(
      sprintf("- Gene list file: `%s`", gene_list_file),
      sprintf("- Requested genes: `%s`", length(genes)),
      sprintf("- Genes displayed: `%s`", length(selected_genes)),
      sprintf("- Projects evaluated: `%s`", length(projects_order)),
      sprintf("- Correlation method: `%s`", method),
      sprintf("- Sample filter: `%s`", paste(sample_filter, collapse = ", ")),
      sprintf("- Heatmap transposed: `%s`", transpose_heatmap),
      sprintf("- Color scheme: `%s`", color_scheme),
      "",
      "## Outputs",
      sprintf("- `plots/%s_gene_correlation_heatmap.pdf`", target_gene),
      sprintf("- `plots/%s_gene_correlation_heatmap.png`", target_gene),
      sprintf("- `results/%s_gene_correlation_matrix.csv`", target_gene),
      sprintf("- `results/%s_gene_correlation_pvalues.csv`", target_gene),
      sprintf("- `results/%s_gene_correlation_long.csv`", target_gene)
    )
  )

  list(correlation = cor_mat, p_values = p_mat, n = n_mat, long = long_df)
}
