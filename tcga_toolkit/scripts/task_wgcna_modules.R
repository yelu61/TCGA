task_wgcna_modules <- function(config, ctx) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    fail("wgcna_modules requires the WGCNA package (install.packages('WGCNA')).")
  }
  project <- config$project %||% fail("wgcna_modules requires a project.")
  n_top_genes <- as.integer(config$n_top_genes %||% 5000L)
  min_module_size <- as.integer(config$min_module_size %||% 30L)
  power <- config$soft_power
  merge_cut_height <- as.numeric(config$merge_cut_height %||% 0.25)
  sample_filter <- unlist(config$sample_filter %||% character())
  trait_columns <- unlist(config$trait_columns %||% c("Group", "age_at_index", "ajcc_pathologic_stage"))

  data <- standardize_project_data(project)
  clinical <- data$clinical
  expr <- log2(data$tpm + 1)
  rownames(expr) <- data$gene_map$gene_name %||% rownames(expr)
  expr <- expr[!is.na(rownames(expr)) & nzchar(rownames(expr)), , drop = FALSE]
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

  if (length(sample_filter) && "Group" %in% colnames(clinical)) {
    keep_ids <- clinical$sample_id[as.character(clinical$Group) %in% sample_filter]
    if (length(keep_ids) >= 20L) {
      expr <- expr[, intersect(colnames(expr), keep_ids), drop = FALSE]
    }
  }

  variances <- apply(expr, 1, stats::var, na.rm = TRUE)
  keep <- order(variances, decreasing = TRUE)[seq_len(min(n_top_genes, length(variances)))]
  expr_top <- expr[keep, , drop = FALSE]
  datExpr <- t(expr_top)
  if (anyNA(datExpr)) {
    info("Replacing %s NA values in WGCNA input with column medians.", sum(is.na(datExpr)))
    datExpr <- apply(datExpr, 2, function(v) {
      if (any(is.na(v))) v[is.na(v)] <- stats::median(v, na.rm = TRUE)
      v
    })
  }

  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }

  if (is.null(power)) {
    powers <- c(1:10, seq(12, 20, by = 2))
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
    sft_tab <- as.data.frame(sft$fitIndices)
    candidate <- sft_tab$Power[sft_tab$SFT.R.sq >= 0.85 & sft_tab$slope < 0]
    power <- if (length(candidate)) min(candidate) else sft$powerEstimate %||% 6L
    write_table_safe(sft_tab, file.path(ctx$results_dir, sprintf("%s_wgcna_soft_threshold.csv", project)))
    info("Selected soft-thresholding power: %s", power)
  } else {
    power <- as.integer(power)
  }

  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    TOMType = "unsigned",
    minModuleSize = min_module_size,
    mergeCutHeight = merge_cut_height,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0,
    maxBlockSize = ncol(datExpr) + 1L
  )

  module_colors <- WGCNA::labels2colors(net$colors)
  module_df <- data.frame(
    gene = colnames(datExpr),
    module_label = net$colors,
    module_color = module_colors,
    stringsAsFactors = FALSE
  )
  module_sizes <- as.data.frame(table(module_df$module_color), stringsAsFactors = FALSE)
  names(module_sizes) <- c("module_color", "n_genes")
  module_sizes <- module_sizes[order(-module_sizes$n_genes), , drop = FALSE]

  module_eigengenes <- WGCNA::moduleEigengenes(datExpr, colors = module_colors)$eigengenes
  module_eigengenes$sample_id <- rownames(datExpr)
  module_eigengenes <- module_eigengenes[, c("sample_id", setdiff(colnames(module_eigengenes), "sample_id"))]

  trait_columns <- intersect(trait_columns, colnames(clinical))
  trait_corr <- NULL
  if (length(trait_columns)) {
    traits <- clinical[match(rownames(datExpr), clinical$sample_id), trait_columns, drop = FALSE]
    me_only <- module_eigengenes[, setdiff(colnames(module_eigengenes), "sample_id"), drop = FALSE]
    corr_rows <- list()
    for (tc in trait_columns) {
      v <- traits[[tc]]
      numeric_v <- suppressWarnings(as.numeric(as.character(v)))
      if (sum(!is.na(numeric_v)) < 20) {
        lv <- as.factor(v)
        if (length(levels(lv)) < 2L) next
        numeric_v <- as.integer(lv)
      }
      for (me in colnames(me_only)) {
        x <- me_only[[me]]
        ok <- !is.na(x) & !is.na(numeric_v)
        if (sum(ok) < 10) next
        test <- tryCatch(stats::cor.test(x[ok], numeric_v[ok], method = "spearman"),
                         error = function(e) NULL)
        if (is.null(test)) next
        corr_rows[[length(corr_rows) + 1L]] <- data.frame(
          module = me,
          trait = tc,
          correlation = unname(test$estimate),
          p_value = test$p.value,
          n = sum(ok),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(corr_rows)) {
      trait_corr <- do.call(rbind, corr_rows)
      trait_corr$padj <- stats::p.adjust(trait_corr$p_value, method = "BH")
      trait_corr <- trait_corr[order(trait_corr$padj), , drop = FALSE]
      write_table_safe(trait_corr, file.path(ctx$results_dir, sprintf("%s_wgcna_module_trait_correlation.csv", project)))
    }
  }

  write_table_safe(module_df, file.path(ctx$results_dir, sprintf("%s_wgcna_module_assignment.csv", project)))
  write_table_safe(module_sizes, file.path(ctx$results_dir, sprintf("%s_wgcna_module_sizes.csv", project)))
  write_table_safe(module_eigengenes, file.path(ctx$results_dir, sprintf("%s_wgcna_module_eigengenes.csv", project)))
  saveRDS(list(network = net, power = power, module_df = module_df),
          file.path(ctx$objects_dir, sprintf("%s_wgcna_network.rds", project)))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    sizes_df <- module_sizes[module_sizes$module_color != "grey", , drop = FALSE]
    if (nrow(sizes_df)) {
      sizes_df$module_color <- factor(sizes_df$module_color, levels = rev(sizes_df$module_color))
      p <- ggplot2::ggplot(sizes_df, ggplot2::aes(x = module_color, y = n_genes, fill = module_color)) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = sprintf("%s WGCNA module sizes", project), x = NULL, y = "Genes") +
        theme_tcga()
      save_plot(p, sprintf("%s_wgcna_module_sizes", project), ctx$plots_dir, width = 7, height = 6)
    }
  }

  write_report(
    ctx,
    sprintf("WGCNA modules for %s", project),
    c(
      sprintf("- Samples used: `%s`", nrow(datExpr)),
      sprintf("- Top-variance genes used: `%s`", ncol(datExpr)),
      sprintf("- Soft-thresholding power: `%s`", power),
      sprintf("- Min module size: `%s`", min_module_size),
      sprintf("- Merge cut height: `%s`", merge_cut_height),
      sprintf("- Modules detected (incl. grey): `%s`", nrow(module_sizes)),
      sprintf("- Largest non-grey module: `%s (%s genes)`",
              module_sizes$module_color[module_sizes$module_color != "grey"][1],
              module_sizes$n_genes[module_sizes$module_color != "grey"][1]),
      sprintf("- Module-trait correlations stored: `%s`", !is.null(trait_corr)),
      "",
      "## Outputs",
      sprintf("- `results/%s_wgcna_module_assignment.csv`", project),
      sprintf("- `results/%s_wgcna_module_sizes.csv`", project),
      sprintf("- `results/%s_wgcna_module_eigengenes.csv`", project),
      sprintf("- `objects/%s_wgcna_network.rds`", project)
    )
  )

  list(module_df = module_df, sizes = module_sizes, eigengenes = module_eigengenes,
       trait_correlation = trait_corr, power = power)
}
