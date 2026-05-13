task_run_enrichment <- function(config, ctx) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    fail("run_enrichment requires the clusterProfiler package.")
  }

  input_deg <- config$input_deg %||% fail("run_enrichment requires input_deg.")
  mode <- config$mode %||% "ora"
  gmt_paths <- config$gmt_paths %||% list(
    file.path(project_root(), "1-Input", "Homo-HALLMARK-gmts.txt"),
    file.path(project_root(), "1-Input", "Homo-GOCC-gmts.txt"),
    file.path(project_root(), "1-Input", "Homo-GOMF-gmts.txt")
  )

  deg <- utils::read.csv(input_deg, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("gene_name" %in% colnames(deg))) {
    fail("DEG table must contain gene_name column.")
  }

  results <- list()
  for (gmt_path in unlist(gmt_paths)) {
    if (!file.exists(gmt_path)) {
      next
    }
    gmt <- clusterProfiler::read.gmt(gmt_path)
    stem <- tools::file_path_sans_ext(basename(gmt_path))

    if (identical(mode, "gsea")) {
      ranked <- deg$log2FoldChange
      names(ranked) <- deg$gene_name
      ranked <- ranked[!is.na(ranked)]
      ranked <- sort(ranked, decreasing = TRUE)
      fit <- suppressWarnings(clusterProfiler::GSEA(geneList = ranked, TERM2GENE = gmt[, 1:2], pvalueCutoff = 1))
      out <- if (nrow(as.data.frame(fit))) as.data.frame(fit) else data.frame()
    } else {
      padj_cutoff <- as.numeric(config$padj_cutoff %||% 0.05)
      log2fc_cutoff <- as.numeric(config$log2fc_cutoff %||% 1)
      selected <- unique(deg$gene_name[!is.na(deg$padj) & deg$padj < padj_cutoff & abs(deg$log2FoldChange) >= log2fc_cutoff])
      fit <- suppressWarnings(clusterProfiler::enricher(gene = selected, TERM2GENE = gmt[, 1:2], pvalueCutoff = 1))
      out <- if (nrow(as.data.frame(fit))) as.data.frame(fit) else data.frame()
    }

    out_path <- file.path(ctx$results_dir, sprintf("%s_%s_enrichment.csv", stem, mode))
    write_table_safe(out, out_path)
    results[[stem]] <- out

    p_dot <- plot_dotplot_enrichment(out, title = sprintf("%s (%s)", stem, mode))
    if (!is.null(p_dot)) {
      save_plot(p_dot, sprintf("%s_%s_dotplot", stem, mode), ctx$plots_dir)
    }
  }

  write_report(
    ctx,
    sprintf("Enrichment (%s)", mode),
    c(
      sprintf("- DEG input: `%s`", input_deg),
      sprintf("- GMT collections processed: `%s`", length(results)),
      "",
      "## Outputs",
      "- One CSV per GMT collection in `results/`"
    )
  )

  results
}
