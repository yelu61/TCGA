task_ssgsea_score <- function(config, ctx) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    fail("ssgsea_score requires the GSVA Bioconductor package.")
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    fail("ssgsea_score uses clusterProfiler::read.gmt to parse GMT files.")
  }

  project <- config$project %||% fail("ssgsea_score requires a project.")
  gmt_paths <- unlist(config$gmt_paths %||% list(
    file.path(project_root(), "1-Input", "Homo-HALLMARK-gmts.txt"),
    file.path(project_root(), "1-Input", "Homo-GOCC-gmts.txt"),
    file.path(project_root(), "1-Input", "Homo-GOMF-gmts.txt")
  ))
  method <- tolower(config$method %||% "ssgsea")
  if (!(method %in% c("ssgsea", "gsva", "zscore", "plage"))) {
    fail("ssgsea_score.method must be one of ssgsea/gsva/zscore/plage, got: %s", method)
  }
  min_size <- as.integer(config$min_size %||% 10L)
  max_size <- as.integer(config$max_size %||% 500L)

  data <- standardize_project_data(project)
  expr <- log2(data$tpm + 1)
  rownames(expr) <- data$gene_map$gene_name %||% rownames(expr)
  expr <- expr[!is.na(rownames(expr)) & nzchar(rownames(expr)), , drop = FALSE]
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

  all_scores <- list()
  per_collection <- list()

  for (gmt_path in gmt_paths) {
    if (!file.exists(gmt_path)) {
      warning(sprintf("Skipping missing GMT: %s", gmt_path), call. = FALSE)
      next
    }
    stem <- tools::file_path_sans_ext(basename(gmt_path))
    gs_df <- clusterProfiler::read.gmt(gmt_path)
    if (!nrow(gs_df)) next
    sets <- split(as.character(gs_df$gene), as.character(gs_df$term))
    sets <- sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
    if (!length(sets)) {
      warning(sprintf("No gene sets passed size filter for %s", stem), call. = FALSE)
      next
    }

    info("ssgsea_score: %s sets from %s", length(sets), stem)
    score_mat <- tryCatch(
      {
        if (utils::packageVersion("GSVA") >= "1.50.0") {
          param <- switch(
            method,
            ssgsea = GSVA::ssgseaParam(expr, sets, minSize = min_size, maxSize = max_size),
            gsva = GSVA::gsvaParam(expr, sets, minSize = min_size, maxSize = max_size),
            zscore = GSVA::zscoreParam(expr, sets, minSize = min_size, maxSize = max_size),
            plage = GSVA::plageParam(expr, sets, minSize = min_size, maxSize = max_size)
          )
          GSVA::gsva(param, verbose = FALSE)
        } else {
          GSVA::gsva(expr, sets, method = method, kcdf = "Gaussian",
                     min.sz = min_size, max.sz = max_size, verbose = FALSE)
        }
      },
      error = function(e) {
        warning(sprintf("GSVA failed for %s: %s", stem, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(score_mat)) next

    score_df <- as.data.frame(t(score_mat), check.names = FALSE)
    score_df$sample_id <- rownames(score_df)
    score_df <- score_df[, c("sample_id", setdiff(colnames(score_df), "sample_id")), drop = FALSE]
    out_path <- file.path(ctx$results_dir, sprintf("%s_%s_%s.csv", project, stem, method))
    write_table_safe(score_df, out_path)
    per_collection[[stem]] <- score_df
    all_scores[[length(all_scores) + 1L]] <- list(collection = stem, n_sets = nrow(score_mat))
  }

  if (!length(all_scores)) {
    fail("No GMT collections produced scores. Check gmt_paths and min/max size filters.")
  }

  summary_lines <- c(
    sprintf("- Method: `%s`", method),
    sprintf("- Sample count: `%s`", ncol(expr)),
    sprintf("- Collections scored: `%s`", length(all_scores)),
    "",
    "## Outputs",
    paste0("- `results/", project, "_<collection>_", method, ".csv` (one per collection)")
  )
  for (item in all_scores) {
    summary_lines <- c(summary_lines,
                       sprintf("- %s: `%s` gene sets", item$collection, item$n_sets))
  }

  write_report(
    ctx,
    sprintf("ssGSEA / GSVA scoring for %s", project),
    summary_lines
  )

  list(collections = per_collection, summary = all_scores)
}
