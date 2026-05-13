task_immune_phenotype <- function(config, ctx) {
  project <- config$project %||% fail("immune_phenotype requires a project.")
  data <- standardize_project_data(project)
  clinical <- data$clinical

  tme_methods <- unlist(config$tme_methods %||% c("estimate", "cibersort"))
  immune_subtype_file <- config$immune_subtype_file %||% file.path(project_root(), "0-Data", "PanCancer_subtypes.rda")

  tme_results <- list()
  if (requireNamespace("IOBR", quietly = TRUE)) {
    counts <- round(data$counts)
    dat <- tryCatch(
      IOBR::count2tpm(countMat = counts, idType = "Symbol", org = "hsa", source = "local"),
      error = function(e) {
        warning(sprintf("count2tpm failed: %s", conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (!is.null(dat)) {
      for (method in tme_methods) {
        fit <- tryCatch({
          if (identical(method, "cibersort")) {
            IOBR::deconvo_tme(eset = dat, method = method, arrays = FALSE, perm = as.integer(config$perm %||% 100))
          } else {
            IOBR::deconvo_tme(eset = dat, method = method, arrays = FALSE)
          }
        }, error = function(e) {
          warning(sprintf("%s failed: %s", method, conditionMessage(e)), call. = FALSE)
          NULL
        })
        if (!is.null(fit)) tme_results[[method]] <- as.data.frame(fit)
      }
    }
  } else {
    info("IOBR not installed; falling back to ssGSEA on immune gene sets.")
  }

  if (!length(tme_results) && requireNamespace("GSVA", quietly = TRUE) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    immune_gmt <- config$immune_gmt %||% file.path(project_root(), "1-Input", "Homo-HALLMARK-gmts.txt")
    if (file.exists(immune_gmt)) {
      gs_df <- clusterProfiler::read.gmt(immune_gmt)
      sets <- split(as.character(gs_df$gene), as.character(gs_df$term))
      keep <- grepl("immune|interferon|inflammatory|cytokine|TNF|complement|IL_", names(sets), ignore.case = TRUE)
      if (any(keep)) {
        sets <- sets[keep]
        expr <- log2(data$tpm + 1)
        rownames(expr) <- data$gene_map$gene_name %||% rownames(expr)
        score_mat <- tryCatch(
          if (utils::packageVersion("GSVA") >= "1.50.0") {
            GSVA::gsva(GSVA::ssgseaParam(expr, sets, minSize = 5, maxSize = 500), verbose = FALSE)
          } else {
            GSVA::gsva(expr, sets, method = "ssgsea", min.sz = 5, max.sz = 500, verbose = FALSE)
          },
          error = function(e) NULL
        )
        if (!is.null(score_mat)) {
          tme_results$immune_ssgsea <- data.frame(
            sample_id = colnames(score_mat),
            t(score_mat),
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (!length(tme_results)) {
    fail("immune_phenotype could not produce any score (IOBR not installed and no immune GMT available).")
  }

  joined <- tme_results[[1]]
  if (!"sample_id" %in% colnames(joined) && "ID" %in% colnames(joined)) {
    names(joined)[names(joined) == "ID"] <- "sample_id"
  }
  for (m in setdiff(names(tme_results), names(tme_results)[1])) {
    df <- tme_results[[m]]
    if ("ID" %in% colnames(df)) names(df)[names(df) == "ID"] <- "sample_id"
    if ("sample_id" %in% colnames(df)) {
      joined <- merge(joined, df, by = "sample_id", all = TRUE, suffixes = c("", sprintf(".%s", m)))
    }
  }

  immune_subtype <- NULL
  if (file.exists(immune_subtype_file)) {
    obj <- tryCatch(load_rdata_file(immune_subtype_file), error = function(e) NULL)
    if (!is.null(obj) && is.data.frame(obj) && "Immune.subtype" %in% colnames(obj)) {
      df <- obj[, c("pan.samplesID", "Immune.subtype"), drop = FALSE]
      df$sample_id_short <- substr(as.character(df$pan.samplesID), 1, 15)
      joined$sample_id_short <- substr(as.character(joined$sample_id), 1, 15)
      joined <- merge(joined, df[, c("sample_id_short", "Immune.subtype")], by = "sample_id_short", all.x = TRUE)
      joined$sample_id_short <- NULL
      immune_subtype <- table(joined$Immune.subtype, useNA = "ifany")
    }
  }

  numeric_cols <- vapply(joined, is.numeric, logical(1))
  numeric_cols["sample_id"] <- FALSE
  summary_df <- data.frame(
    metric = names(joined)[numeric_cols],
    median = vapply(joined[, numeric_cols, drop = FALSE], stats::median, FUN.VALUE = 0, na.rm = TRUE),
    mean = vapply(joined[, numeric_cols, drop = FALSE], mean, FUN.VALUE = 0, na.rm = TRUE),
    sd = vapply(joined[, numeric_cols, drop = FALSE], stats::sd, FUN.VALUE = 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  write_table_safe(joined, file.path(ctx$results_dir, sprintf("%s_immune_phenotype_per_sample.csv", project)))
  write_table_safe(summary_df, file.path(ctx$results_dir, sprintf("%s_immune_phenotype_summary.csv", project)))
  if (!is.null(immune_subtype)) {
    write_json(as.list(immune_subtype),
               file.path(ctx$results_dir, sprintf("%s_immune_subtype_counts.json", project)))
  }

  if (requireNamespace("ggplot2", quietly = TRUE) && "Immune.subtype" %in% colnames(joined)) {
    plot_df <- joined[!is.na(joined$Immune.subtype), , drop = FALSE]
    if (nrow(plot_df) && "StromalScore" %in% colnames(plot_df) && "ImmuneScore" %in% colnames(plot_df)) {
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = StromalScore, y = ImmuneScore, colour = Immune.subtype)) +
        ggplot2::geom_point(alpha = 0.7, size = 2) +
        ggplot2::labs(title = sprintf("%s immune phenotype", project)) +
        theme_tcga()
      save_plot(p, sprintf("%s_immune_phenotype_scatter", project), ctx$plots_dir, width = 8, height = 6)
    }
  }

  write_report(
    ctx,
    sprintf("Immune phenotype for %s", project),
    c(
      sprintf("- TME methods: `%s`", paste(names(tme_results), collapse = ", ")),
      sprintf("- Samples scored: `%s`", nrow(joined)),
      sprintf("- Immune subtype joined: `%s`", !is.null(immune_subtype)),
      "",
      "## Outputs",
      sprintf("- `results/%s_immune_phenotype_per_sample.csv`", project),
      sprintf("- `results/%s_immune_phenotype_summary.csv`", project)
    )
  )

  list(per_sample = joined, summary = summary_df, immune_subtype = immune_subtype)
}
