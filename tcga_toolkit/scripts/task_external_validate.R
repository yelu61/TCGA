task_external_validate <- function(config, ctx) {
  cohort_name <- config$cohort_name %||% fail("external_validate requires cohort_name (e.g. CGGA, METABRIC, GSE12345).")
  expression_file <- config$expression_file %||% fail("external_validate requires expression_file (gene x sample matrix).")
  clinical_file <- config$clinical_file %||% fail("external_validate requires clinical_file with survival columns.")
  if (!file.exists(expression_file)) fail("expression_file not found: %s", expression_file)
  if (!file.exists(clinical_file)) fail("clinical_file not found: %s", clinical_file)

  features <- unlist(config$feature_genes %||% character())
  if (!length(features) && !is.null(config$signature_file)) {
    sig <- read_marker_table(config$signature_file)
    features <- unique(sig$gene)
  }
  # Locked model metadata from a prognostic_model run: features, coefficients,
  # transform and cutoff are reused as-is instead of being re-derived on the
  # external cohort (re-deriving invalidates the validation claim).
  model_meta <- NULL
  if (!is.null(config$model_meta_file) && file.exists(config$model_meta_file)) {
    model_meta <- jsonlite::read_json(config$model_meta_file, simplifyVector = TRUE)
  }
  weights <- NULL
  if (!is.null(config$weight_file) && nzchar(config$weight_file) && file.exists(config$weight_file)) {
    w <- utils::read.csv(config$weight_file, stringsAsFactors = FALSE, check.names = FALSE)
    gene_col <- intersect(c("gene_name", "gene", "Gene", "feature_id"), colnames(w))[1]
    wt_col <- intersect(c("coefficient", "weight", "estimate"), colnames(w))[1]
    if (!is.na(gene_col) && !is.na(wt_col)) {
      weights <- stats::setNames(as.numeric(w[[wt_col]]), as.character(w[[gene_col]]))
      features <- unique(c(features, names(weights)))
    }
  }
  if (is.null(weights) && !is.null(model_meta) && length(model_meta$coefficients)) {
    weights <- as.numeric(unlist(model_meta$coefficients))
    names(weights) <- names(model_meta$coefficients)
  }
  if (!length(features) && !is.null(model_meta) && length(model_meta$features)) {
    features <- unlist(model_meta$features)
  }
  if (!length(features)) {
    fail("external_validate requires feature_genes, signature_file, weight_file, or model_meta_file with features.")
  }

  read_matrix <- function(path) {
    e <- tolower(tools::file_ext(path))
    if (e %in% c("rda", "rdata")) {
      x <- load_rdata_file(path)
    } else if (e == "tsv") {
      x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    } else if (e == "csv") {
      x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    } else if (e == "gz") {
      x <- utils::read.delim(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    } else {
      fail("Unsupported expression format: %s", e)
    }
    if (is.data.frame(x)) x <- as.matrix(x)
    x
  }
  expr <- read_matrix(expression_file)
  expr <- apply(expr, 2, as.numeric)
  rownames(expr) <- rownames(read_matrix(expression_file))
  if (anyNA(rownames(expr))) {
    fail("expression_file must use gene symbols (or IDs) as row names.")
  }
  # Transform: when the model metadata declares a locked transform, honour it;
  # otherwise fall back to the scale heuristic and say so in the report.
  transform_note <- "heuristic (max>100 -> log2(x+1))"
  if (!is.null(model_meta) && identical(as.character(model_meta$transform), "log2p1")) {
    if (max(expr, na.rm = TRUE) > 100) {
      info("Applying locked transform log2(x+1) from model metadata.")
      expr <- log2(expr + 1)
      transform_note <- "locked from training: log2(x+1)"
    } else {
      info("Model metadata declares log2(x+1); external matrix already looks log-scale, using as-is.")
      transform_note <- "locked from training: log2(x+1) (external already log-scale)"
    }
  } else if (max(expr, na.rm = TRUE) > 100) {
    info("Expression looks like raw TPM/counts; applying log2(x+1).")
    expr <- log2(expr + 1)
  }

  e <- tolower(tools::file_ext(clinical_file))
  clinical <- if (e %in% c("tsv", "txt")) {
    utils::read.delim(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (e == "csv") {
    utils::read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (e %in% c("rda", "rdata")) {
    obj <- load_rdata_file(clinical_file)
    if (!is.data.frame(obj)) fail("clinical_file rda must contain a data.frame.")
    obj
  } else {
    fail("Unsupported clinical format: %s", e)
  }

  id_col <- intersect(c("sample_id", "Sample", "sampleID", "barcode", "Patient", "patient_id"), colnames(clinical))[1]
  if (is.na(id_col)) {
    fail("clinical_file must have one of: sample_id, Sample, sampleID, barcode, Patient.")
  }
  time_col <- intersect(c("OS.time", "os_time", "survival_time", "days_to_event", "time"),
                       colnames(clinical))[1]
  event_col <- intersect(c("OS", "os", "event", "vital_status", "status", "survival_event"),
                        colnames(clinical))[1]
  if (is.na(time_col) || is.na(event_col)) {
    fail("clinical_file must include a survival time column (e.g. OS.time) and event column (e.g. OS).")
  }

  clinical$sample_id <- as.character(clinical[[id_col]])
  clinical$survival_time <- suppressWarnings(as.numeric(clinical[[time_col]]))
  clinical$survival_event <- coerce_status(clinical[[event_col]])

  matched <- intersect(features, rownames(expr))
  if (length(matched) < 2L) {
    fail("Only %s of %s features matched the external expression matrix; aborting.", length(matched), length(features))
  }
  info("Matched %s / %s features in external cohort %s.", length(matched), length(features), cohort_name)

  expr_sub <- expr[matched, , drop = FALSE]
  if (!is.null(weights)) {
    w_vec <- weights[matched]
    w_vec[is.na(w_vec)] <- 0
    risk <- as.numeric(colSums(expr_sub * w_vec, na.rm = TRUE))
  } else {
    expr_z <- t(scale(t(expr_sub)))
    expr_z[is.na(expr_z)] <- 0
    risk <- colMeans(expr_z, na.rm = TRUE)
  }

  risk_df <- data.frame(
    sample_id = colnames(expr),
    risk_score = risk,
    stringsAsFactors = FALSE
  )
  merged <- merge(risk_df, clinical, by = "sample_id")
  merged <- merged[!is.na(merged$survival_time) & !is.na(merged$survival_event), , drop = FALSE]
  if (nrow(merged) < 20L) {
    fail("Only %s samples have both risk score and survival; aborting.", nrow(merged))
  }
  # Risk grouping cutoff: locked from training (config$cutoff or model
  # metadata) whenever the risk score is on the training scale; re-deriving a
  # cutoff on the external cohort turns validation into exploration.
  locked_cutoff <- suppressWarnings(as.numeric(config$cutoff %||% model_meta$cutoff %||% NA))
  if (!is.na(locked_cutoff) && is.null(weights)) {
    info("Ignoring locked cutoff: no locked coefficients supplied, risk score is a re-derived z-score mean.")
    locked_cutoff <- NA_real_
  }
  if (!is.na(locked_cutoff)) {
    cutoff <- locked_cutoff
    cutoff_source <- "locked from training"
  } else {
    cutoff <- stats::median(merged$risk_score, na.rm = TRUE)
    cutoff_source <- "external-cohort median (EXPLORATORY: no locked cutoff supplied)"
  }
  merged$risk_group <- ifelse(merged$risk_score >= cutoff, "High", "Low")

  cox <- NULL
  km_summary <- NULL
  if (requireNamespace("survival", quietly = TRUE)) {
    surv <- survival::Surv(merged$survival_time, merged$survival_event)
    fit <- tryCatch(survival::coxph(surv ~ risk_score, data = merged), error = function(e) NULL)
    if (!is.null(fit)) {
      sm <- summary(fit)
      cox <- data.frame(
        cohort = cohort_name,
        n = nrow(merged),
        events = sum(merged$survival_event == 1, na.rm = TRUE),
        hazard_ratio = unname(sm$coefficients[, "exp(coef)"]),
        conf_low = unname(sm$conf.int[, "lower .95"]),
        conf_high = unname(sm$conf.int[, "upper .95"]),
        p_value = unname(sm$coefficients[, "Pr(>|z|)"]),
        concordance = unname(sm$concordance[1]),
        stringsAsFactors = FALSE
      )
    }
    km_fit <- tryCatch(
      survival::survfit(survival::Surv(survival_time, survival_event) ~ risk_group, data = merged),
      error = function(e) NULL
    )
    if (!is.null(km_fit)) {
      km_summary <- list(
        cutoff = cutoff,
        cutoff_source = cutoff_source,
        n = nrow(merged),
        groups = as.list(table(merged$risk_group))
      )
    }
  }

  write_table_safe(merged, file.path(ctx$results_dir, sprintf("%s_external_risk_scores.csv", cohort_name)))
  if (!is.null(cox)) {
    write_table_safe(cox, file.path(ctx$results_dir, sprintf("%s_external_cox.csv", cohort_name)))
  }
  if (!is.null(km_summary)) {
    write_json(km_summary, file.path(ctx$results_dir, sprintf("%s_external_km.json", cohort_name)))
  }

  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("survival", quietly = TRUE)) {
    p_km <- tryCatch(
      plot_km_curve(
        score_df = data.frame(sample_id = merged$sample_id, score = merged$risk_score, stringsAsFactors = FALSE),
        clinical = data.frame(
          sample_id = merged$sample_id,
          vital_status = ifelse(merged$survival_event == 1, "dead", "alive"),
          days_to_death = ifelse(merged$survival_event == 1, merged$survival_time, NA),
          days_to_last_follow_up = ifelse(merged$survival_event == 0, merged$survival_time, NA),
          stringsAsFactors = FALSE
        ),
        title = sprintf("External validation in %s", cohort_name)
      ),
      error = function(e) NULL
    )
    if (!is.null(p_km)) {
      save_plot(p_km, sprintf("%s_external_km", cohort_name), ctx$plots_dir, width = 7, height = 6)
    }
  }

  write_report(
    ctx,
    sprintf("External validation in %s", cohort_name),
    c(
      sprintf("- Expression matrix: `%s`", expression_file),
      sprintf("- Clinical: `%s`", clinical_file),
      sprintf("- Features requested / matched: `%s / %s`", length(features), length(matched)),
      sprintf("- Risk derived via: `%s`", ifelse(is.null(weights), "z-score mean (re-derived on external cohort)", "locked weighted sum")),
      sprintf("- Transform: `%s`", transform_note),
      sprintf("- Risk cutoff: `%s` (%s)", signif(cutoff, 4), cutoff_source),
      sprintf("- Samples with survival: `%s` (events=%s)", nrow(merged), sum(merged$survival_event == 1, na.rm = TRUE)),
      sprintf("- Cox HR (95%% CI): `%s (%s - %s)` p=`%s`",
              if (!is.null(cox)) signif(cox$hazard_ratio, 3) else "NA",
              if (!is.null(cox)) signif(cox$conf_low, 3) else "NA",
              if (!is.null(cox)) signif(cox$conf_high, 3) else "NA",
              if (!is.null(cox)) signif(cox$p_value, 3) else "NA"),
      "",
      "## Outputs",
      sprintf("- `results/%s_external_risk_scores.csv`", cohort_name),
      sprintf("- `results/%s_external_cox.csv`", cohort_name)
    )
  )

  list(risk = merged, cox = cox, matched_features = matched, cohort = cohort_name)
}
