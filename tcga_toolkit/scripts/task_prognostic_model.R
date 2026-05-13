task_prognostic_model <- function(config, ctx) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    fail("prognostic_model requires the survival package.")
  }
  method <- tolower(config$method %||% "lasso_cox")
  if (!(method %in% c("lasso_cox", "rsf"))) {
    fail("prognostic_model.method must be lasso_cox or rsf, got: %s", method)
  }
  if (method == "lasso_cox" && !requireNamespace("glmnet", quietly = TRUE)) {
    fail("lasso_cox requires the glmnet package.")
  }
  if (method == "rsf" && !requireNamespace("randomForestSRC", quietly = TRUE)) {
    fail("rsf requires the randomForestSRC package.")
  }

  project <- config$project %||% fail("prognostic_model requires a project.")
  data <- standardize_project_data(project)
  clinical <- data$clinical

  features <- unique(unlist(config$feature_genes %||% character()))
  if (!length(features) && !is.null(config$signature_file)) {
    sig <- read_marker_table(config$signature_file)
    features <- unique(sig$gene)
  }
  if (!length(features) && !is.null(config$gene_list_file)) {
    gl <- read_gene_list_file(
      config$gene_list_file,
      sheet = config$sheet,
      gene_column = config$gene_column
    )
    features <- unique(gl$gene)
  }
  if (!length(features)) {
    fail("prognostic_model requires feature_genes, signature_file or gene_list_file.")
  }

  feature_ids <- vapply(features, function(g) {
    fid <- resolve_feature_id(g, data$tpm, data$gene_map)
    if (is.null(fid)) NA_character_ else fid
  }, character(1))
  feature_ids <- feature_ids[!is.na(feature_ids)]
  if (length(feature_ids) < 2L) {
    fail("Need at least 2 features that resolve in %s; resolved %s.", project, length(feature_ids))
  }
  mat <- log2(data$tpm[feature_ids, , drop = FALSE] + 1)
  if (any(is.na(mat))) mat[is.na(mat)] <- 0
  feature_df <- as.data.frame(t(mat), stringsAsFactors = FALSE)
  feature_df$sample_id <- rownames(feature_df)

  cols <- survival_columns(clinical)
  if (is.na(cols$status %||% NA)) {
    fail("Clinical table for %s has no vital_status / status column.", project)
  }

  time_to_death <- suppressWarnings(as.numeric(clinical[[cols$time_death]]))
  if (!is.na(cols$time_followup) && cols$time_followup %in% colnames(clinical)) {
    followup <- suppressWarnings(as.numeric(clinical[[cols$time_followup]]))
    time_to_death <- ifelse(is.na(time_to_death), followup, time_to_death)
  }
  surv_df <- data.frame(
    sample_id = clinical$sample_id,
    survival_time = time_to_death,
    survival_event = coerce_status(clinical[[cols$status]]),
    stringsAsFactors = FALSE
  )

  covariates <- unlist(config$covariates %||% character())
  covariates <- covariates[covariates %in% colnames(clinical)]
  if (length(covariates)) {
    cov_df <- clinical[, c("sample_id", covariates), drop = FALSE]
    surv_df <- merge(surv_df, cov_df, by = "sample_id", all.x = TRUE)
  }

  if ("Group" %in% colnames(clinical)) {
    tumor_ids <- clinical$sample_id[as.character(clinical$Group) %in% c("Tumor")]
    if (length(tumor_ids) >= 30L) {
      surv_df <- surv_df[surv_df$sample_id %in% tumor_ids, , drop = FALSE]
    }
  }

  dat <- merge(surv_df, feature_df, by = "sample_id")
  dat <- dat[!is.na(dat$survival_time) & !is.na(dat$survival_event) & dat$survival_time > 0, , drop = FALSE]
  if (nrow(dat) < 30L) {
    fail("Only %s samples have complete survival + expression data; need at least 30.", nrow(dat))
  }

  seed <- as.integer(config$seed %||% 20260512L)
  set.seed(seed)
  test_fraction <- as.numeric(config$test_fraction %||% 0.3)
  if (test_fraction <= 0 || test_fraction >= 1) test_fraction <- 0.3
  n <- nrow(dat)
  idx_test <- sample.int(n, size = max(5L, floor(n * test_fraction)))
  train <- dat[-idx_test, , drop = FALSE]
  test <- dat[idx_test, , drop = FALSE]

  surv_train <- survival::Surv(train$survival_time, train$survival_event)
  surv_test <- survival::Surv(test$survival_time, test$survival_event)

  feature_cols <- intersect(colnames(feature_df), colnames(train))
  feature_cols <- setdiff(feature_cols, "sample_id")
  x_train <- as.matrix(train[, feature_cols, drop = FALSE])
  x_test <- as.matrix(test[, feature_cols, drop = FALSE])

  coef_df <- NULL
  risk_train <- NULL
  risk_test <- NULL
  model_obj <- NULL
  method_summary <- character()

  if (method == "lasso_cox") {
    cvfit <- glmnet::cv.glmnet(x_train, surv_train, family = "cox", nfolds = 10)
    chosen_lambda <- cvfit$lambda.min
    coefs <- as.matrix(stats::coef(cvfit, s = chosen_lambda))
    coef_df <- data.frame(
      feature_id = rownames(coefs),
      coefficient = unname(coefs[, 1]),
      stringsAsFactors = FALSE
    )
    coef_df <- coef_df[coef_df$coefficient != 0, , drop = FALSE]
    coef_df <- merge(coef_df, data$gene_map[, c("feature_id", "gene_id", "gene_name")], by = "feature_id", all.x = TRUE)
    coef_df <- coef_df[order(-abs(coef_df$coefficient)), , drop = FALSE]

    risk_train <- as.numeric(stats::predict(cvfit, newx = x_train, s = chosen_lambda, type = "link"))
    risk_test <- as.numeric(stats::predict(cvfit, newx = x_test, s = chosen_lambda, type = "link"))
    model_obj <- cvfit
    method_summary <- c(
      sprintf("- Method: Lasso-Cox (cv.glmnet, nfolds=10)"),
      sprintf("- lambda.min: `%s`", signif(chosen_lambda, 4)),
      sprintf("- Non-zero features: `%s / %s`", nrow(coef_df), length(feature_cols))
    )
  } else {
    train_df <- train[, c("survival_time", "survival_event", feature_cols), drop = FALSE]
    rsf <- randomForestSRC::rfsrc(
      survival::Surv(survival_time, survival_event) ~ .,
      data = train_df,
      ntree = as.integer(config$ntree %||% 500L),
      nodesize = as.integer(config$nodesize %||% 15L),
      importance = TRUE
    )
    risk_train <- rsf$predicted.oob
    pred_test <- stats::predict(rsf, newdata = test[, feature_cols, drop = FALSE])
    risk_test <- pred_test$predicted
    importance <- as.numeric(rsf$importance)
    coef_df <- data.frame(
      feature_id = feature_cols,
      importance = importance,
      stringsAsFactors = FALSE
    )
    coef_df <- merge(coef_df, data$gene_map[, c("feature_id", "gene_id", "gene_name")], by = "feature_id", all.x = TRUE)
    coef_df <- coef_df[order(-coef_df$importance), , drop = FALSE]
    model_obj <- rsf
    method_summary <- c(
      sprintf("- Method: Random Survival Forest (ntree=%s)", rsf$ntree),
      sprintf("- OOB error rate: `%s`", signif(rsf$err.rate[length(rsf$err.rate)], 4))
    )
  }

  risk_df <- rbind(
    data.frame(sample_id = train$sample_id, split = "train", risk_score = risk_train,
               survival_time = train$survival_time, survival_event = train$survival_event,
               stringsAsFactors = FALSE),
    data.frame(sample_id = test$sample_id, split = "test", risk_score = risk_test,
               survival_time = test$survival_time, survival_event = test$survival_event,
               stringsAsFactors = FALSE)
  )

  cindex <- function(risk, time, event) {
    if (sum(event == 1, na.rm = TRUE) < 5) return(NA_real_)
    fit <- tryCatch(
      survival::coxph(survival::Surv(time, event) ~ risk),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    unname(summary(fit)$concordance[1])
  }
  cindex_train <- cindex(risk_train, train$survival_time, train$survival_event)
  cindex_test <- cindex(risk_test, test$survival_time, test$survival_event)

  perf_df <- data.frame(
    split = c("train", "test"),
    n = c(nrow(train), nrow(test)),
    events = c(sum(train$survival_event, na.rm = TRUE), sum(test$survival_event, na.rm = TRUE)),
    cindex = c(cindex_train, cindex_test),
    stringsAsFactors = FALSE
  )

  if (requireNamespace("timeROC", quietly = TRUE) && method == "lasso_cox") {
    times <- as.numeric(config$time_points %||% c(365, 1095, 1825))
    times <- times[times > 0 & times < max(dat$survival_time, na.rm = TRUE)]
    if (length(times)) {
      roc_test <- tryCatch(
        timeROC::timeROC(
          T = test$survival_time, delta = test$survival_event, marker = risk_test,
          cause = 1, times = times, iid = FALSE
        ),
        error = function(e) NULL
      )
      if (!is.null(roc_test)) {
        auc_df <- data.frame(
          time_days = times,
          AUC_test = as.numeric(roc_test$AUC),
          stringsAsFactors = FALSE
        )
        write_table_safe(auc_df, file.path(ctx$results_dir, sprintf("%s_prognostic_timeROC_auc.csv", project)))
        method_summary <- c(method_summary,
                            sprintf("- timeROC AUC (test): %s",
                                    paste(sprintf("t=%s d, AUC=%s", times, signif(auc_df$AUC_test, 3)),
                                          collapse = "; ")))
      }
    }
  }

  write_table_safe(coef_df, file.path(ctx$results_dir, sprintf("%s_prognostic_coefficients.csv", project)))
  write_table_safe(risk_df, file.path(ctx$results_dir, sprintf("%s_prognostic_risk_scores.csv", project)))
  write_table_safe(perf_df, file.path(ctx$results_dir, sprintf("%s_prognostic_performance.csv", project)))
  saveRDS(model_obj, file.path(ctx$objects_dir, sprintf("%s_prognostic_model.rds", project)))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    risk_df$risk_group <- ifelse(risk_df$risk_score >= stats::median(risk_df$risk_score, na.rm = TRUE), "High", "Low")
    if (requireNamespace("survival", quietly = TRUE)) {
      surv <- survival::survfit(survival::Surv(survival_time, survival_event) ~ risk_group + split, data = risk_df)
      p_km <- tryCatch(
        plot_km_curve(
          score_df = data.frame(sample_id = risk_df$sample_id, score = risk_df$risk_score, stringsAsFactors = FALSE),
          clinical = data.frame(
            sample_id = risk_df$sample_id,
            vital_status = ifelse(risk_df$survival_event == 1, "dead", "alive"),
            days_to_death = ifelse(risk_df$survival_event == 1, risk_df$survival_time, NA),
            days_to_last_follow_up = ifelse(risk_df$survival_event == 0, risk_df$survival_time, NA),
            stringsAsFactors = FALSE
          ),
          title = sprintf("Risk-score KM in %s", project)
        ),
        error = function(e) NULL
      )
      if (!is.null(p_km)) {
        save_plot(p_km, sprintf("%s_prognostic_risk_km", project), ctx$plots_dir, width = 7, height = 6)
      }
    }
  }

  write_report(
    ctx,
    sprintf("Prognostic Model for %s", project),
    c(
      method_summary,
      sprintf("- Train samples: `%s` (events=%s)", nrow(train), sum(train$survival_event, na.rm = TRUE)),
      sprintf("- Test samples: `%s` (events=%s)", nrow(test), sum(test$survival_event, na.rm = TRUE)),
      sprintf("- C-index train/test: `%s / %s`", signif(cindex_train, 3), signif(cindex_test, 3)),
      "",
      "## Outputs",
      sprintf("- `results/%s_prognostic_coefficients.csv`", project),
      sprintf("- `results/%s_prognostic_risk_scores.csv`", project),
      sprintf("- `results/%s_prognostic_performance.csv`", project),
      sprintf("- `objects/%s_prognostic_model.rds`", project)
    )
  )

  list(
    coefficients = coef_df,
    risk_scores = risk_df,
    performance = perf_df,
    method = method
  )
}
