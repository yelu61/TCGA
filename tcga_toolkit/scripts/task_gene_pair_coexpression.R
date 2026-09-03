task_gene_pair_coexpression <- function(config, ctx) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    fail("gene_pair_coexpression requires the survival package.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    fail("gene_pair_coexpression requires the ggplot2 package.")
  }

  genes <- as.character(unlist(config$genes %||% character()))
  genes <- genes[nzchar(genes)]
  if (length(genes) != 2L) {
    fail("gene_pair_coexpression requires exactly two genes in config$genes.")
  }

  projects <- if (is.null(config$projects)) {
    available_projects(include_target = FALSE)
  } else {
    as.character(unlist(config$projects))
  }
  projects <- projects[nzchar(projects)]
  if (!length(projects)) {
    fail("gene_pair_coexpression requires at least one project.")
  }

  cutpoint <- config$cutpoint %||% "median"
  if (!identical(cutpoint, "median")) {
    fail("gene_pair_coexpression currently supports cutpoint = 'median'.")
  }

  sample_filter <- as.character(unlist(config$sample_filter %||% c("Tumor")))
  covariates <- as.character(unlist(config$covariates %||% c("age_at_index", "gender", "ajcc_pathologic_stage")))
  min_survival_n <- as.integer(config$min_survival_n %||% 30L)
  min_events <- as.integer(config$min_events %||% 10L)
  min_group_n <- as.integer(config$min_group_n %||% 5L)
  min_expression_tumor_n <- as.integer(config$min_expression_tumor_n %||% 10L)
  min_expression_normal_n <- as.integer(config$min_expression_normal_n %||% 5L)
  expression_p_adjust_method <- as.character(config$expression_p_adjust_method %||% "BH")
  sensitivity_thresholds <- as.character(unlist(config$sensitivity_thresholds %||% c("z_gt_0", "z_ge_2")))
  correlation_method <- as.character(config$correlation_method %||% "spearman")
  gtex_augment_normal <- isTRUE(config$gtex_augment_normal %||% FALSE)
  apply_cohort_filter <- isTRUE(config$apply_cohort_filter %||% FALSE)
  filter_min_tumor_n <- as.integer(config$filter_min_tumor_n %||% 100L)
  filter_min_normal_n <- as.integer(config$filter_min_normal_n %||% 10L)
  degenerate_cutoff_epsilon <- as.numeric(config$degenerate_cutoff_epsilon %||% 0.1)
  single_gene_binary_cox <- isTRUE(config$single_gene_binary_cox %||% TRUE)

  if (is.na(min_expression_tumor_n) || min_expression_tumor_n < 1L) {
    fail("gene_pair_coexpression requires min_expression_tumor_n >= 1.")
  }
  if (is.na(min_expression_normal_n) || min_expression_normal_n < 1L) {
    fail("gene_pair_coexpression requires min_expression_normal_n >= 1.")
  }
  tryCatch(
    stats::p.adjust(0.5, method = expression_p_adjust_method),
    error = function(e) fail("Invalid expression_p_adjust_method: %s", expression_p_adjust_method)
  )
  if (is.na(filter_min_tumor_n) || filter_min_tumor_n < 1L) {
    fail("gene_pair_coexpression requires filter_min_tumor_n >= 1.")
  }
  if (is.na(filter_min_normal_n) || filter_min_normal_n < 1L) {
    fail("gene_pair_coexpression requires filter_min_normal_n >= 1.")
  }
  if (is.na(degenerate_cutoff_epsilon) || degenerate_cutoff_epsilon < 0) {
    fail("gene_pair_coexpression requires degenerate_cutoff_epsilon >= 0.")
  }

  label <- paste(genes, collapse = "_")
  file_label <- gsub("[^A-Za-z0-9_-]+", "_", label)
  group_levels <- c("DoubleLow", paste0(genes[[1L]], "HighOnly"), paste0(genes[[2L]], "HighOnly"), "DoubleHigh")
  group_colors <- stats::setNames(c("#6B7280", "#2563EB", "#D97706", "#B91C1C"), group_levels)
  warnings_rows <- list()
  samples_rows <- list()
  expression_long_rows <- list()
  summary_rows <- list()
  group_count_rows <- list()
  correlation_rows <- list()
  survival_rows <- list()
  model_status_rows <- list()
  km_plot_rows <- list()
  audit_rows <- list()
  degenerate_flag_rows <- list()

  add_warning <- function(project, message) {
    warnings_rows <<- c(warnings_rows, list(data.frame(
      project = project,
      message = message,
      stringsAsFactors = FALSE
    )))
  }

  clean_missing_string <- function(x) {
    out <- as.character(x)
    bad <- is.na(out) | trimws(tolower(out)) %in% c(
      "", "na", "n/a", "not reported", "unknown", "not available",
      "not applicable", "[not available]", "[not applicable]", "--"
    )
    out[bad] <- NA_character_
    out
  }

  safe_z <- function(x, center = NULL, scale = NULL) {
    center <- center %||% mean(x, na.rm = TRUE)
    scale <- scale %||% stats::sd(x, na.rm = TRUE)
    if (is.na(scale) || scale == 0) {
      return(rep(NA_real_, length(x)))
    }
    (x - center) / scale
  }

  assign_group <- function(high1, high2) {
    out <- rep(NA_character_, length(high1))
    ok <- !is.na(high1) & !is.na(high2)
    out[ok & !high1 & !high2] <- "DoubleLow"
    out[ok & high1 & !high2] <- group_levels[[2L]]
    out[ok & !high1 & high2] <- group_levels[[3L]]
    out[ok & high1 & high2] <- "DoubleHigh"
    factor(out, levels = group_levels)
  }

  format_p_label <- function(p) {
    if (is.na(p)) {
      return("NA")
    }
    if (p < 0.001) {
      return(format(p, digits = 2, scientific = TRUE))
    }
    sprintf("%.3f", p)
  }

  p_to_stars <- function(p) {
    if (is.na(p)) {
      return("NA")
    }
    if (p < 0.001) {
      return("***")
    }
    if (p < 0.01) {
      return("**")
    }
    if (p < 0.05) {
      return("*")
    }
    "ns"
  }

  rank_biserial_effect <- function(tumor, normal) {
    tumor <- tumor[!is.na(tumor)]
    normal <- normal[!is.na(normal)]
    n_tumor <- length(tumor)
    n_normal <- length(normal)
    if (!n_tumor || !n_normal) {
      return(NA_real_)
    }
    ranks <- rank(c(tumor, normal), ties.method = "average")
    u_tumor <- sum(ranks[seq_len(n_tumor)]) - n_tumor * (n_tumor + 1) / 2
    2 * (u_tumor / (n_tumor * n_normal)) - 1
  }

  format_effect_value <- function(x) {
    if (is.na(x)) {
      return("NA")
    }
    if (abs(x) < 0.01) {
      x <- 0
    }
    sprintf("%.2f", x)
  }

  bind_rows_fill <- function(rows) {
    if (!length(rows)) {
      return(data.frame())
    }
    all_cols <- unique(unlist(lapply(rows, names)))
    rows <- lapply(rows, function(x) {
      missing_cols <- setdiff(all_cols, names(x))
      for (col in missing_cols) {
        x[[col]] <- NA
      }
      x[, all_cols, drop = FALSE]
    })
    do.call(rbind, rows)
  }

  gtex <- list(available = FALSE, long = data.frame(), tissues_by_project = list(), ensg = character())
  if (gtex_augment_normal) {
    assert_gtex_materialized <- function(path) {
      if (!file.exists(path)) {
        fail("GTEx file not found: %s", path)
      }
      listing <- tryCatch(
        system2("ls", c("-lO", shQuote(path)), stdout = TRUE, stderr = FALSE),
        error = function(e) character(),
        warning = function(w) character()
      )
      if (length(listing) && any(grepl("dataless", listing, fixed = TRUE))) {
        fail(
          "GTEx file '%s' is iCloud-evicted (dataless). Materialize it first, e.g.: brctl download '%s'  (or open it once in Finder), then re-run.",
          path,
          dirname(path)
        )
      }
      invisible(TRUE)
    }

    map_symbol_to_ensg <- function(symbol) {
      probemap <- file.path(project_root(), "GTEX", "gencode.v23.annotation.gene.probemap")
      if (file.exists(probemap)) {
        probe <- utils::read.delim(probemap, stringsAsFactors = FALSE, check.names = FALSE)
        gene_col <- intersect(c("gene", "gene_name", "symbol"), colnames(probe))
        id_col <- intersect(c("id", "gene_id"), colnames(probe))
        if (length(gene_col) && length(id_col)) {
          hit <- probe[[id_col[[1L]]]][probe[[gene_col[[1L]]]] == symbol]
          if (length(hit)) {
            return(hit[[1L]])
          }
        }
      }
      NA_character_
    }

    assert_gtex_materialized(file.path(project_root(), "GTEX", "gtex_RSEM_gene_tpm.gz"))
    assert_gtex_materialized(file.path(project_root(), "GTEX", "tissue_gtex.csv"))

    gtex_samples <- read_gtex_samples()
    gtex_long_rows <- list()
    ensg_map <- character()
    for (gene in genes) {
      ensg <- map_symbol_to_ensg(gene)
      if (is.na(ensg)) {
        add_warning("GTEx", sprintf("no ENSG mapping found for %s; GTEx augmentation skipped for this gene", gene))
        next
      }
      ensg_map[[gene]] <- ensg
      vec <- extract_gtex_gene_vector(ensg)
      if (is.null(vec)) {
        add_warning("GTEx", sprintf("gene %s (%s) was not found in the GTEx matrix", gene, ensg))
        next
      }
      gtex_long_rows <- c(gtex_long_rows, list(data.frame(
        sample_id = names(vec),
        gene = gene,
        expression = log2(pmax(as.numeric(vec), 0) + 1),
        stringsAsFactors = FALSE
      )))
    }
    if (length(gtex_long_rows)) {
      gtex_long <- bind_rows_fill(gtex_long_rows)
      gtex_long <- merge(gtex_long, gtex_samples[, c("sample_id", "tissue"), drop = FALSE], by = "sample_id", all.x = TRUE)
      gtex_long <- gtex_long[!is.na(gtex_long$tissue), , drop = FALSE]
      gtex <- list(
        available = TRUE,
        long = gtex_long,
        tissues_by_project = stats::setNames(lapply(projects, resolve_gtex_tissues), projects),
        ensg = ensg_map
      )
      info("GTEx augmentation enabled: %d genes extracted, %d samples with tissue annotation", length(gtex_long_rows), length(unique(gtex_long$sample_id)))
    } else {
      add_warning("GTEx", "no GTEx expression could be extracted for the requested genes; running TCGA-only")
    }
  }

  extract_survival_frame <- function(df) {
    cols <- survival_columns(df)
    if (is.na(cols$status) || !(cols$status %in% colnames(df))) {
      return(NULL)
    }

    time <- rep(NA_real_, nrow(df))
    if (!is.na(cols$time_death) && cols$time_death %in% colnames(df)) {
      time <- suppressWarnings(as.numeric(df[[cols$time_death]]))
    }
    if (!is.na(cols$time_followup) && cols$time_followup %in% colnames(df) && !all(is.na(df[[cols$time_followup]]))) {
      followup <- suppressWarnings(as.numeric(df[[cols$time_followup]]))
      time <- ifelse(is.na(time), followup, time)
    }
    status <- coerce_status(df[[cols$status]])
    ok <- !is.na(time) & time >= 0 & !is.na(status)
    if (!any(ok)) {
      return(NULL)
    }
    out <- df[ok, , drop = FALSE]
    out$survival_time <- time[ok]
    out$survival_event <- status[ok]
    out
  }

  prepare_adjustment <- function(df) {
    out <- df
    terms <- character()
    labels <- character()

    usable_factor <- function(values) {
      values <- droplevels(factor(values))
      ok <- !is.na(values) & !is.na(df$survival_event)
      if (sum(ok) < min_survival_n) {
        return(FALSE)
      }
      values <- droplevels(values[ok])
      if (nlevels(values) < 2L) {
        return(FALSE)
      }
      counts <- table(values)
      events <- tapply(df$survival_event[ok] == 1, values, sum)
      all(counts >= min_group_n) && all(events >= 1)
    }

    for (cov in covariates) {
      if (!(cov %in% colnames(out))) {
        next
      }

      if (cov == "ajcc_pathologic_stage") {
        term <- "cov_stage"
        values <- merge_clinical_stage(out[[cov]])
        if (usable_factor(values)) {
          out[[term]] <- values
          terms <- c(terms, term)
          labels <- c(labels, "ajcc_pathologic_stage_merged")
        }
        next
      }

      term <- paste0("cov_", make.names(cov))
      raw <- out[[cov]]
      numeric_value <- suppressWarnings(as.numeric(as.character(raw)))
      non_missing <- sum(!is.na(numeric_value))
      if (non_missing >= min_survival_n && stats::sd(numeric_value, na.rm = TRUE) > 0) {
        out[[term]] <- numeric_value
        terms <- c(terms, term)
        labels <- c(labels, cov)
      } else {
        values <- factor(clean_missing_string(raw))
        if (usable_factor(values)) {
          out[[term]] <- values
          terms <- c(terms, term)
          labels <- c(labels, cov)
        }
      }
    }

    list(data = out, terms = terms, labels = labels)
  }

  record_model_status <- function(project, project_code, model, feature, adjusted, status, reason, n, events, covariates_used) {
    model_status_rows <<- c(model_status_rows, list(data.frame(
      project = project,
      project_code = project_code,
      model = model,
      feature = feature,
      adjusted = adjusted,
      status = status,
      reason = reason,
      n = n,
      events = events,
      covariates = paste(covariates_used, collapse = ";"),
      stringsAsFactors = FALSE
    )))
  }

  run_cox <- function(df, rhs, required_vars, project, project_code, model, feature, adjusted = FALSE, covariates_used = character()) {
    needed <- unique(c("survival_time", "survival_event", required_vars))
    missing_vars <- setdiff(needed, colnames(df))
    if (length(missing_vars)) {
      record_model_status(project, project_code, model, feature, adjusted, "skipped", paste("missing variables:", paste(missing_vars, collapse = ",")), 0L, 0L, covariates_used)
      return(NULL)
    }

    dat <- df[, needed, drop = FALSE]
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    dat <- droplevels(dat)
    n <- nrow(dat)
    events <- sum(dat$survival_event == 1, na.rm = TRUE)
    if (n < min_survival_n) {
      record_model_status(project, project_code, model, feature, adjusted, "skipped", "insufficient complete survival samples", n, events, covariates_used)
      return(NULL)
    }
    if (events < min_events) {
      record_model_status(project, project_code, model, feature, adjusted, "skipped", "insufficient events", n, events, covariates_used)
      return(NULL)
    }
    for (var in required_vars) {
      if (is.factor(dat[[var]]) && nlevels(dat[[var]]) < 2L) {
        record_model_status(project, project_code, model, feature, adjusted, "skipped", sprintf("factor %s has fewer than two levels", var), n, events, covariates_used)
        return(NULL)
      }
      if (is.numeric(dat[[var]]) && (is.na(stats::sd(dat[[var]], na.rm = TRUE)) || stats::sd(dat[[var]], na.rm = TRUE) == 0)) {
        record_model_status(project, project_code, model, feature, adjusted, "skipped", sprintf("numeric %s has zero variance", var), n, events, covariates_used)
        return(NULL)
      }
    }

    formula <- stats::as.formula(sprintf("survival::Surv(survival_time, survival_event) ~ %s", rhs))
    fit <- tryCatch(
      survival::coxph(formula, data = dat),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      record_model_status(project, project_code, model, feature, adjusted, "failed", conditionMessage(fit), n, events, covariates_used)
      return(NULL)
    }

    cox <- summary(fit)
    if (!nrow(cox$coefficients)) {
      record_model_status(project, project_code, model, feature, adjusted, "skipped", "no coefficients returned", n, events, covariates_used)
      return(NULL)
    }
    record_model_status(project, project_code, model, feature, adjusted, "ok", "", n, events, covariates_used)
    data.frame(
      project = project,
      project_code = project_code,
      model = model,
      feature = feature,
      adjusted = adjusted,
      term = rownames(cox$coefficients),
      hazard_ratio = unname(cox$coefficients[, "exp(coef)"]),
      conf_low = if ("lower .95" %in% colnames(cox$conf.int)) unname(cox$conf.int[, "lower .95"]) else NA_real_,
      conf_high = if ("upper .95" %in% colnames(cox$conf.int)) unname(cox$conf.int[, "upper .95"]) else NA_real_,
      p_value = unname(cox$coefficients[, "Pr(>|z|)"]),
      n = n,
      events = events,
      covariates = paste(covariates_used, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  plot_km_group <- function(df, group_col, title, basename, palette = NULL) {
    if (!requireNamespace("survminer", quietly = TRUE)) {
      return(NULL)
    }
    if (!(group_col %in% colnames(df))) {
      return(NULL)
    }
    dat <- df[!is.na(df[[group_col]]) & !is.na(df$survival_time) & !is.na(df$survival_event), , drop = FALSE]
    dat$surv_group <- droplevels(factor(dat[[group_col]]))
    if (nrow(dat) < 10L || length(unique(dat$surv_group)) < 2L) {
      return(NULL)
    }
    group_sizes <- table(dat$surv_group)
    if (any(group_sizes < min_group_n)) {
      return(NULL)
    }
    fit <- survival::survfit(survival::Surv(survival_time, survival_event) ~ surv_group, data = dat)
    palette_use <- palette
    if (!is.null(palette_use)) {
      palette_use <- unname(palette_use[levels(dat$surv_group)])
    }
    p <- survminer::ggsurvplot(
      fit,
      data = dat,
      pval = TRUE,
      risk.table = TRUE,
      title = title,
      xlab = "Time (days)",
      ylab = "Overall survival probability",
      legend.title = group_col,
      legend.labs = levels(dat$surv_group),
      palette = palette_use,
      risk.table.y.text.col = TRUE,
      risk.table.y.text = FALSE
    )
    p$plot <- p$plot + theme_tcga(base_size = 10) + ggplot2::theme(legend.position = "bottom")
    saved <- save_plot(p, basename, ctx$plots_dir, width = 7, height = 6)
    data.frame(
      plot_type = group_col,
      pdf = saved$pdf,
      png = saved$png,
      stringsAsFactors = FALSE
    )
  }

  plot_single_gene_forest <- function(df, title) {
    df <- df[
      !is.na(df$hazard_ratio) & df$hazard_ratio > 0 &
        !is.na(df$conf_low) & df$conf_low > 0 & !is.na(df$conf_high),
      ,
      drop = FALSE
    ]
    if (!nrow(df)) {
      return(NULL)
    }
    df$log2_hr <- log2(df$hazard_ratio)
    df$log2_conf_low <- log2(df$conf_low)
    df$log2_conf_high <- log2(df$conf_high)
    unadj <- df[!df$adjusted, , drop = FALSE]
    order_levels <- unadj$project_code[order(unadj$log2_hr, na.last = TRUE)]
    order_levels <- unique(c(order_levels, df$project_code))
    df$project_code <- factor(df$project_code, levels = order_levels)
    df$model_label <- factor(ifelse(df$adjusted, "Adjusted", "Unadjusted"), levels = c("Unadjusted", "Adjusted"))
    ggplot2::ggplot(df, ggplot2::aes(x = log2_hr, y = project_code, color = model_label, shape = model_label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = log2_conf_low, xmax = log2_conf_high),
        width = 0.18,
        orientation = "y",
        position = ggplot2::position_dodge(width = 0.5)
      ) +
      ggplot2::geom_point(size = 2.4, position = ggplot2::position_dodge(width = 0.5)) +
      ggplot2::scale_color_manual(values = c(Unadjusted = "#64748B", Adjusted = "#B91C1C")) +
      ggplot2::scale_shape_manual(values = c(Unadjusted = 16, Adjusted = 17)) +
      ggplot2::labs(title = title, x = expression(log[2] ~ HR), y = NULL, color = "Model", shape = "Model") +
      theme_tcga(base_size = 10) +
      ggplot2::theme(legend.position = "bottom")
  }

  for (project in projects) {
    project_code <- tcga_project_code(project)
    info("Processing gene pair coexpression: %s", project)

    result <- tryCatch({
      data <- standardize_project_data(project)
      gene1 <- resolve_expression_vector(data, genes[[1L]])
      gene2 <- resolve_expression_vector(data, genes[[2L]])
      gene1 <- gene1[, c("sample_id", "value"), drop = FALSE]
      gene2 <- gene2[, c("sample_id", "value"), drop = FALSE]
      names(gene1)[2L] <- "gene1_expression"
      names(gene2)[2L] <- "gene2_expression"
      merged <- merge(gene1, gene2, by = "sample_id")
      merged <- merge(merged, data$clinical, by = "sample_id", all.x = TRUE)
      merged$project <- project
      merged$project_code <- project_code
      merged$dataset_group <- ifelse(as.character(merged$Group) == "Normal", "Normal", ifelse(as.character(merged$Group) %in% sample_filter, "Tumor", "Other"))

      analysis <- merged[
        as.character(merged$Group) %in% sample_filter &
          !is.na(merged$gene1_expression) &
          !is.na(merged$gene2_expression),
        ,
        drop = FALSE
      ]
      if (nrow(analysis) < 10L) {
        add_warning(project, "skipped: fewer than 10 analysis samples with both genes")
        next
      }

      n_tumor <- sum(merged$dataset_group == "Tumor", na.rm = TRUE)
      n_tcga_normal <- sum(merged$dataset_group == "Normal", na.rm = TRUE)
      gtex_tissues <- if (gtex$available && project %in% names(gtex$tissues_by_project)) gtex$tissues_by_project[[project]] else character()
      if (is.null(gtex_tissues)) {
        gtex_tissues <- character()
      }
      n_gtex_normal <- 0L
      if (gtex$available && length(gtex_tissues)) {
        n_gtex_normal <- length(unique(gtex$long$sample_id[gtex$long$tissue %in% gtex_tissues]))
      }
      n_normal_total <- n_tcga_normal + n_gtex_normal
      filter_reasons <- character()
      if (n_tumor < filter_min_tumor_n) {
        filter_reasons <- c(filter_reasons, sprintf("tumor n %d < %d", n_tumor, filter_min_tumor_n))
      }
      if (n_normal_total < filter_min_normal_n) {
        filter_reasons <- c(filter_reasons, sprintf("normal total %d < %d", n_normal_total, filter_min_normal_n))
      }
      audit_rows <- c(audit_rows, list(data.frame(
        project = project,
        project_code = project_code,
        n_tumor = n_tumor,
        n_tcga_normal = n_tcga_normal,
        n_gtex_normal = n_gtex_normal,
        n_normal_total = n_normal_total,
        gtex_tissues = paste(gtex_tissues, collapse = ";"),
        included = if (apply_cohort_filter) !length(filter_reasons) else TRUE,
        reason = if (apply_cohort_filter) paste(filter_reasons, collapse = "; ") else "",
        stringsAsFactors = FALSE
      )))
      if (apply_cohort_filter && length(filter_reasons)) {
        add_warning(project, sprintf("excluded by cohort filter: %s", paste(filter_reasons, collapse = "; ")))
        next
      }

      gene1_cutoff <- stats::median(analysis$gene1_expression, na.rm = TRUE)
      gene2_cutoff <- stats::median(analysis$gene2_expression, na.rm = TRUE)
      gene1_mean <- mean(analysis$gene1_expression, na.rm = TRUE)
      gene2_mean <- mean(analysis$gene2_expression, na.rm = TRUE)
      gene1_sd <- stats::sd(analysis$gene1_expression, na.rm = TRUE)
      gene2_sd <- stats::sd(analysis$gene2_expression, na.rm = TRUE)

      gene1_cutoff_near_zero <- is.finite(gene1_cutoff) && gene1_cutoff <= degenerate_cutoff_epsilon
      gene2_cutoff_near_zero <- is.finite(gene2_cutoff) && gene2_cutoff <= degenerate_cutoff_epsilon
      degenerate_split <- gene1_cutoff_near_zero || gene2_cutoff_near_zero
      if (degenerate_split) {
        degenerate_flag_rows <- c(degenerate_flag_rows, list(data.frame(
          project = project,
          project_code = project_code,
          gene = genes,
          median_cutoff = c(gene1_cutoff, gene2_cutoff),
          near_zero = c(gene1_cutoff_near_zero, gene2_cutoff_near_zero),
          stringsAsFactors = FALSE
        )))
        add_warning(project, sprintf(
          "degenerate median split (%s cutoff=%.3f, %s cutoff=%.3f, epsilon=%.2f): coexpression groups may collapse",
          genes[[1L]], gene1_cutoff, genes[[2L]], gene2_cutoff, degenerate_cutoff_epsilon
        ))
      }

      merged$gene1 <- genes[[1L]]
      merged$gene2 <- genes[[2L]]
      merged$gene1_z_tumor_ref <- safe_z(merged$gene1_expression, gene1_mean, gene1_sd)
      merged$gene2_z_tumor_ref <- safe_z(merged$gene2_expression, gene2_mean, gene2_sd)
      merged$gene1_high_median <- merged$gene1_expression >= gene1_cutoff
      merged$gene2_high_median <- merged$gene2_expression >= gene2_cutoff
      merged$coexpression_group_median <- assign_group(merged$gene1_high_median, merged$gene2_high_median)
      merged$gene1_high_z_gt_0 <- merged$gene1_z_tumor_ref > 0
      merged$gene2_high_z_gt_0 <- merged$gene2_z_tumor_ref > 0
      merged$coexpression_group_z_gt_0 <- assign_group(merged$gene1_high_z_gt_0, merged$gene2_high_z_gt_0)
      merged$gene1_high_z_ge_2 <- merged$gene1_z_tumor_ref >= 2
      merged$gene2_high_z_ge_2 <- merged$gene2_z_tumor_ref >= 2
      merged$coexpression_group_z_ge_2 <- assign_group(merged$gene1_high_z_ge_2, merged$gene2_high_z_ge_2)

      analysis <- merged[merged$sample_id %in% analysis$sample_id, , drop = FALSE]
      analysis$gene1_z <- safe_z(analysis$gene1_expression)
      analysis$gene2_z <- safe_z(analysis$gene2_expression)
      analysis$gene1_group <- factor(ifelse(analysis$gene1_high_median, "High", "Low"), levels = c("Low", "High"))
      analysis$gene2_group <- factor(ifelse(analysis$gene2_high_median, "High", "Low"), levels = c("Low", "High"))
      analysis$coexpression_group <- factor(as.character(analysis$coexpression_group_median), levels = group_levels)

      samples_rows <- c(samples_rows, list(merged))

      expression_long_rows <- c(expression_long_rows, list(data.frame(
        project = project,
        project_code = project_code,
        sample_id = merged$sample_id,
        dataset_group = merged$dataset_group,
        gene = genes[[1L]],
        expression = merged$gene1_expression,
        source = "TCGA",
        stringsAsFactors = FALSE
      )))
      expression_long_rows <- c(expression_long_rows, list(data.frame(
        project = project,
        project_code = project_code,
        sample_id = merged$sample_id,
        dataset_group = merged$dataset_group,
        gene = genes[[2L]],
        expression = merged$gene2_expression,
        source = "TCGA",
        stringsAsFactors = FALSE
      )))
      if (gtex$available && length(gtex_tissues)) {
        gtex_sub <- gtex$long[gtex$long$tissue %in% gtex_tissues & !is.na(gtex$long$expression), , drop = FALSE]
        if (nrow(gtex_sub)) {
          expression_long_rows <- c(expression_long_rows, list(data.frame(
            project = project,
            project_code = project_code,
            sample_id = gtex_sub$sample_id,
            dataset_group = "Normal",
            gene = gtex_sub$gene,
            expression = gtex_sub$expression,
            source = "GTEx",
            stringsAsFactors = FALSE
          )))
        }
      }

      threshold_groups <- c("median", sensitivity_thresholds)
      threshold_groups <- unique(threshold_groups[threshold_groups %in% c("median", "z_gt_0", "z_ge_2")])
      for (threshold in threshold_groups) {
        group_col <- paste0("coexpression_group_", threshold)
        group_values <- factor(as.character(analysis[[group_col]]), levels = group_levels)
        counts <- as.integer(table(group_values))
        total <- sum(counts)
        group_count_rows <- c(group_count_rows, list(data.frame(
          project = project,
          project_code = project_code,
          threshold = threshold,
          group = group_levels,
          n = counts,
          proportion = if (total > 0) counts / total else NA_real_,
          stringsAsFactors = FALSE
        )))
      }

      cor_result <- data.frame(
        project = project,
        project_code = project_code,
        method = correlation_method,
        n = nrow(analysis),
        correlation = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      if (nrow(analysis) >= 10L && stats::sd(analysis$gene1_expression) > 0 && stats::sd(analysis$gene2_expression) > 0) {
        ct <- suppressWarnings(stats::cor.test(analysis$gene1_expression, analysis$gene2_expression, method = correlation_method))
        cor_result$correlation <- unname(ct$estimate)
        cor_result$p_value <- ct$p.value
      }
      correlation_rows <- c(correlation_rows, list(cor_result))

      surv_dat <- extract_survival_frame(analysis)
      n_survival <- 0L
      events <- 0L
      if (!is.null(surv_dat)) {
        n_survival <- nrow(surv_dat)
        events <- sum(surv_dat$survival_event == 1, na.rm = TRUE)
        adjusted_data <- prepare_adjustment(surv_dat)

        gene1_rows <- run_cox(surv_dat, "gene1_z", "gene1_z", project, project_code, "single_gene_continuous", genes[[1L]])
        gene2_rows <- run_cox(surv_dat, "gene2_z", "gene2_z", project, project_code, "single_gene_continuous", genes[[2L]])
        pair_rows <- run_cox(surv_dat, "gene1_z * gene2_z", c("gene1_z", "gene2_z"), project, project_code, "gene_pair_interaction", paste(genes, collapse = "+"))
        group_rows <- run_cox(surv_dat, "coexpression_group", "coexpression_group", project, project_code, "coexpression_group", paste(genes, collapse = "+"))
        gene1_bin_rows <- NULL
        gene2_bin_rows <- NULL
        if (single_gene_binary_cox) {
          gene1_bin_rows <- run_cox(surv_dat, "gene1_group", "gene1_group", project, project_code, "single_gene_binary", genes[[1L]])
          gene2_bin_rows <- run_cox(surv_dat, "gene2_group", "gene2_group", project, project_code, "single_gene_binary", genes[[2L]])
        }
        for (cox_rows in list(gene1_rows, gene2_rows, pair_rows, group_rows, gene1_bin_rows, gene2_bin_rows)) {
          if (!is.null(cox_rows)) survival_rows <- c(survival_rows, list(cox_rows))
        }

        if (length(adjusted_data$terms)) {
          adj_rhs <- paste(c("gene1_z", adjusted_data$terms), collapse = " + ")
          gene1_adj <- run_cox(adjusted_data$data, adj_rhs, c("gene1_z", adjusted_data$terms), project, project_code, "single_gene_continuous", genes[[1L]], adjusted = TRUE, covariates_used = adjusted_data$labels)
          adj_rhs <- paste(c("gene2_z", adjusted_data$terms), collapse = " + ")
          gene2_adj <- run_cox(adjusted_data$data, adj_rhs, c("gene2_z", adjusted_data$terms), project, project_code, "single_gene_continuous", genes[[2L]], adjusted = TRUE, covariates_used = adjusted_data$labels)
          adj_rhs <- paste(c("gene1_z * gene2_z", adjusted_data$terms), collapse = " + ")
          pair_adj <- run_cox(adjusted_data$data, adj_rhs, c("gene1_z", "gene2_z", adjusted_data$terms), project, project_code, "gene_pair_interaction", paste(genes, collapse = "+"), adjusted = TRUE, covariates_used = adjusted_data$labels)
          adj_rhs <- paste(c("coexpression_group", adjusted_data$terms), collapse = " + ")
          group_adj <- run_cox(adjusted_data$data, adj_rhs, c("coexpression_group", adjusted_data$terms), project, project_code, "coexpression_group", paste(genes, collapse = "+"), adjusted = TRUE, covariates_used = adjusted_data$labels)
          gene1_bin_adj <- NULL
          gene2_bin_adj <- NULL
          if (single_gene_binary_cox) {
            adj_rhs <- paste(c("gene1_group", adjusted_data$terms), collapse = " + ")
            gene1_bin_adj <- run_cox(adjusted_data$data, adj_rhs, c("gene1_group", adjusted_data$terms), project, project_code, "single_gene_binary", genes[[1L]], adjusted = TRUE, covariates_used = adjusted_data$labels)
            adj_rhs <- paste(c("gene2_group", adjusted_data$terms), collapse = " + ")
            gene2_bin_adj <- run_cox(adjusted_data$data, adj_rhs, c("gene2_group", adjusted_data$terms), project, project_code, "single_gene_binary", genes[[2L]], adjusted = TRUE, covariates_used = adjusted_data$labels)
          }
          for (cox_rows in list(gene1_adj, gene2_adj, pair_adj, group_adj, gene1_bin_adj, gene2_bin_adj)) {
            if (!is.null(cox_rows)) survival_rows <- c(survival_rows, list(cox_rows))
          }
        } else {
          record_model_status(project, project_code, "adjusted_models", paste(genes, collapse = "+"), TRUE, "skipped", "no usable adjustment covariates", n_survival, events, character())
        }

        km1 <- plot_km_group(
          surv_dat,
          "gene1_group",
          sprintf("%s %s high/low OS", project_code, genes[[1L]]),
          sprintf("%s_%s_high_low_km", project, genes[[1L]]),
          palette = c(Low = "#4069B0", High = "#D1495B")
        )
        km2 <- plot_km_group(
          surv_dat,
          "gene2_group",
          sprintf("%s %s high/low OS", project_code, genes[[2L]]),
          sprintf("%s_%s_high_low_km", project, genes[[2L]]),
          palette = c(Low = "#4069B0", High = "#D1495B")
        )
        km_pair <- plot_km_group(
          surv_dat,
          "coexpression_group",
          sprintf("%s %s coexpression OS", project_code, label),
          sprintf("%s_%s_coexpression_group_km", project, file_label),
          palette = group_colors
        )
        for (km in list(km1, km2, km_pair)) {
          if (!is.null(km)) {
            km$project <- project
            km$project_code <- project_code
            km_plot_rows <- c(km_plot_rows, list(km))
          }
        }
      } else {
        add_warning(project, "no usable OS time/event data after filtering")
      }

      median_groups <- table(factor(as.character(analysis$coexpression_group), levels = group_levels))
      z0_groups <- table(factor(as.character(analysis$coexpression_group_z_gt_0), levels = group_levels))
      z2_groups <- table(factor(as.character(analysis$coexpression_group_z_ge_2), levels = group_levels))
      summary_rows <- c(summary_rows, list(data.frame(
        project = project,
        project_code = project_code,
        n_samples_all = nrow(merged),
        n_analysis = nrow(analysis),
        n_survival = n_survival,
        events = events,
        gene1 = genes[[1L]],
        gene2 = genes[[2L]],
        gene1_median_cutoff = gene1_cutoff,
        gene2_median_cutoff = gene2_cutoff,
        gene1_cutoff_near_zero = gene1_cutoff_near_zero,
        gene2_cutoff_near_zero = gene2_cutoff_near_zero,
        degenerate_split = degenerate_split,
        median_double_high_n = unname(median_groups[["DoubleHigh"]]),
        median_double_high_prop = unname(median_groups[["DoubleHigh"]]) / nrow(analysis),
        median_concordant_prop = unname(median_groups[["DoubleHigh"]] + median_groups[["DoubleLow"]]) / nrow(analysis),
        z_gt_0_double_high_prop = unname(z0_groups[["DoubleHigh"]]) / nrow(analysis),
        z_ge_2_double_high_prop = unname(z2_groups[["DoubleHigh"]]) / nrow(analysis),
        spearman_correlation = cor_result$correlation,
        spearman_p_value = cor_result$p_value,
        stringsAsFactors = FALSE
      )))

      TRUE
    }, error = function(e) {
      add_warning(project, conditionMessage(e))
      audit_rows <<- c(audit_rows, list(data.frame(
        project = project,
        project_code = project_code,
        n_tumor = NA_integer_,
        n_tcga_normal = NA_integer_,
        n_gtex_normal = NA_integer_,
        n_normal_total = NA_integer_,
        gtex_tissues = "",
        included = FALSE,
        reason = paste("project data load failed:", conditionMessage(e)),
        stringsAsFactors = FALSE
      )))
      NULL
    })

    if (is.null(result)) {
      next
    }
  }

  if (!length(samples_rows)) {
    details <- if (length(warnings_rows)) {
      paste(vapply(warnings_rows, function(x) paste(x$project, x$message, sep = ": "), character(1)), collapse = " | ")
    } else {
      "no project-level warning was captured"
    }
    fail("No gene pair coexpression results were generated. %s", details)
  }

  samples_df <- bind_rows_fill(samples_rows)
  expression_long <- bind_rows_fill(expression_long_rows)
  project_summary <- bind_rows_fill(summary_rows)
  expression_compare <- expression_long[
    expression_long$dataset_group %in% c("Tumor", "Normal") & !is.na(expression_long$expression),
    ,
    drop = FALSE
  ]
  expression_group_summary <- if (nrow(expression_compare)) {
    do.call(rbind, lapply(split(
      expression_compare,
      interaction(expression_compare$project, expression_compare$gene, expression_compare$dataset_group, drop = TRUE)
    ), function(x) {
      data.frame(
        project = x$project[[1L]],
        project_code = x$project_code[[1L]],
        gene = x$gene[[1L]],
        dataset_group = x$dataset_group[[1L]],
        n = nrow(x),
        mean_expression = mean(x$expression, na.rm = TRUE),
        median_expression = stats::median(x$expression, na.rm = TRUE),
        q1_expression = unname(stats::quantile(x$expression, 0.25, na.rm = TRUE)),
        q3_expression = unname(stats::quantile(x$expression, 0.75, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame()
  }
  expression_tn_stats <- if (nrow(expression_long)) {
    expression_keys <- unique(expression_long[, c("project", "project_code", "gene"), drop = FALSE])
    do.call(rbind, lapply(seq_len(nrow(expression_keys)), function(i) {
      key <- expression_keys[i, , drop = FALSE]
      dat <- expression_compare[
        expression_compare$project == key$project & expression_compare$gene == key$gene,
        ,
        drop = FALSE
      ]
      tumor <- dat$expression[dat$dataset_group == "Tumor"]
      normal <- dat$expression[dat$dataset_group == "Normal"]
      tumor_n <- length(tumor)
      normal_n <- length(normal)
      normal_n_tcga <- sum(dat$dataset_group == "Normal" & dat$source == "TCGA", na.rm = TRUE)
      normal_n_gtex <- sum(dat$dataset_group == "Normal" & dat$source == "GTEx", na.rm = TRUE)
      normal_tcga <- dat$expression[dat$dataset_group == "Normal" & dat$source == "TCGA"]
      median_tumor <- if (tumor_n) stats::median(tumor, na.rm = TRUE) else NA_real_
      median_normal <- if (normal_n) stats::median(normal, na.rm = TRUE) else NA_real_
      reasons <- character()
      if (tumor_n < min_expression_tumor_n) {
        reasons <- c(reasons, sprintf("tumour n < %s", min_expression_tumor_n))
      }
      if (normal_n < min_expression_normal_n) {
        reasons <- c(reasons, sprintf("normal n < %s", min_expression_normal_n))
      }
      status <- if (length(reasons)) "skipped_low_n" else "tested"
      reason <- paste(reasons, collapse = "; ")
      p_value <- NA_real_
      effect <- NA_real_
      if (identical(status, "tested")) {
        wt <- tryCatch(
          stats::wilcox.test(tumor, normal, exact = FALSE),
          error = function(e) e
        )
        if (inherits(wt, "error")) {
          status <- "failed"
          reason <- conditionMessage(wt)
        } else {
          p_value <- wt$p.value
          effect <- rank_biserial_effect(tumor, normal)
        }
      }
      p_value_tcga_only <- NA_real_
      effect_tcga_only <- NA_real_
      if (tumor_n >= min_expression_tumor_n && normal_n_tcga >= min_expression_normal_n) {
        wt_tcga <- tryCatch(
          stats::wilcox.test(tumor, normal_tcga, exact = FALSE),
          error = function(e) e
        )
        if (!inherits(wt_tcga, "error")) {
          p_value_tcga_only <- wt_tcga$p.value
          effect_tcga_only <- rank_biserial_effect(tumor, normal_tcga)
        }
      }
      data.frame(
        project = key$project,
        project_code = key$project_code,
        gene = key$gene,
        tumor_n = tumor_n,
        normal_n = normal_n,
        normal_n_tcga = normal_n_tcga,
        normal_n_gtex = normal_n_gtex,
        status = status,
        reason = reason,
        median_tumor = median_tumor,
        median_normal = median_normal,
        median_diff_tumor_minus_normal = median_tumor - median_normal,
        rank_biserial_tumor_vs_normal = effect,
        p_value = p_value,
        rank_biserial_tcga_only = effect_tcga_only,
        p_value_tcga_only = p_value_tcga_only,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame()
  }
  if (nrow(expression_tn_stats)) {
    expression_tn_stats$fdr <- NA_real_
    expression_tn_stats$fdr_tcga_only <- NA_real_
    for (gene in unique(expression_tn_stats$gene)) {
      idx <- which(expression_tn_stats$gene == gene & expression_tn_stats$status == "tested" & !is.na(expression_tn_stats$p_value))
      if (length(idx)) {
        expression_tn_stats$fdr[idx] <- stats::p.adjust(expression_tn_stats$p_value[idx], method = expression_p_adjust_method)
      }
      idx_tcga <- which(expression_tn_stats$gene == gene & !is.na(expression_tn_stats$p_value_tcga_only))
      if (length(idx_tcga)) {
        expression_tn_stats$fdr_tcga_only[idx_tcga] <- stats::p.adjust(expression_tn_stats$p_value_tcga_only[idx_tcga], method = expression_p_adjust_method)
      }
    }
    expression_tn_stats$significance <- vapply(expression_tn_stats$fdr, p_to_stars, character(1))
  }
  group_counts <- bind_rows_fill(group_count_rows)
  correlation_df <- bind_rows_fill(correlation_rows)
  if (nrow(correlation_df) && any(!is.na(correlation_df$p_value))) {
    correlation_df$fdr <- NA_real_
    idx <- which(!is.na(correlation_df$p_value))
    correlation_df$fdr[idx] <- stats::p.adjust(correlation_df$p_value[idx], method = "BH")
  }

  survival_df <- bind_rows_fill(survival_rows)
  if (nrow(survival_df)) {
    survival_df$fdr <- NA_real_
    keys <- paste(survival_df$model, survival_df$feature, survival_df$adjusted, survival_df$term, sep = "|")
    for (key in unique(keys)) {
      idx <- which(keys == key & !is.na(survival_df$p_value))
      if (length(idx)) {
        survival_df$fdr[idx] <- stats::p.adjust(survival_df$p_value[idx], method = "BH")
      }
    }
  }

  model_status <- bind_rows_fill(model_status_rows)
  warnings_df <- if (length(warnings_rows)) bind_rows_fill(warnings_rows) else data.frame(project = character(), message = character())
  km_plots <- bind_rows_fill(km_plot_rows)
  cohort_audit <- bind_rows_fill(audit_rows)
  degenerate_flags <- bind_rows_fill(degenerate_flag_rows)
  degenerate_codes <- character()
  if (nrow(project_summary) && "degenerate_split" %in% colnames(project_summary)) {
    degenerate_codes <- unique(project_summary$project_code[!is.na(project_summary$degenerate_split) & project_summary$degenerate_split])
  }

  double_high_term <- "coexpression_groupDoubleHigh"
  double_high <- if (nrow(survival_df)) {
    survival_df[survival_df$model == "coexpression_group" & survival_df$term == double_high_term, , drop = FALSE]
  } else {
    data.frame()
  }
  if (nrow(double_high)) {
    double_high <- do.call(rbind, lapply(split(double_high, double_high$project), function(x) {
      adjusted <- x[x$adjusted, , drop = FALSE]
      if (nrow(adjusted)) {
        return(adjusted[1L, , drop = FALSE])
      }
      x[1L, , drop = FALSE]
    }))
  }

  write_table_safe(samples_df, file.path(ctx$results_dir, sprintf("%s_sample_expression.csv", file_label)))
  write_table_safe(expression_long, file.path(ctx$results_dir, sprintf("%s_expression_long.csv", file_label)))
  write_table_safe(expression_group_summary, file.path(ctx$results_dir, sprintf("%s_expression_group_summary.csv", file_label)))
  write_table_safe(expression_tn_stats, file.path(ctx$results_dir, sprintf("%s_tumor_normal_expression_stats.csv", file_label)))
  write_table_safe(project_summary, file.path(ctx$results_dir, sprintf("%s_project_summary.csv", file_label)))
  write_table_safe(group_counts, file.path(ctx$results_dir, sprintf("%s_coexpression_group_counts.csv", file_label)))
  write_table_safe(correlation_df, file.path(ctx$results_dir, sprintf("%s_correlation_summary.csv", file_label)))
  write_table_safe(survival_df, file.path(ctx$results_dir, sprintf("%s_survival_models.csv", file_label)))
  write_table_safe(double_high, file.path(ctx$results_dir, sprintf("%s_double_high_survival_summary.csv", file_label)))
  write_table_safe(model_status, file.path(ctx$results_dir, sprintf("%s_model_status.csv", file_label)))
  write_table_safe(km_plots, file.path(ctx$results_dir, sprintf("%s_km_plots.csv", file_label)))
  write_table_safe(warnings_df, file.path(ctx$results_dir, sprintf("%s_warnings.csv", file_label)))
  write_table_safe(cohort_audit, file.path(ctx$results_dir, sprintf("%s_cohort_filter_audit.csv", file_label)))
  write_table_safe(degenerate_flags, file.path(ctx$results_dir, sprintf("%s_degenerate_group_flags.csv", file_label)))

  saveRDS(
    list(
      samples = samples_df,
      expression_long = expression_long,
      expression_group_summary = expression_group_summary,
      tumor_normal_expression = expression_tn_stats,
      project_summary = project_summary,
      group_counts = group_counts,
      correlation = correlation_df,
      survival = survival_df,
      model_status = model_status,
      cohort_filter_audit = cohort_audit,
      degenerate_flags = degenerate_flags,
      warnings = warnings_df
    ),
    file.path(ctx$objects_dir, sprintf("%s_gene_pair_coexpression.rds", file_label))
  )

  write_json(
    list(
      backend = "tcga-toolkit",
      toolkit_version = ctx$toolkit_version,
      task = ctx$task,
      task_id = ctx$task_id,
      genes = as.list(genes),
      projects = as.list(projects),
      input_scale = "log2(TPM+1)",
      sample_filter = as.list(sample_filter),
      cutpoint = cutpoint,
      expression_tumor_normal_filter = list(
        min_tumor_n = min_expression_tumor_n,
        min_normal_n = min_expression_normal_n,
        p_adjust_method = expression_p_adjust_method
      ),
      sensitivity_thresholds = as.list(sensitivity_thresholds),
      survival_endpoint = "overall_survival",
      survival_time_rule = "days_to_death; otherwise days_to_last_follow_up",
      cox_adjustment_defaults = as.list(covariates),
      single_gene_binary_cox = single_gene_binary_cox,
      gtex_augment_normal = gtex_augment_normal,
      gtex_genes_ensg = as.list(gtex$ensg),
      gtex_tissues_by_project = if (gtex$available) gtex$tissues_by_project else list(),
      cohort_filter = list(
        enabled = apply_cohort_filter,
        min_tumor_n = filter_min_tumor_n,
        min_normal_n = filter_min_normal_n,
        n_projects_included = if (nrow(cohort_audit)) sum(cohort_audit$included %in% TRUE) else 0L,
        n_projects_excluded = if (nrow(cohort_audit)) sum(!is.na(cohort_audit$included) & !cohort_audit$included) else 0L
      ),
      degenerate_cutoff_epsilon = degenerate_cutoff_epsilon,
      degenerate_projects = as.list(degenerate_codes),
      cross_study_batch_caveat = if (gtex$available) {
        "TCGA tumour TPM (GDC STAR-Counts) and GTEx normal TPM (RSEM, GTEx v8) come from different quantification pipelines and are NOT batch-corrected. Tumour-vs-augmented-normal differences partly reflect this technical offset; TCGA-only normal p-values (columns p_value_tcga_only, fdr_tcga_only in the tumor_normal_expression_stats table) are provided as a sensitivity analysis. Interpret GTEx-augmented comparisons with caution."
      } else {
        NULL
      },
      cbiportal_note = "z_gt_0 and z_ge_2 are tumour-reference z-score sensitivity analyses from local TPM data; they do not use cBioPortal diploid CNA reference populations."
    ),
    file.path(ctx$results_dir, sprintf("%s_analysis_manifest.json", file_label))
  )

  plot_expression <- expression_long[expression_long$dataset_group %in% c("Tumor", "Normal"), , drop = FALSE]
  if (nrow(plot_expression)) {
    p_expr <- ggplot2::ggplot(
      plot_expression,
      ggplot2::aes(x = project_code, y = expression, fill = gene)
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge(width = 0.8), alpha = 0.8) +
      ggplot2::facet_wrap(~dataset_group, ncol = 1, scales = "free_y") +
      ggplot2::labs(
        title = sprintf("%s paired pan-cancer expression", label),
        x = "TCGA cancer type",
        y = "log2(TPM + 1)",
        caption = if (gtex$available) "Normal panel combines TCGA adjacent normal and GTEx matched-tissue samples (not batch-corrected)." else NULL
      ) +
      ggplot2::scale_fill_manual(values = stats::setNames(c("#2563EB", "#D97706"), genes)) +
      theme_tcga(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save_plot(p_expr, sprintf("%s_paired_expression_boxplot", file_label), ctx$plots_dir, width = max(12, length(unique(plot_expression$project_code)) * 0.45), height = 8)
  }

  for (gene in genes) {
    gene_stats <- expression_tn_stats[
      expression_tn_stats$gene == gene & expression_tn_stats$status == "tested",
      ,
      drop = FALSE
    ]
    gene_plot <- plot_expression[
      plot_expression$gene == gene & plot_expression$project %in% gene_stats$project,
      ,
      drop = FALSE
    ]
    if (!nrow(gene_stats) || !nrow(gene_plot)) {
      next
    }
    project_levels <- sort(unique(gene_stats$project_code))
    gene_plot$project_code <- factor(gene_plot$project_code, levels = project_levels)
    gene_plot$dataset_group <- factor(gene_plot$dataset_group, levels = c("Normal", "Tumor"))
    gene_stats <- gene_stats[match(project_levels, gene_stats$project_code), , drop = FALSE]
    gene_stats$project_code <- factor(gene_stats$project_code, levels = project_levels)
    y_positions <- stats::aggregate(expression ~ project_code, data = gene_plot, FUN = max)
    y_range <- range(gene_plot$expression, na.rm = TRUE)
    y_pad <- if (all(is.finite(y_range)) && diff(y_range) > 0) diff(y_range) * 0.08 else 0.5
    gene_stats$label_y <- y_positions$expression[match(as.character(gene_stats$project_code), as.character(y_positions$project_code))] + y_pad
    gene_stats$label <- sprintf(
      "%s\nr=%s",
      gene_stats$significance,
      vapply(gene_stats$rank_biserial_tumor_vs_normal, format_effect_value, character(1))
    )
    gene_file <- gsub("[^A-Za-z0-9_-]+", "_", gene)
    p_gene <- ggplot2::ggplot(
      gene_plot,
      ggplot2::aes(x = project_code, y = expression, fill = dataset_group)
    ) +
      ggplot2::geom_boxplot(
        outlier.shape = NA,
        position = ggplot2::position_dodge(width = 0.75),
        width = 0.65,
        alpha = 0.82
      )
    if (gtex$available) {
      p_gene <- p_gene +
        ggplot2::geom_jitter(
          ggplot2::aes(color = dataset_group, alpha = source),
          position = ggplot2::position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75),
          size = 0.45
        ) +
        ggplot2::scale_alpha_manual(values = c(TCGA = 0.28, GTEx = 0.1), name = "Normal\nsource")
    } else {
      p_gene <- p_gene +
        ggplot2::geom_jitter(
          ggplot2::aes(color = dataset_group),
          position = ggplot2::position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75),
          alpha = 0.22,
          size = 0.45,
          show.legend = FALSE
        )
    }
    caption_text <- sprintf(
      "Wilcoxon Tumor vs Normal; stars use %s-adjusted FDR. r > 0 means Tumor > Normal.",
      expression_p_adjust_method
    )
    if (gtex$available) {
      caption_text <- sprintf(
        paste0(
          "Wilcoxon Tumor vs combined Normal (TCGA adjacent + GTEx matched tissue, faint points); stars use %s-adjusted FDR. r > 0 means Tumor > Normal.\n",
          "TCGA-only normal sensitivity p-values: see results/%s_tumor_normal_expression_stats.csv (p_value_tcga_only, fdr_tcga_only).\n",
          "TCGA (GDC STAR-TPM) and GTEx (RSEM-TPM) use different quantification pipelines and are not batch-corrected; interpret augmented comparisons with caution."
        ),
        expression_p_adjust_method,
        file_label
      )
    }
    p_gene <- p_gene +
      ggplot2::geom_text(
        data = gene_stats,
        ggplot2::aes(x = project_code, y = label_y, label = label),
        inherit.aes = FALSE,
        size = 2.5,
        lineheight = 0.85
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.16))) +
      ggplot2::scale_fill_manual(values = c(Normal = "#2563EB", Tumor = "#B91C1C"), drop = FALSE) +
      ggplot2::scale_color_manual(values = c(Normal = "#2563EB", Tumor = "#B91C1C"), drop = FALSE, guide = "none") +
      ggplot2::labs(
        title = sprintf("%s tumour-normal expression by cancer type", gene),
        x = "TCGA cancer type",
        y = "log2(TPM + 1)",
        fill = "Sample type",
        caption = caption_text
      ) +
      theme_tcga(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.caption = ggplot2::element_text(hjust = 0, size = 8)
      )
    save_plot(
      p_gene,
      sprintf("%s_tumor_normal_boxplot", gene_file),
      ctx$plots_dir,
      width = max(10, length(project_levels) * 0.45),
      height = 6.2
    )
  }

  scatter_df <- samples_df[samples_df$dataset_group == "Tumor", , drop = FALSE]
  if (nrow(scatter_df)) {
    annot <- correlation_df
    annot$label <- sprintf("rho=%s\nP=%s\nn=%s", sprintf("%.2f", annot$correlation), vapply(annot$p_value, format_p_label, character(1)), annot$n)
    scatter_df$project_label <- ifelse(scatter_df$project_code %in% degenerate_codes, paste0(scatter_df$project_code, " *"), scatter_df$project_code)
    annot$project_label <- ifelse(annot$project_code %in% degenerate_codes, paste0(annot$project_code, " *"), annot$project_code)
    p_scatter <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = gene1_expression, y = gene2_expression)) +
      ggplot2::geom_point(alpha = 0.45, size = 0.8, color = "#334155") +
      ggplot2::facet_wrap(~project_label, scales = "free", ncol = 6) +
      ggplot2::geom_text(
        data = annot,
        ggplot2::aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.15,
        size = 2.4
      ) +
      ggplot2::labs(title = sprintf("%s tumour coexpression", label), x = sprintf("%s log2(TPM + 1)", genes[[1L]]), y = sprintf("%s log2(TPM + 1)", genes[[2L]])) +
      theme_tcga(base_size = 9)
    if (length(degenerate_codes)) {
      p_scatter <- p_scatter + ggplot2::labs(caption = "* median cutoff near 0 for at least one gene; coexpression groups degenerate (one gene effectively expressed/absent).")
    }
    save_plot(p_scatter, sprintf("%s_tumor_scatter_facets", file_label), ctx$plots_dir, width = 12, height = max(8, ceiling(length(unique(scatter_df$project_label)) / 6) * 2.2))
  }

  main_counts <- group_counts[group_counts$threshold == "median", , drop = FALSE]
  if (nrow(main_counts)) {
    main_counts$project_label <- ifelse(main_counts$project_code %in% degenerate_codes, paste0(main_counts$project_code, " *"), main_counts$project_code)
    p_groups <- ggplot2::ggplot(main_counts, ggplot2::aes(x = project_label, y = proportion, fill = group)) +
      ggplot2::geom_col(width = 0.85) +
      ggplot2::scale_y_continuous(labels = function(x) sprintf("%d%%", round(100 * x))) +
      ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
      ggplot2::labs(title = sprintf("%s median-split coexpression groups", label), x = "TCGA cancer type", y = "Tumour sample proportion", fill = "Group") +
      theme_tcga(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
    if (length(degenerate_codes)) {
      p_groups <- p_groups + ggplot2::labs(caption = "* median cutoff near 0 for at least one gene; coexpression groups degenerate (one gene effectively expressed/absent).")
    }
    save_plot(p_groups, sprintf("%s_coexpression_group_proportions", file_label), ctx$plots_dir, width = max(11, length(unique(main_counts$project_label)) * 0.35), height = 5.5)
  }

  if (nrow(double_high)) {
    double_high$log2_hr <- log2(double_high$hazard_ratio)
    double_high$log2_conf_low <- log2(double_high$conf_low)
    double_high$log2_conf_high <- log2(double_high$conf_high)
    double_high$project_code <- factor(double_high$project_code, levels = double_high$project_code[order(double_high$log2_hr, na.last = TRUE)])
    p_forest <- ggplot2::ggplot(double_high, ggplot2::aes(x = log2_hr, y = project_code)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = log2_conf_low, xmax = log2_conf_high), width = 0.18, orientation = "y", color = "#334155") +
      ggplot2::geom_point(ggplot2::aes(color = adjusted), size = 2.6) +
      ggplot2::scale_color_manual(values = c("FALSE" = "#64748B", "TRUE" = "#B91C1C"), labels = c("Unadjusted", "Adjusted")) +
      ggplot2::labs(title = "Double-high vs double-low OS risk by cancer type", x = expression(log[2] ~ HR), y = NULL, color = "Model") +
      theme_tcga(base_size = 10) +
      ggplot2::theme(legend.position = "bottom")
    save_plot(p_forest, sprintf("%s_double_high_survival_forest", file_label), ctx$plots_dir, width = 7, height = max(6, nrow(double_high) * 0.25))

    bubble_df <- double_high[!is.na(double_high$fdr) & !is.na(double_high$hazard_ratio), , drop = FALSE]
    if (nrow(bubble_df)) {
      bubble_df$neg_log10_fdr <- -log10(pmax(bubble_df$fdr, .Machine$double.xmin))
      p_bubble <- ggplot2::ggplot(bubble_df, ggplot2::aes(x = log2(hazard_ratio), y = neg_log10_fdr, size = events, color = adjusted)) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::geom_text(ggplot2::aes(label = project_code), size = 2.6, vjust = -0.7, check_overlap = TRUE) +
        ggplot2::scale_color_manual(values = c("FALSE" = "#64748B", "TRUE" = "#B91C1C"), labels = c("Unadjusted", "Adjusted")) +
        ggplot2::labs(title = "Double-high OS risk significance", x = expression(log[2] ~ HR), y = expression(-log[10] ~ FDR), size = "Events", color = "Model") +
        theme_tcga(base_size = 10) +
        ggplot2::theme(legend.position = "bottom")
      save_plot(p_bubble, sprintf("%s_double_high_risk_bubble", file_label), ctx$plots_dir, width = 8, height = 6)
    }
  }

  if (nrow(survival_df)) {
    for (i in seq_along(genes)) {
      gene <- genes[[i]]
      gene_file <- gsub("[^A-Za-z0-9_-]+", "_", gene)
      forest_specs <- list(
        list(model = "single_gene_continuous", term = sprintf("gene%d_z", i), tag = "continuous", desc = "continuous z-score, per-SD HR"),
        list(model = "single_gene_binary", term = sprintf("gene%d_groupHigh", i), tag = "binary", desc = "median-split High vs Low")
      )
      for (spec in forest_specs) {
        sub <- survival_df[
          survival_df$model == spec$model & survival_df$feature == gene & survival_df$term == spec$term,
          ,
          drop = FALSE
        ]
        if (!nrow(sub)) {
          next
        }
        p_sf <- plot_single_gene_forest(sub, sprintf("%s OS by cancer type (%s)", gene, spec$desc))
        if (!is.null(p_sf)) {
          save_plot(
            p_sf,
            sprintf("%s_single_gene_%s_survival_forest", gene_file, spec$tag),
            ctx$plots_dir,
            width = 7.5,
            height = max(6, length(unique(sub$project_code)) * 0.32)
          )
        }
      }
    }
  }

  if (nrow(survival_df)) {
    primary_terms <- data.frame(
      model = c("single_gene_continuous", "single_gene_continuous", "gene_pair_interaction", "coexpression_group"),
      feature = c(genes[[1L]], genes[[2L]], paste(genes, collapse = "+"), paste(genes, collapse = "+")),
      term = c("gene1_z", "gene2_z", "gene1_z:gene2_z", double_high_term),
      comparison = c(genes[[1L]], genes[[2L]], "Interaction", "DoubleHigh"),
      stringsAsFactors = FALSE
    )
    heatmap_df <- merge(
      survival_df[!survival_df$adjusted, , drop = FALSE],
      primary_terms,
      by = c("model", "feature", "term")
    )
    heatmap_df <- heatmap_df[!is.na(heatmap_df$fdr) & !is.na(heatmap_df$hazard_ratio), , drop = FALSE]
    if (nrow(heatmap_df)) {
      heatmap_df$signed_neg_log10_fdr <- sign(log2(heatmap_df$hazard_ratio)) * pmin(-log10(pmax(heatmap_df$fdr, .Machine$double.xmin)), 6)
      p_heat <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = comparison, y = project_code, fill = signed_neg_log10_fdr)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.2) +
        ggplot2::scale_fill_gradient2(low = "#2563EB", mid = "white", high = "#B91C1C", midpoint = 0, name = "signed\n-log10 FDR") +
        ggplot2::labs(title = sprintf("%s OS association map", label), x = NULL, y = "TCGA cancer type") +
        theme_tcga(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
      save_plot(p_heat, sprintf("%s_survival_fdr_heatmap", file_label), ctx$plots_dir, width = 7, height = max(6, length(unique(heatmap_df$project_code)) * 0.25))
    }
  }

  report_bullets <- c(
    sprintf("- Genes: `%s` and `%s`", genes[[1L]], genes[[2L]]),
    sprintf("- Projects requested: `%s`", length(projects)),
    sprintf("- Projects with expression results: `%s`", length(unique(project_summary$project))),
    sprintf("- Input scale: `log2(TPM + 1)`"),
    sprintf("- Main grouping: tumour samples split within each cancer type by each gene's median expression."),
    sprintf(
      "- Tumour-normal expression tests require at least `%s` tumour and `%s` normal samples per cancer type; skipped projects remain in the stats table and are not removed from tumour-only coexpression or survival analyses.",
      min_expression_tumor_n,
      min_expression_normal_n
    ),
    sprintf("- Sensitivity grouping: `%s`", paste(sensitivity_thresholds, collapse = ", ")),
    sprintf("- Survival endpoint: overall survival, per cancer type only."),
    sprintf("- Adjusted Cox covariates attempted: `%s`", paste(covariates, collapse = ", ")),
    sprintf("- Single-gene binary Cox (median-split High vs Low, matching KM grouping): `%s`", single_gene_binary_cox),
    sprintf("- cBioPortal note: `z_ge_2` is a local tumour-reference z-score sensitivity analysis, not a diploid-CNA cBioPortal z-score.")
  )
  if (apply_cohort_filter) {
    report_bullets <- c(report_bullets, sprintf(
      "- Cohort filter ON: projects need tumour n >= `%s` and combined normal n >= `%s`; included `%s`, excluded `%s` (see `results/%s_cohort_filter_audit.csv`).",
      filter_min_tumor_n,
      filter_min_normal_n,
      sum(cohort_audit$included %in% TRUE),
      sum(!is.na(cohort_audit$included) & !cohort_audit$included),
      file_label
    ))
  } else {
    report_bullets <- c(report_bullets, "- Cohort filter OFF: all projects analysed; per-project counts in `results/` cohort filter audit table.")
  }
  if (gtex$available) {
    gtex_pairs <- vapply(
      names(gtex$tissues_by_project),
      function(p) if (length(gtex$tissues_by_project[[p]])) sprintf("%s=%s", tcga_project_code(p), paste(gtex$tissues_by_project[[p]], collapse = "+")) else NA_character_,
      character(1)
    )
    gtex_pairs <- gtex_pairs[!is.na(gtex_pairs)]
    report_bullets <- c(
      report_bullets,
      sprintf("- GTEx normal augmentation ON: Normal group = TCGA adjacent normal + GTEx matched tissue (%s).", paste(gtex_pairs, collapse = ", ")),
      "- **Cross-study batch caveat**: TCGA tumour TPM (GDC STAR-Counts) and GTEx normal TPM (RSEM, GTEx v8) use different quantification pipelines and are NOT batch-corrected. Tumour-vs-augmented-normal differences partly reflect this technical offset; TCGA-only sensitivity p-values are in `p_value_tcga_only` / `fdr_tcga_only` of the tumour-normal stats table. Interpret GTEx-augmented comparisons with caution."
    )
  } else {
    report_bullets <- c(report_bullets, "- GTEx normal augmentation OFF: Normal group = TCGA adjacent normal only.")
  }
  if (length(degenerate_codes)) {
    report_bullets <- c(report_bullets, sprintf(
      "- Degenerate median split (cutoff <= `%s`, at least one gene effectively expressed/absent) in: `%s`. Group proportions/KM/Cox for these cancer types should be interpreted with caution; see `results/%s_degenerate_group_flags.csv`. Marked with `*` on plots.",
      degenerate_cutoff_epsilon,
      paste(degenerate_codes, collapse = ", "),
      file_label
    ))
  }
  report_bullets <- c(
    report_bullets,
    "",
    "## Outputs",
    sprintf("- `results/%s_sample_expression.csv`", file_label),
    sprintf("- `results/%s_expression_long.csv`", file_label),
    sprintf("- `results/%s_expression_group_summary.csv`", file_label),
    sprintf("- `results/%s_tumor_normal_expression_stats.csv`", file_label),
    sprintf("- `results/%s_project_summary.csv`", file_label),
    sprintf("- `results/%s_coexpression_group_counts.csv`", file_label),
    sprintf("- `results/%s_correlation_summary.csv`", file_label),
    sprintf("- `results/%s_survival_models.csv`", file_label),
    sprintf("- `results/%s_double_high_survival_summary.csv`", file_label),
    sprintf("- `results/%s_cohort_filter_audit.csv`", file_label),
    sprintf("- `results/%s_degenerate_group_flags.csv`", file_label),
    sprintf("- `plots/%s_paired_expression_boxplot.pdf/png`", file_label),
    sprintf("- `plots/%s_tumor_normal_boxplot.pdf/png`", gsub("[^A-Za-z0-9_-]+", "_", genes[[1L]])),
    sprintf("- `plots/%s_tumor_normal_boxplot.pdf/png`", gsub("[^A-Za-z0-9_-]+", "_", genes[[2L]])),
    sprintf("- `plots/%s_tumor_scatter_facets.pdf/png`", file_label),
    sprintf("- `plots/%s_coexpression_group_proportions.pdf/png`", file_label),
    sprintf("- `plots/%s_double_high_survival_forest.pdf/png`", file_label),
    sprintf("- `plots/%s_survival_fdr_heatmap.pdf/png`", file_label),
    sprintf("- `plots/%s_double_high_risk_bubble.pdf/png`", file_label),
    sprintf("- `plots/%s_single_gene_continuous_survival_forest.pdf/png`", gsub("[^A-Za-z0-9_-]+", "_", genes[[1L]])),
    sprintf("- `plots/%s_single_gene_binary_survival_forest.pdf/png`", gsub("[^A-Za-z0-9_-]+", "_", genes[[1L]])),
    sprintf("- `plots/%s_single_gene_continuous_survival_forest.pdf/png`", gsub("[^A-Za-z0-9_-]+", "_", genes[[2L]])),
    sprintf("- `plots/%s_single_gene_binary_survival_forest.pdf/png`", gsub("[^A-Za-z0-9_-]+", "_", genes[[2L]]))
  )

  write_report(
    ctx,
    sprintf("Gene Pair Coexpression for %s", label),
    report_bullets
  )

  list(
    samples = samples_df,
    expression_group_summary = expression_group_summary,
    tumor_normal_expression = expression_tn_stats,
    project_summary = project_summary,
    group_counts = group_counts,
    correlation = correlation_df,
    survival = survival_df,
    double_high = double_high,
    cohort_filter_audit = cohort_audit,
    degenerate_flags = degenerate_flags,
    warnings = warnings_df
  )
}
