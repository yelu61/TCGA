#!/usr/bin/env Rscript

# Post-processing script: compute pan-cancer pooled HR (95% CI) for each gene
# using fixed-effect meta-analysis across TCGA projects.
# Usage: Rscript enhance_fa_survival.R <run_dir>

args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]

if (length(args) < 1) {
  stop("Usage: Rscript enhance_fa_survival.R <run_dir>", call. = FALSE)
}

suppressPackageStartupMessages({
  library(survival)
})

script_flag <- grep("^--file=", commandArgs(), value = TRUE)[1]
if (!is.na(script_flag)) {
  script_dir <- normalizePath(dirname(sub("^--file=", "", script_flag)), mustWork = TRUE)
} else {
  script_dir <- normalizePath("tcga_toolkit/scripts", mustWork = TRUE)
}
options(tcga_toolkit.root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
source(file.path(script_dir, "common.R"))

results_dir <- file.path(run_dir, "results")
plots_dir <- file.path(run_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Read gene list from the correlation long file
long_df <- utils::read.csv(file.path(results_dir, "MYBL2_gene_correlation_long.csv"), stringsAsFactors = FALSE)
genes <- unique(long_df$gene)
projects <- available_projects(include_target = FALSE)

hr_rows <- list()
per_project_rows <- list()

for (gene in genes) {
  info("Processing survival for gene: %s", gene)
  project_hrs <- list()

  for (project in projects) {
    data <- tryCatch(standardize_project_data(project), error = function(e) NULL)
    if (is.null(data)) next

    feature_id <- resolve_feature_id(gene, data$tpm, data$gene_map)
    if (is.null(feature_id)) {
      next
    }

    expr <- data.frame(
      sample_id = colnames(data$tpm),
      score = log2(as.numeric(data$tpm[feature_id, ]) + 1),
      stringsAsFactors = FALSE
    )

    merged <- merge(expr, data$clinical, by = "sample_id", all.x = TRUE)
    cols <- survival_columns(merged)
    if (all(is.na(unlist(cols)))) next

    time <- suppressWarnings(as.numeric(merged[[cols$time_death]]))
    if (!is.na(cols$time_followup) && cols$time_followup %in% colnames(merged) && !all(is.na(merged[[cols$time_followup]]))) {
      followup <- suppressWarnings(as.numeric(merged[[cols$time_followup]]))
      time <- ifelse(is.na(time), followup, time)
    }
    status <- coerce_status(merged[[cols$status]])
    ok <- !is.na(time) & !is.na(status) & !is.na(merged$score) & is.finite(merged$score)
    if (sum(ok) < 10) next

    dat <- merged[ok, , drop = FALSE]
    dat$survival_time <- time[ok]
    dat$survival_event <- status[ok]
    cutoff <- stats::median(dat$score, na.rm = TRUE)
    dat$expr_group <- ifelse(dat$score >= cutoff, "High", "Low")
    dat$expr_group <- factor(dat$expr_group, levels = c("Low", "High"))
    if (length(unique(dat$expr_group)) < 2 || any(table(dat$expr_group) < 3)) next
    if (sum(dat$survival_event == 1, na.rm = TRUE) < 5) next

    fit <- tryCatch(
      survival::coxph(survival::Surv(time = survival_time, event = survival_event) ~ expr_group, data = dat),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    cox <- summary(fit)

    if (nrow(cox$conf.int) > 0 && is.finite(cox$conf.int[1, 1]) && is.finite(cox$conf.int[1, 3]) && is.finite(cox$conf.int[1, 4])) {
      project_hrs[[length(project_hrs) + 1]] <- data.frame(
        gene = gene,
        project = project,
        hr = unname(cox$conf.int[1, 1]),
        conf_low = unname(cox$conf.int[1, 3]),
        conf_high = unname(cox$conf.int[1, 4]),
        p = unname(cox$coefficients[1, 5]),
        n = sum(ok),
        events = sum(dat$survival_event == 1),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(project_hrs) >= 3) {
    df <- do.call(rbind, project_hrs)
    per_project_rows[[length(per_project_rows) + 1]] <- df
    log_hr <- log(df$hr)
    se <- (log(df$conf_high) - log(df$conf_low)) / (2 * 1.96)
    valid <- is.finite(log_hr) & is.finite(se) & se > 0 & df$conf_low > 0 & df$conf_high > 0

    if (sum(valid) >= 3) {
      log_hr <- log_hr[valid]
      se <- se[valid]
      weight <- 1 / se^2
      log_hr_pooled <- sum(weight * log_hr) / sum(weight)
      se_pooled <- sqrt(1 / sum(weight))
      hr_pooled <- exp(log_hr_pooled)
      ci_low <- exp(log_hr_pooled - 1.96 * se_pooled)
      ci_high <- exp(log_hr_pooled + 1.96 * se_pooled)
      z <- log_hr_pooled / se_pooled
      p_pooled <- 2 * (1 - stats::pnorm(abs(z)))

      Q <- sum(weight * (log_hr - log_hr_pooled)^2)
      df_q <- sum(valid) - 1
      i_squared <- max(0, (Q - df_q) / Q) * 100

      hr_rows[[length(hr_rows) + 1]] <- data.frame(
        gene = gene,
        hr = hr_pooled,
        ci_low = ci_low,
        ci_high = ci_high,
        p = p_pooled,
        n_projects = sum(valid),
        n_total = sum(df$n[valid]),
        events_total = sum(df$events[valid]),
        i_squared = i_squared,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(hr_rows) > 0) {
  result <- do.call(rbind, hr_rows)
  result <- result[order(-result$hr), ]
  out_path <- file.path(results_dir, "MYBL2_fa_hr_meta.csv")
  utils::write.csv(result, out_path, row.names = FALSE)

  if (length(per_project_rows) > 0) {
    per_project <- do.call(rbind, per_project_rows)
    utils::write.csv(per_project, file.path(results_dir, "MYBL2_fa_hr_per_project.csv"), row.names = FALSE)
  }

  message(sprintf("HR meta-analysis complete. Genes with pooled HR: %s", nrow(result)))
  message(sprintf("Top 5 HR genes: %s", paste(head(result$gene, 5), collapse = ", ")))
  message(sprintf("Output: %s", out_path))
} else {
  message("No genes had valid pooled HR after meta-analysis.")
}
