task_tmb_analysis <- function(config, ctx) {
  project <- config$project %||% fail("tmb_analysis requires a project.")
  tmb_source <- tolower(config$tmb_source %||% "auto")
  if (!(tmb_source %in% c("auto", "maf", "pancancer"))) {
    fail("tmb_analysis.tmb_source must be one of auto/maf/pancancer.")
  }
  exome_size_mb <- as.numeric(config$exome_size_mb %||% 38)

  data <- standardize_project_data(project)
  clinical <- data$clinical

  tmb_df <- NULL
  source_used <- NA_character_

  if (tmb_source %in% c("auto", "pancancer")) {
    pc_path <- file.path(project_root(), "0-Data", "PanCancer_TMB.Rdata")
    if (file.exists(pc_path)) {
      pc <- tryCatch(load_rdata_file(pc_path), error = function(e) NULL)
      candidates <- if (is.data.frame(pc)) list(pc) else if (is.list(pc)) Filter(is.data.frame, pc) else list()
      for (df in candidates) {
        if (any(grepl("TMB|nonsilent|mutation", colnames(df), ignore.case = TRUE))) {
          df <- as.data.frame(df, stringsAsFactors = FALSE)
          sid_col <- intersect(c("Tumor_Sample_Barcode", "sample_id", "Sample", "samplesID", "pan.samplesID", "patient"), colnames(df))[1]
          if (is.na(sid_col)) next
          tmb_col <- grep("TMB|nonsilent|mutation", colnames(df), ignore.case = TRUE, value = TRUE)[1]
          if (is.na(tmb_col)) next
          tmp <- data.frame(
            sample_id = substr(as.character(df[[sid_col]]), 1, 15),
            patient_short = substr(as.character(df[[sid_col]]), 1, 12),
            tmb = suppressWarnings(as.numeric(df[[tmb_col]])),
            stringsAsFactors = FALSE
          )
          tmp <- tmp[!is.na(tmp$tmb), , drop = FALSE]
          if (nrow(tmp)) {
            project_filter <- if ("cancer.type" %in% colnames(df)) {
              code <- tcga_project_code(project)
              keep <- as.character(df$cancer.type) %in% c(code, project)
              tmp <- tmp[keep, , drop = FALSE]
            } else tmp
            if (nrow(tmp)) {
              tmb_df <- tmp
              source_used <- "PanCancer_TMB.Rdata"
              break
            }
          }
        }
      }
    }
  }

  if (is.null(tmb_df) && tmb_source %in% c("auto", "maf")) {
    files <- project_files(project)
    if (!file.exists(files$maf)) {
      fail("No PanCancer TMB record and no MAF file at %s", files$maf)
    }
    if (!requireNamespace("maftools", quietly = TRUE)) {
      fail("Computing TMB from MAF requires the maftools package.")
    }
    maf_obj <- load_rdata_file(files$maf)
    if (!inherits(maf_obj, "MAF")) {
      maf_obj <- tryCatch(maftools::read.maf(maf_obj, verbose = FALSE), error = function(e) NULL)
    }
    if (is.null(maf_obj) || !inherits(maf_obj, "MAF")) {
      fail("Could not parse MAF for %s.", project)
    }
    tmb_obj <- maftools::tmb(maf = maf_obj, captureSize = exome_size_mb, logScale = FALSE)
    tmb_obj <- as.data.frame(tmb_obj, stringsAsFactors = FALSE)
    tmb_col <- intersect(c("total_perMB", "total"), colnames(tmb_obj))[1]
    sid_col <- intersect(c("Tumor_Sample_Barcode", "samples"), colnames(tmb_obj))[1]
    tmb_df <- data.frame(
      sample_id = substr(as.character(tmb_obj[[sid_col]]), 1, 15),
      patient_short = substr(as.character(tmb_obj[[sid_col]]), 1, 12),
      tmb = suppressWarnings(as.numeric(tmb_obj[[tmb_col]])),
      stringsAsFactors = FALSE
    )
    source_used <- "MAF"
  }

  if (is.null(tmb_df) || !nrow(tmb_df)) {
    fail("No TMB data could be assembled for %s.", project)
  }

  clinical$patient_short <- clinical$patient_short %||% substr(as.character(clinical$patient_barcode %||% clinical$sample_id), 1, 12)
  merged <- merge(tmb_df, clinical, by = "patient_short", all.x = TRUE, suffixes = c("", ".clin"))
  merged <- merged[!is.na(merged$tmb), , drop = FALSE]

  threshold <- if (!is.null(config$tmb_threshold)) as.numeric(config$tmb_threshold)
               else as.numeric(stats::quantile(merged$tmb, 0.75, na.rm = TRUE))
  merged$tmb_group <- ifelse(merged$tmb >= threshold, "TMB-High", "TMB-Low")

  summary_df <- data.frame(
    n_samples = nrow(merged),
    median_tmb = stats::median(merged$tmb, na.rm = TRUE),
    mean_tmb = mean(merged$tmb, na.rm = TRUE),
    threshold = threshold,
    n_high = sum(merged$tmb_group == "TMB-High"),
    n_low = sum(merged$tmb_group == "TMB-Low"),
    source = source_used,
    stringsAsFactors = FALSE
  )

  write_table_safe(merged[, c("sample_id", "patient_short", "tmb", "tmb_group"), drop = FALSE],
                   file.path(ctx$results_dir, sprintf("%s_tmb_per_sample.csv", project)))
  write_table_safe(summary_df, file.path(ctx$results_dir, sprintf("%s_tmb_summary.csv", project)))

  survival <- run_survival_models(
    score_df = data.frame(sample_id = merged$sample_id.clin %||% merged$sample_id, score = merged$tmb, stringsAsFactors = FALSE),
    clinical = data$clinical
  )
  if (!is.null(survival)) {
    write_table_safe(survival$cox, file.path(ctx$results_dir, sprintf("%s_tmb_survival_cox.csv", project)))
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_hist <- ggplot2::ggplot(merged, ggplot2::aes(x = log10(tmb + 0.01))) +
      ggplot2::geom_histogram(bins = 40, fill = "#4C72B0", colour = "white") +
      ggplot2::geom_vline(xintercept = log10(threshold + 0.01), linetype = "dashed", colour = "red") +
      ggplot2::labs(title = sprintf("%s TMB distribution (source=%s)", project, source_used),
                    x = "log10(mutations / Mb)", y = "Samples") +
      theme_tcga()
    save_plot(p_hist, sprintf("%s_tmb_histogram", project), ctx$plots_dir, width = 8, height = 5)

    if ("tmb_group" %in% colnames(merged) && !is.null(survival)) {
      p_km <- tryCatch(
        plot_km_curve(
          score_df = data.frame(sample_id = merged$sample_id.clin %||% merged$sample_id, score = merged$tmb, stringsAsFactors = FALSE),
          clinical = data$clinical,
          title = sprintf("Survival by TMB in %s", project)
        ),
        error = function(e) NULL
      )
      if (!is.null(p_km)) {
        save_plot(p_km, sprintf("%s_tmb_km", project), ctx$plots_dir, width = 7, height = 6)
      }
    }
  }

  write_report(
    ctx,
    sprintf("TMB analysis for %s", project),
    c(
      sprintf("- Source: `%s`", source_used),
      sprintf("- Samples with TMB: `%s`", nrow(merged)),
      sprintf("- Median TMB: `%s`", signif(summary_df$median_tmb, 4)),
      sprintf("- Threshold (mut/Mb): `%s`", signif(threshold, 4)),
      sprintf("- TMB-High / TMB-Low: `%s / %s`", summary_df$n_high, summary_df$n_low),
      sprintf("- Cox model: `%s`", !is.null(survival)),
      "",
      "## Outputs",
      sprintf("- `results/%s_tmb_per_sample.csv`", project),
      sprintf("- `results/%s_tmb_summary.csv`", project)
    )
  )

  list(per_sample = merged, summary = summary_df, survival = survival, source = source_used)
}
