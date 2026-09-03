task_maf_summary <- function(config, ctx) {
  if (!requireNamespace("maftools", quietly = TRUE)) {
    fail("task_maf_summary requires the maftools package.")
  }

  project <- config$project %||% fail("maf_summary requires a project.")
  files <- project_files(project)
  if (!file.exists(files$maf)) {
    fail("Missing MAF file for %s: %s", project, files$maf)
  }

  maf_df <- load_rdata_file(files$maf)
  if (!is.data.frame(maf_df)) {
    fail("MAF object for %s is not a data.frame.", project)
  }

  required_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")
  missing_cols <- setdiff(required_cols, colnames(maf_df))
  if (length(missing_cols)) {
    fail("MAF data for %s is missing required columns: %s", project, paste(missing_cols, collapse = ", "))
  }

  maf_df <- as.data.frame(maf_df)
  # top_n / lof_top_n accept a vector to emit one oncoplot per value.
  top_n <- as.integer(config$top_n %||% 20)

  # --- optional clinical annotation / LoF focus parameters -------------------
  clinical_columns <- config$clinical_columns %||% NULL
  if (length(clinical_columns)) {
    clinical_columns <- as.character(clinical_columns)
  }
  group_column <- config$group_column %||% (if (length(clinical_columns)) clinical_columns[[1L]] else NULL)
  lof_focus <- isTRUE(config$lof_focus %||% FALSE)
  lof_classes <- as.character(config$lof_classes %||%
    c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"))
  lof_top_n <- as.integer(config$lof_top_n %||% top_n)
  age_breaks <- as.numeric(config$age_breaks %||% c(50, 65))
  # How to order samples WITHIN each group_column block on the annotated
  # oncoplot: "tmb" sorts groups by descending median mutation burden and
  # samples by descending burden inside each group (clean waterfall look);
  # "annotation" (default) keeps maftools' sortByAnnotation ordering.
  within_group_order <- tolower(as.character(config$within_group_order %||% "annotation"))

  # Derived clinical columns built on top of the raw patient table:
  # stage_merged (Stage I-IV), grade_merged (G1-G4), age_group (binned by
  # age_breaks), molecular_subtype (from the local Pan-Cancer subtype table).
  add_derived_clinical <- function(patient) {
    if ("ajcc_pathologic_tumor_stage" %in% colnames(patient)) {
      patient$stage_merged <- merge_clinical_stage(patient$ajcc_pathologic_tumor_stage)
    }
    if ("tumor_grade" %in% colnames(patient)) {
      patient$grade_merged <- merge_clinical_grade(patient$tumor_grade)
    }
    if ("age_at_initial_pathologic_diagnosis" %in% colnames(patient)) {
      patient$age_group <- bin_age_group(patient$age_at_initial_pathologic_diagnosis, age_breaks)
    }
    subtypes <- get_project_subtypes(project)
    if (!is.null(subtypes)) {
      patient$molecular_subtype <- subtypes$molecular_subtype[
        match(patient$bcr_patient_barcode, subtypes$bcr_patient_barcode)]
    }
    patient
  }

  # Short display labels for oncoplot side annotation rows.
  display_label <- function(cols) {
    known <- c(
      molecular_subtype = "Subtype", stage_merged = "Stage",
      grade_merged = "Grade", age_group = "Age", gender = "Gender",
      ajcc_pathologic_tumor_stage = "Stage", tumor_grade = "Grade",
      age_at_initial_pathologic_diagnosis = "Age",
      h_pylori_infection = "H.pylori", histologic_diagnosis = "Histology",
      ajcc_tumor_pathologic_pt = "T", ajcc_nodes_pathologic_pn = "N",
      ajcc_metastasis_pathologic_pm = "M", vital_status = "Vital"
    )
    out <- cols
    hit <- cols %in% names(known)
    out[hit] <- known[cols[hit]]
    out
  }

  # Build a Tumor_Sample_Barcode-keyed clinical annotation table for maftools.
  # Returns NULL when no clinical_columns were requested.
  build_clinical_annotation <- function(maf_data) {
    if (!length(clinical_columns)) {
      return(NULL)
    }
    patient <- get_patient_clinical(project)
    if (!("bcr_patient_barcode" %in% colnames(patient))) {
      info("Clinical patient table for %s lacks bcr_patient_barcode; skipping annotation.", project)
      return(NULL)
    }
    # Drop non-barcode header artefacts that GDC TSVs sometimes carry.
    patient <- patient[grepl("^TCGA-", patient$bcr_patient_barcode), , drop = FALSE]
    patient <- add_derived_clinical(patient)
    avail <- intersect(clinical_columns, colnames(patient))
    missing_clin <- setdiff(clinical_columns, colnames(patient))
    if (length(missing_clin)) {
      info("Clinical columns not present for %s and skipped: %s", project, paste(missing_clin, collapse = ", "))
    }
    if (!length(avail)) {
      return(NULL)
    }
    samples <- sort(unique(maf_data$Tumor_Sample_Barcode))
    # read.maf(isTCGA = TRUE) truncates barcodes to the 12-char patient code,
    # so the clinical annotation must be keyed by that same patient barcode.
    patient_barcodes <- sort(unique(substr(samples, 1L, 12L)))
    matched <- patient[match(patient_barcodes, patient$bcr_patient_barcode), avail, drop = FALSE]
    # Convert bracketed GDC placeholders to NA so they do not pollute the legend.
    matched[] <- lapply(matched, function(col) {
      if (is.factor(col)) {
        return(col)
      }
      col <- as.character(col)
      col[grepl("^\\[", col)] <- NA_character_
      col
    })
    ann <- data.frame(Tumor_Sample_Barcode = patient_barcodes, matched,
                      stringsAsFactors = FALSE, check.names = FALSE)
    # Put the grouping column first so sortByAnnotation orders by it.
    if (!is.null(group_column) && group_column %in% avail) {
      ann <- ann[, c("Tumor_Sample_Barcode", group_column, setdiff(avail, group_column)), drop = FALSE]
    }
    colnames(ann) <- c("Tumor_Sample_Barcode", display_label(colnames(ann)[-1L]))
    rownames(ann) <- NULL
    ann
  }

  # Compute an explicit sample order for the annotated oncoplot: groups
  # ordered by descending median TMB, samples by descending TMB within each
  # group. `clinical_ann` is the display-labelled annotation table; the TMB
  # per 12-char patient barcode comes from the MAF. Returns a character
  # vector of barcodes, or NULL to fall back to sortByAnnotation.
  compute_tmb_sample_order <- function(maf_data, clinical_ann) {
    if (is.null(clinical_ann) || is.null(group_column)) {
      return(NULL)
    }
    group_label <- display_label(group_column)
    if (!(group_label %in% colnames(clinical_ann))) {
      return(NULL)
    }
    tmb <- as.data.frame(table(patient = substr(maf_data$Tumor_Sample_Barcode, 1L, 12L)),
                         stringsAsFactors = FALSE)
    colnames(tmb) <- c("patient", "total")
    ann <- clinical_ann
    ann$total <- tmb$total[match(ann$Tumor_Sample_Barcode, tmb$patient)]
    ann$total[is.na(ann$total)] <- 0
    grp <- as.character(ann[[group_label]])
    # Order known groups by descending median TMB; NA-group samples go last.
    med <- tapply(ann$total[!is.na(grp)], grp[!is.na(grp)], stats::median, na.rm = TRUE)
    grp_levels <- c(names(sort(med, decreasing = TRUE)), NA)
    ann$grp_key <- factor(grp, levels = grp_levels, exclude = NULL)
    ann <- ann[order(ann$grp_key, -ann$total, ann$Tumor_Sample_Barcode), , drop = FALSE]
    ann$Tumor_Sample_Barcode
  }

  draw_oncoplot <- function(maf_obj, basename, n, clinical_features = NULL, sample_order = NULL) {
    pdf_path <- file.path(ctx$plots_dir, sprintf("%s.pdf", basename))
    png_path <- file.path(ctx$plots_dir, sprintf("%s.png", basename))
    height_pdf <- max(6, n * 0.3 + if (length(clinical_features)) 0.4 * length(clinical_features) else 0)
    height_png <- max(900, n * 40 + if (length(clinical_features)) 30 * length(clinical_features) else 0)
    # Wide canvas so the side legends and annotation labels are not clipped.
    width_pdf <- max(14, if (length(clinical_features)) 16 else 12)
    width_png <- max(2400, if (length(clinical_features)) 3200 else 2000)
    do_plot <- function() {
      if (length(clinical_features)) {
        if (length(sample_order)) {
          # Explicit order: groups by median TMB, samples by TMB within group.
          maftools::oncoplot(maf = maf_obj, top = n,
                             clinicalFeatures = clinical_features,
                             sampleOrder = sample_order,
                             sortByAnnotation = FALSE)
        } else {
          maftools::oncoplot(maf = maf_obj, top = n,
                             clinicalFeatures = clinical_features,
                             sortByAnnotation = TRUE)
        }
      } else {
        maftools::oncoplot(maf = maf_obj, top = n)
      }
    }
    grDevices::pdf(pdf_path, width = width_pdf, height = height_pdf)
    do_plot()
    grDevices::dev.off()
    grDevices::png(png_path, width = width_png, height = height_png, res = 150)
    do_plot()
    grDevices::dev.off()
    invisible(list(pdf = pdf_path, png = png_path))
  }

  # --- global summary --------------------------------------------------------
  maf_obj <- maftools::read.maf(maf = maf_df, isTCGA = TRUE)

  summary_df <- as.data.frame(maftools::getSampleSummary(maf_obj))
  write_table_safe(summary_df, file.path(ctx$results_dir, sprintf("%s_maf_summary.csv", project)))

  gene_summary_df <- as.data.frame(maftools::getGeneSummary(maf_obj))
  write_table_safe(gene_summary_df, file.path(ctx$results_dir, sprintf("%s_gene_mutation_frequency.csv", project)))

  pdf_path <- file.path(ctx$plots_dir, sprintf("%s_maf_summary.pdf", project))
  png_path <- file.path(ctx$plots_dir, sprintf("%s_maf_summary.png", project))
  grDevices::pdf(pdf_path, width = 10, height = 8)
  maftools::plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
  grDevices::dev.off()
  grDevices::png(png_path, width = 1600, height = 1200, res = 150)
  maftools::plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
  grDevices::dev.off()

  # Global (unannotated) oncoplot, one per top_n value.
  for (n in top_n) {
    draw_oncoplot(maf_obj, sprintf("%s_oncoplot_top%s", project, n), n)
  }

  report_bullets <- c(
    sprintf("- Samples in MAF: `%s`", nrow(summary_df)),
    sprintf("- Unique mutated genes: `%s`", nrow(gene_summary_df)),
    sprintf("- Top N oncoplot: `%s`", paste(top_n, collapse = ", ")),
    "",
    "## Outputs",
    sprintf("- `results/%s_maf_summary.csv`", project),
    sprintf("- `results/%s_gene_mutation_frequency.csv`", project),
    sprintf("- `plots/%s_maf_summary.pdf/png`", project),
    sprintf("- `plots/%s_oncoplot_top%s.pdf/png`", project, paste(top_n, collapse = ", "))
  )

  # --- clinical-annotated oncoplot -------------------------------------------
  clinical_ann <- build_clinical_annotation(maf_df)
  if (!is.null(clinical_ann)) {
    maf_clin <- maftools::read.maf(maf = maf_df, clinicalData = clinical_ann, isTCGA = TRUE)
    clin_features <- colnames(clinical_ann)[-1L]
    sample_order <- if (within_group_order == "tmb") compute_tmb_sample_order(maf_df, clinical_ann) else NULL
    for (n in top_n) {
      draw_oncoplot(maf_clin, sprintf("%s_oncoplot_clinical_top%s", project, n),
                    n, clinical_features = clin_features, sample_order = sample_order)
    }
    write_table_safe(clinical_ann, file.path(ctx$results_dir, sprintf("%s_oncoplot_clinical_annotation.csv", project)))
    order_desc <- if (length(sample_order)) "groups by median TMB, samples by TMB within group" else "sortByAnnotation"
    report_bullets <- c(
      report_bullets,
      sprintf("- `plots/%s_oncoplot_clinical_top%s.pdf/png` (annotated by: %s; grouped by: %s; order: %s)",
              project, paste(top_n, collapse = ", "), paste(clin_features, collapse = ", "), group_column, order_desc),
      sprintf("- `results/%s_oncoplot_clinical_annotation.csv`", project)
    )
  }

  # --- LoF (loss-of-function) focus ------------------------------------------
  lof_gene_summary_df <- NULL
  if (lof_focus) {
    lof_df <- maf_df[maf_df$Variant_Classification %in% lof_classes, , drop = FALSE]
    if (!nrow(lof_df)) {
      info("No LoF variants found for %s with classes: %s", project, paste(lof_classes, collapse = ", "))
    } else {
      lof_clin <- if (!is.null(clinical_ann)) clinical_ann else NULL
      maf_lof <- if (!is.null(lof_clin)) {
        maftools::read.maf(maf = lof_df, clinicalData = lof_clin, isTCGA = TRUE)
      } else {
        maftools::read.maf(maf = lof_df, isTCGA = TRUE)
      }
      lof_gene_summary_df <- as.data.frame(maftools::getGeneSummary(maf_lof))
      write_table_safe(lof_gene_summary_df, file.path(ctx$results_dir, sprintf("%s_lof_gene_mutation_frequency.csv", project)))
      clin_features <- if (!is.null(lof_clin)) colnames(lof_clin)[-1L] else NULL
      # Order by LoF burden (not total burden) when requested.
      lof_sample_order <- if (within_group_order == "tmb" && !is.null(lof_clin)) {
        compute_tmb_sample_order(lof_df, lof_clin)
      } else {
        NULL
      }
      for (n in lof_top_n) {
        draw_oncoplot(maf_lof, sprintf("%s_oncoplot_LOF_top%s", project, n),
                      n, clinical_features = clin_features, sample_order = lof_sample_order)
      }
      n_lof_samples <- length(unique(lof_df$Tumor_Sample_Barcode))
      top_lof_genes <- head(lof_gene_summary_df$Hugo_Symbol, max(lof_top_n))
      report_bullets <- c(
        report_bullets,
        "",
        sprintf("## Loss-of-function focus (%s)", paste(lof_classes, collapse = ", ")),
        sprintf("- LoF variants: `%s` across `%s` samples", nrow(lof_df), n_lof_samples),
        sprintf("- Top LoF genes: %s", paste(top_lof_genes, collapse = ", ")),
        sprintf("- `results/%s_lof_gene_mutation_frequency.csv`", project),
        sprintf("- `plots/%s_oncoplot_LOF_top%s.pdf/png`", project, paste(lof_top_n, collapse = ", "))
      )
    }
  }

  write_report(ctx, sprintf("MAF Summary for %s", project), report_bullets)

  list(
    summary = summary_df,
    gene_summary = gene_summary_df,
    lof_gene_summary = lof_gene_summary_df,
    maf_obj = maf_obj
  )
}
