task_audit_data <- function(config, ctx) {
  root <- project_root()
  data_dir <- file.path(root, "0-Data")
  raw_dir <- file.path(root, "GDCdata")
  projects_from_raw <- if (dir.exists(raw_dir)) sort(list.dirs(raw_dir, recursive = FALSE, full.names = FALSE)) else character()
  projects_from_data <- sort(unique(gsub("_(mrna|clinical|maf)\\.rda$", "", list.files(data_dir, pattern = "_(mrna|clinical|maf)\\.rda$"))))
  projects <- sort(unique(c(projects_from_raw, projects_from_data)))

  gdc_release_version <- config$gdc_release_version %||% "45.0"
  gdc_release_date <- as.POSIXct(config$gdc_release_date %||% "2025-12-04", tz = "UTC")

  audit_rows <- lapply(projects, function(project) {
    files <- project_files(project)
    raw_files <- if (dir.exists(files$raw_dir)) list.files(files$raw_dir, recursive = TRUE, full.names = TRUE) else character()
    mrna_info <- safe_file_info(files$mrna)
    clinical_info <- safe_file_info(files$clinical)
    maf_info <- safe_file_info(files$maf)
    deg_input_info <- safe_file_info(files$deg_input)

    data.frame(
      project = project,
      has_raw_dir = dir.exists(files$raw_dir),
      raw_file_count = length(raw_files),
      has_mrna = mrna_info$exists,
      mrna_modified_at = mrna_info$modified_at,
      mrna_stale_vs_gdc45 = if (mrna_info$exists) as.POSIXct(mrna_info$modified_at, tz = "UTC") < gdc_release_date else NA,
      has_clinical = clinical_info$exists,
      clinical_modified_at = clinical_info$modified_at,
      has_maf = maf_info$exists,
      maf_modified_at = maf_info$modified_at,
      has_deg_input = deg_input_info$exists,
      deg_input_modified_at = deg_input_info$modified_at,
      stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, audit_rows)

  gtex_files <- c(
    exp_gtex_tpm = file.path(root, "GTEX", "exp_gtex.tpm.csv"),
    gtex_rsem = file.path(root, "GTEX", "gtex_RSEM_gene_tpm.gz"),
    tissue_gtex = file.path(root, "GTEX", "tissue_gtex.csv"),
    samplepair = file.path(root, "GTEX", "samplepair.txt"),
    probemap = file.path(root, "GTEX", "gencode.v23.annotation.gene.probemap")
  )
  gtex_audit <- do.call(rbind, lapply(names(gtex_files), function(name) {
    info <- safe_file_info(gtex_files[[name]])
    data.frame(
      file_key = name,
      path = gtex_files[[name]],
      exists = info$exists,
      size_bytes = info$size_bytes,
      modified_at = info$modified_at,
      stringsAsFactors = FALSE
    )
  }))

  write_table_safe(audit, file.path(ctx$results_dir, "project_audit.csv"))
  write_table_safe(gtex_audit, file.path(ctx$results_dir, "gtex_audit.csv"))
  write_json(
    list(
      gdc_release_version = gdc_release_version,
      gdc_release_date = format(gdc_release_date, "%Y-%m-%d"),
      project_count = nrow(audit),
      projects_with_mrna = sum(audit$has_mrna),
      projects_with_raw = sum(audit$has_raw_dir)
    ),
    file.path(ctx$results_dir, "summary.json")
  )

  write_report(
    ctx,
    "TCGA Data Audit",
    c(
      sprintf("- GDC reference release used for staleness flag: `%s` (%s)", gdc_release_version, format(gdc_release_date, "%Y-%m-%d")),
      sprintf("- Projects detected: `%s`", nrow(audit)),
      sprintf("- Projects with local expression RDA: `%s`", sum(audit$has_mrna)),
      sprintf("- Projects with raw GDC directories: `%s`", sum(audit$has_raw_dir)),
      sprintf("- GTEx files present: `%s / %s`", sum(gtex_audit$exists), nrow(gtex_audit)),
      "",
      "## Outputs",
      "- `results/project_audit.csv`",
      "- `results/gtex_audit.csv`",
      "- `results/summary.json`"
    )
  )

  list(project_audit = audit, gtex_audit = gtex_audit)
}
