task_run_tme <- function(config, ctx) {
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    fail("run_tme requires the IOBR package.")
  }

  project <- config$project %||% fail("run_tme requires a project.")
  methods <- unlist(config$methods %||% c("estimate"))
  data <- standardize_project_data(project)
  counts <- round(data$counts)

  dat <- IOBR::count2tpm(countMat = counts, idType = "Symbol", org = "hsa", source = "local")
  outputs <- list()

  for (method in methods) {
    info("Running TME method: %s", method)
    if (identical(method, "cibersort")) {
      fit <- IOBR::deconvo_tme(eset = dat, method = method, arrays = FALSE, perm = as.integer(config$perm %||% 100))
    } else {
      fit <- IOBR::deconvo_tme(eset = dat, method = method, arrays = FALSE)
    }
    out_path <- file.path(ctx$results_dir, sprintf("%s_%s_tme.csv", project, method))
    write_table_safe(as.data.frame(fit), out_path)
    outputs[[method]] <- fit
  }

  write_report(
    ctx,
    sprintf("TME Deconvolution for %s", project),
    c(
      sprintf("- Methods executed: `%s`", paste(methods, collapse = ", ")),
      sprintf("- Samples used: `%s`", ncol(counts)),
      "",
      "## Outputs",
      "- One CSV per TME method in `results/`"
    )
  )

  outputs
}
