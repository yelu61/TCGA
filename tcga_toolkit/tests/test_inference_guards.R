#!/usr/bin/env Rscript

# Adversarial fixtures for the inference guards added to common.R,
# task_run_deg.R, task_prognostic_model.R and task_external_validate.R:
# assay fallback must fail, non-integer counts must fail, survival models must
# dedup patients and report PH checks, and external validation must honour a
# locked training cutoff.

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), mustWork = TRUE)
toolkit_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
scripts_dir <- file.path(toolkit_root, "scripts")

options(tcga_toolkit.root = toolkit_root)

passes <- 0L
failures <- 0L
log <- function(name, ok, detail = "") {
  if (ok) {
    passes <<- passes + 1L
    cat(sprintf("  PASS  %s\n", name))
  } else {
    failures <<- failures + 1L
    cat(sprintf("  FAIL  %s\n", name))
    if (nzchar(detail)) cat(sprintf("        %s\n", detail))
  }
}

source(file.path(scripts_dir, "common.R"))
source(file.path(scripts_dir, "task_external_validate.R"))

throws <- function(expr, pattern) {
  err <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  !is.null(err) && grepl(pattern, err)
}

# ---- extract_assay: specified-but-missing assay must fail -------------------
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm_unstrand = matrix(1:6, nrow = 2,
                                        dimnames = list(paste0("g", 1:2), paste0("s", 1:3))))
  )
  log("extract_assay fails on missing requested assay",
      throws(extract_assay(se, "unstranded"), "Assay 'unstranded' not found"))
  log("extract_assay returns requested assay when present",
      identical(dim(extract_assay(se, "tpm_unstrand")), c(2L, 3L)))
  log("extract_assay falls back to first assay only when unspecified",
      identical(dim(extract_assay(se, NULL)), c(2L, 3L)))
} else {
  log("extract_assay tests (SummarizedExperiment installed)", FALSE, "package missing")
}

# ---- assert_integer_counts: TPM-like input must fail ------------------------
tpm_like <- matrix(c(10.5, 20.25, 30.1, 40.0), nrow = 2)
int_like <- matrix(c(10, 20, 30, 40), nrow = 2)
log("assert_integer_counts rejects TPM-like values",
    throws(assert_integer_counts(tpm_like, "tpm fixture"), "non-integer"))
log("assert_integer_counts accepts integer-like values",
    isTRUE(assert_integer_counts(int_like, "int fixture")))
log("assert_integer_counts rejects NA",
    throws(assert_integer_counts(matrix(c(1, NA), nrow = 1), "na fixture"), "NA"))

# ---- run_survival_models: patient dedup + PH check + meta -------------------
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(1)
  n_pat <- 20
  # Two aliquots for the first 5 patients -> 25 rows, 20 patients.
  score_df <- data.frame(
    sample_id = c(sprintf("P%02dA", 1:n_pat), sprintf("P%02dB", 1:5)),
    score = rnorm(n_pat + 5),
    stringsAsFactors = FALSE
  )
  clinical <- data.frame(
    sample_id = score_df$sample_id,
    patient_short = c(sprintf("P%02d", 1:n_pat), sprintf("P%02d", 1:5)),
    vital_status = rep(c("dead", "alive"), length.out = n_pat + 5),
    days_to_death = ifelse(rep(c(TRUE, FALSE), length.out = n_pat + 5),
                           sample(100:900, n_pat + 5, replace = TRUE), NA),
    days_to_last_follow_up = sample(200:1500, n_pat + 5, replace = TRUE),
    stringsAsFactors = FALSE
  )
  res <- run_survival_models(score_df, clinical)
  log("run_survival_models dedups patients",
      !is.null(res) && res$km$n == n_pat,
      detail = if (is.null(res)) "returned NULL" else sprintf("n=%s", res$km$n))
  log("run_survival_models reports PH check",
      !is.null(res$ph) && all(c("term", "ph_p_value") %in% colnames(res$ph)))
  log("run_survival_models declares time origin and event definition",
      !is.null(res$meta) && grepl("days", res$meta$time_origin) && res$meta$patient_dedup == 5L)
} else {
  log("run_survival_models tests (survival installed)", FALSE, "package missing")
}

# ---- external_validate: locked cutoff and transform -------------------------
tmp <- tempfile()
dir.create(file.path(tmp, "results"), recursive = TRUE)
dir.create(file.path(tmp, "plots"), recursive = TRUE)
dir.create(file.path(tmp, "objects"), recursive = TRUE)
ctx <- list(results_dir = file.path(tmp, "results"),
            plots_dir = file.path(tmp, "plots"),
            objects_dir = file.path(tmp, "objects"),
            run_dir = tmp,
            task = "external_validate",
            task_id = "fixture")

set.seed(2)
n_ext <- 40
genes <- c("GENEA", "GENEB", "GENEC")
expr_mat <- matrix(rexp(length(genes) * n_ext, rate = 0.1), nrow = length(genes),
                   dimnames = list(genes, sprintf("EXT%02d", 1:n_ext)))
expr_file <- file.path(tmp, "expr.csv")
write.csv(expr_mat, expr_file, quote = FALSE)

clin <- data.frame(
  sample_id = sprintf("EXT%02d", 1:n_ext),
  OS.time = sample(100:2000, n_ext, replace = TRUE),
  OS = rep(c(1, 0), length.out = n_ext),
  stringsAsFactors = FALSE
)
clin_file <- file.path(tmp, "clinical.csv")
write.csv(clin, clin_file, row.names = FALSE)

locked_cutoff <- 1234.5
meta_file <- file.path(tmp, "meta.json")
jsonlite::write_json(
  list(transform = "log2p1",
       features = as.list(genes),
       coefficients = list(GENEA = 0.8, GENEB = -0.5, GENEC = 0.2),
       cutoff = locked_cutoff),
  meta_file, auto_unbox = TRUE
)

out <- task_external_validate(
  list(cohort_name = "FIXTURE",
       expression_file = expr_file,
       clinical_file = clin_file,
       model_meta_file = meta_file),
  ctx
)

km_json <- jsonlite::read_json(file.path(tmp, "results", "FIXTURE_external_km.json"),
                               simplifyVector = TRUE)
log("external_validate honours locked training cutoff",
    isTRUE(all.equal(km_json$cutoff, locked_cutoff)),
    detail = sprintf("cutoff=%s", km_json$cutoff))
log("external_validate labels cutoff source as locked",
    identical(km_json$cutoff_source, "locked from training"))

# Without locked metadata the cutoff must be flagged as exploratory.
out2 <- task_external_validate(
  list(cohort_name = "FIXTURE2",
       expression_file = expr_file,
       clinical_file = clin_file,
       weight_file = NULL,
       feature_genes = genes),
  ctx
)
km_json2 <- jsonlite::read_json(file.path(tmp, "results", "FIXTURE2_external_km.json"),
                                simplifyVector = TRUE)
log("external_validate flags re-derived cutoff as exploratory",
    grepl("EXPLORATORY", km_json2$cutoff_source))

cat(sprintf("\n%4d passed, %d failed\n", passes, failures))
if (failures > 0) quit(status = 1L)
