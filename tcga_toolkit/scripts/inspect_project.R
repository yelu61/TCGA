#!/usr/bin/env Rscript

# inspect_project.R
# Print a quick summary of a TCGA/TARGET project so an agent can choose group
# columns, subtypes and other config fields without loading the rda by hand.
#
# Usage:
#   Rscript tcga_toolkit/scripts/inspect_project.R --project TCGA-BRCA
#   Rscript tcga_toolkit/scripts/inspect_project.R --project TCGA-BRCA --json
#   Rscript tcga_toolkit/scripts/inspect_project.R --list

args <- commandArgs(trailingOnly = TRUE)
project <- NULL
as_json <- FALSE
list_only <- FALSE
for (i in seq_along(args)) {
  if (args[[i]] == "--project" && i < length(args)) {
    project <- args[[i + 1L]]
  }
  if (args[[i]] == "--json") {
    as_json <- TRUE
  }
  if (args[[i]] == "--list") {
    list_only <- TRUE
  }
}

script_flag <- grep("^--file=", commandArgs(), value = TRUE)[1]
if (is.na(script_flag)) {
  stop("Unable to resolve script path.", call. = FALSE)
}
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_flag), perl = TRUE)
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
options(tcga_toolkit.root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
source(file.path(script_dir, "common.R"))

if (list_only) {
  projs <- available_projects(include_target = TRUE)
  if (as_json) {
    cat(jsonlite::toJSON(projs, auto_unbox = FALSE, pretty = TRUE), "\n")
  } else {
    cat(paste(projs, collapse = "\n"), "\n", sep = "")
  }
  quit(status = 0)
}

if (is.null(project)) {
  stop("Usage: Rscript tcga_toolkit/scripts/inspect_project.R --project <TCGA-XXX> [--json] [--list]", call. = FALSE)
}

files <- project_files(project)
file_status <- list(
  mrna = safe_file_info(files$mrna),
  clinical = safe_file_info(files$clinical),
  maf = safe_file_info(files$maf),
  deg_input = safe_file_info(files$deg_input),
  raw_dir_exists = dir.exists(files$raw_dir)
)

result <- list(project = project, files = file_status)

if (file.exists(files$mrna) || file.exists(files$deg_input)) {
  data <- tryCatch(standardize_project_data(project), error = function(e) {
    message(sprintf("Failed to standardize %s: %s", project, conditionMessage(e)))
    NULL
  })

  if (!is.null(data)) {
    clinical <- data$clinical
    counts <- data$counts

    group_counts <- if ("Group" %in% colnames(clinical)) {
      tab <- table(as.character(clinical$Group), useNA = "ifany")
      as.list(tab)
    } else {
      list()
    }

    candidate_cols <- intersect(
      c("Group", "gender", "race", "ethnicity", "vital_status",
        "age_at_index", "age_at_diagnosis", "ajcc_pathologic_stage",
        "tumor_stage", "primary_diagnosis", "tumor_descriptor",
        "definition", "sample_type", "shortLetterCode",
        "days_to_death", "days_to_last_follow_up"),
      colnames(clinical)
    )

    clinical_cols <- lapply(candidate_cols, function(col) {
      v <- as.character(clinical[[col]])
      uniq <- unique(v[!is.na(v) & nzchar(v)])
      list(
        column = col,
        n_non_missing = sum(!is.na(v) & nzchar(v)),
        unique_values = length(uniq),
        sample_values = utils::head(sort(uniq), 12)
      )
    })

    subtypes <- tryCatch({
      sub <- join_project_subtypes(data, subtype_column = "Subtype_Selected")
      cnt <- table(as.character(sub$subtype), useNA = "ifany")
      list(column = "Subtype_Selected", counts = as.list(cnt))
    }, error = function(e) {
      list(column = "Subtype_Selected", counts = list(), note = conditionMessage(e))
    })

    survival_ready <- "vital_status" %in% colnames(clinical) &&
      (any(c("days_to_death", "days_to_last_follow_up") %in% colnames(clinical)))

    result$samples <- list(
      total = ncol(counts),
      genes = nrow(counts),
      groups = group_counts,
      survival_columns_present = survival_ready
    )
    result$clinical_columns <- clinical_cols
    result$subtypes <- subtypes
    result$data_source <- data$source
  }
}

if (as_json) {
  cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")
} else {
  cat(sprintf("Project: %s\n", result$project))
  cat("\nFile presence:\n")
  for (k in setdiff(names(result$files), "raw_dir_exists")) {
    info <- result$files[[k]]
    cat(sprintf("  %-12s exists=%s modified=%s size=%s\n",
                k,
                ifelse(isTRUE(info$exists), "yes", "no"),
                info$modified_at %||% "NA",
                info$size_bytes %||% "NA"))
  }
  cat(sprintf("  raw_dir      exists=%s\n", ifelse(isTRUE(result$files$raw_dir_exists), "yes", "no")))

  if (!is.null(result$samples)) {
    cat("\nSamples:\n")
    cat(sprintf("  total: %s\n", result$samples$total))
    cat(sprintf("  genes: %s\n", result$samples$genes))
    if (length(result$samples$groups)) {
      cat("  Group counts:\n")
      for (g in names(result$samples$groups)) {
        cat(sprintf("    %s: %s\n", g, result$samples$groups[[g]]))
      }
    }
    cat(sprintf("  survival columns present: %s\n", result$samples$survival_columns_present))
  }

  if (length(result$clinical_columns)) {
    cat("\nClinical columns (non-missing / unique):\n")
    for (c in result$clinical_columns) {
      cat(sprintf("  %-28s n=%s unique=%s sample=%s\n",
                  c$column, c$n_non_missing, c$unique_values,
                  paste(utils::head(c$sample_values, 6), collapse = ", ")))
    }
  }

  if (length(result$subtypes$counts)) {
    cat(sprintf("\nSubtypes (%s):\n", result$subtypes$column))
    for (s in names(result$subtypes$counts)) {
      cat(sprintf("  %s: %s\n", s, result$subtypes$counts[[s]]))
    }
  }
}
