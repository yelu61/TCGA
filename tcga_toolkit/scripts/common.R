suppressPackageStartupMessages({
  library(jsonlite)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

fail <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

info <- function(...) {
  message(sprintf(...))
}

toolkit_root <- function() {
  root <- getOption("tcga_toolkit.root")
  if (is.null(root)) {
    fail("Option tcga_toolkit.root is not set.")
  }
  normalizePath(root, mustWork = TRUE)
}

toolkit_version <- function() {
  cached <- getOption("tcga_toolkit.version")
  if (!is.null(cached)) {
    return(cached)
  }
  version_path <- file.path(toolkit_root(), "VERSION")
  version <- if (file.exists(version_path)) {
    trimws(readLines(version_path, warn = FALSE, n = 1L))
  } else {
    "0.0.0-dev"
  }
  options(tcga_toolkit.version = version)
  version
}

project_root <- function() {
  normalizePath(file.path(toolkit_root(), ".."), mustWork = TRUE)
}

output_root <- function() {
  or <- getOption("tcga_toolkit.output_root")
  if (is.null(or)) {
    return(project_root())
  }
  normalizePath(or, mustWork = TRUE)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

timestamp_id <- function() {
  format(Sys.time(), "%Y%m%d-%H%M%S")
}

normalize_task_id <- function(task, explicit = NULL) {
  slug <- gsub("[^A-Za-z0-9_-]+", "-", explicit %||% "")
  slug <- gsub("(^-+|-+$)", "", slug)
  if (!nzchar(slug)) {
    slug <- task
  }
  sprintf("%s-%s", timestamp_id(), slug)
}

read_text_lines <- function(path, n = -1L) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

read_config <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("json", "jsonc")) {
    return(fromJSON(path, simplifyVector = FALSE))
  }
  if (ext %in% c("yml", "yaml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      fail("Config %s is YAML but package 'yaml' is not installed.", path)
    }
    return(yaml::read_yaml(path))
  }
  fail("Unsupported config format: %s", path)
}

write_json <- function(x, path) {
  write_json_impl <- jsonlite::write_json
  write_json_impl(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(path)
}

write_table_safe <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

load_rdata_file <- function(path) {
  env <- new.env(parent = emptyenv())
  objs <- load(path, envir = env)
  if (!length(objs)) {
    fail("No objects found in %s", path)
  }
  if (length(objs) == 1L) {
    return(env[[objs[[1L]]]])
  }
  structure(
    lapply(objs, function(name) env[[name]]),
    names = objs
  )
}

safe_file_info <- function(path) {
  if (!file.exists(path)) {
    return(list(exists = FALSE, size_bytes = NA_real_, modified_at = NA_character_))
  }
  info <- file.info(path)
  list(
    exists = TRUE,
    size_bytes = unname(info$size[[1L]]),
    modified_at = format(info$mtime[[1L]], "%Y-%m-%d %H:%M:%S")
  )
}

project_files <- function(project) {
  list(
    mrna = file.path(project_root(), "0-Data", sprintf("%s_mrna.rda", project)),
    clinical = file.path(project_root(), "0-Data", sprintf("%s_clinical.rda", project)),
    maf = file.path(project_root(), "0-Data", sprintf("%s_maf.rda", project)),
    deg_input = file.path(project_root(), "1-Input", sprintf("%s_DEGinput.rda", project)),
    raw_dir = file.path(project_root(), "GDCdata", project)
  )
}

available_projects <- function(include_target = TRUE) {
  data_dir <- file.path(project_root(), "0-Data")
  projects <- sort(unique(gsub("_(mrna|clinical|maf)\\.rda$", "", list.files(data_dir, pattern = "_(mrna|clinical|maf)\\.rda$"))))
  if (!include_target) {
    projects <- projects[grepl("^TCGA-", projects)]
  }
  projects
}

derive_group_from_clinical <- function(df) {
  candidates <- c(
    as.character(df$Group %||% character()),
    as.character(df$definition %||% character()),
    as.character(df$sample_type %||% character()),
    as.character(df$shortLetterCode %||% character())
  )
  if (!length(candidates)) {
    return(rep("Unknown", nrow(df)))
  }

  pick <- function(i) {
    values <- c(
      as.character(df$Group %||% NA_character_)[i],
      as.character(df$definition %||% NA_character_)[i],
      as.character(df$sample_type %||% NA_character_)[i],
      as.character(df$shortLetterCode %||% NA_character_)[i]
    )
    values <- values[!is.na(values) & nzchar(values)]
    if (!length(values)) {
      return("Unknown")
    }
    value <- values[[1L]]
    lower <- tolower(value)
    if (grepl("normal", lower, fixed = TRUE) || grepl("solid tissue normal", lower, fixed = TRUE) || lower %in% c("nt")) {
      return("Normal")
    }
    if (grepl("tumor", lower, fixed = TRUE) || grepl("cancer", lower, fixed = TRUE) || lower %in% c("tp", "tb")) {
      return("Tumor")
    }
    value
  }

  vapply(seq_len(nrow(df)), pick, character(1))
}

dedupe_gene_rows <- function(mat, gene_map, reference_mat = NULL) {
  ref <- reference_mat %||% mat
  gene_name <- as.character(gene_map$gene_name %||% rownames(mat))
  gene_name[is.na(gene_name) | gene_name == ""] <- as.character(gene_map$gene_id[is.na(gene_name) | gene_name == ""])
  gene_name[is.na(gene_name) | gene_name == ""] <- rownames(mat)[is.na(gene_name) | gene_name == ""]
  gene_name[is.na(gene_name) | gene_name == ""] <- sprintf("feature_%s", seq_along(gene_name))[is.na(gene_name) | gene_name == ""]

  mean_signal <- rowMeans(ref, na.rm = TRUE)
  order_idx <- order(mean_signal, decreasing = TRUE, na.last = TRUE)
  keep <- !duplicated(gene_name[order_idx])
  keep_idx <- order_idx[keep]

  mat2 <- mat[keep_idx, , drop = FALSE]
  gene_map2 <- gene_map[keep_idx, , drop = FALSE]
  rownames(mat2) <- gene_name[keep_idx]
  gene_map2$feature_id <- rownames(mat2)

  list(matrix = mat2, gene_map = gene_map2)
}

extract_assay <- function(se, preferred) {
  assay_names <- names(SummarizedExperiment::assays(se))
  if (length(preferred) && preferred %in% assay_names) {
    return(SummarizedExperiment::assay(se, preferred))
  }
  if (!length(assay_names)) {
    fail("No assays found in SummarizedExperiment object.")
  }
  SummarizedExperiment::assay(se, assay_names[[1L]])
}

standardize_project_data <- function(project, counts_assay = "unstranded", tpm_assay = "tpm_unstrand") {
  files <- project_files(project)
  source_path <- if (file.exists(files$deg_input)) files$deg_input else files$mrna
  cache_path <- project_data_cache_path(project)

  if (file.exists(cache_path) && file.exists(source_path)) {
    cache_mtime <- file.info(cache_path)$mtime
    source_mtime <- file.info(source_path)$mtime
    if (!is.na(cache_mtime) && !is.na(source_mtime) && cache_mtime > source_mtime) {
      info("Using cached standardized data for %s", project)
      return(readRDS(cache_path))
    }
  }

  result <- NULL

  if (file.exists(files$deg_input)) {
    prepared <- load_rdata_file(files$deg_input)
    if (is.list(prepared) && all(c("Counts", "TPM") %in% names(prepared))) {
      counts <- as.matrix(prepared$Counts)
      tpm <- as.matrix(prepared$TPM)
      clinical <- prepared$clinical_clean %||% prepared$clinical
      clinical <- as.data.frame(clinical, stringsAsFactors = FALSE)
      sample_ids <- colnames(counts)
      candidate_cols <- intersect(c("sample_id", "Sample", "sample", "Patient", "patient", "barcode"), colnames(clinical))
      if (length(candidate_cols)) {
        overlaps <- vapply(candidate_cols, function(col) sum(as.character(clinical[[col]]) %in% sample_ids), numeric(1))
        best_col <- candidate_cols[[which.max(overlaps)]]
        clinical$sample_id <- as.character(clinical[[best_col]])
      } else {
        clinical$sample_id <- sample_ids
      }
      clinical$Sample <- clinical$Sample %||% clinical$sample %||% clinical$sample_id
      clinical$project_id <- clinical$project_id %||% project
      clinical$patient_barcode <- as.character(clinical$Patient %||% clinical$patient %||% clinical$sample_id)
      clinical$patient_short <- substr(clinical$patient_barcode, 1, 12)
      if (!("Group" %in% colnames(clinical))) {
        clinical$Group <- derive_group_from_clinical(clinical)
      }
      rownames(counts) <- rownames(counts) %||% seq_len(nrow(counts))
      rownames(tpm) <- rownames(counts)
      gene_map <- prepared$id2symbol %||% data.frame(
        gene_id = rownames(counts),
        gene_name = rownames(counts),
        gene_type = NA_character_,
        stringsAsFactors = FALSE
      )
      if (!("gene_name" %in% colnames(gene_map))) {
        gene_map$gene_name <- rownames(counts)
      }
      if (nrow(gene_map) != nrow(counts)) {
        gene_key <- if ("gene_name" %in% colnames(gene_map)) gene_map$gene_name else gene_map[[1]]
        gene_map <- gene_map[match(rownames(counts), gene_key), , drop = FALSE]
      }
      if (!("gene_id" %in% colnames(gene_map))) {
        gene_map$gene_id <- rownames(counts)
      }
      if (!("gene_type" %in% colnames(gene_map))) {
        gene_map$gene_type <- NA_character_
      }
      gene_map$feature_id <- rownames(counts)
      result <- list(
        project = project,
        counts = counts,
        tpm = tpm,
        clinical = clinical,
        gene_map = gene_map,
        source = basename(files$deg_input)
      )
    }
  }

  if (is.null(result)) {
    if (!file.exists(files$mrna)) {
      fail("Missing expression file for %s: %s", project, files$mrna)
    }

    obj <- load_rdata_file(files$mrna)
    if (!inherits(obj, "RangedSummarizedExperiment")) {
      fail("Expression file for %s is not a RangedSummarizedExperiment.", project)
    }

    counts <- extract_assay(obj, counts_assay)
    tpm <- extract_assay(obj, tpm_assay)
    cd <- as.data.frame(SummarizedExperiment::colData(obj))
    rd <- as.data.frame(SummarizedExperiment::rowData(obj))

    if (!ncol(counts)) {
      fail("Counts matrix for %s is empty.", project)
    }

    sample_ids <- colnames(obj)
    if (is.null(sample_ids) || !length(sample_ids)) {
      sample_ids <- colnames(counts)
    }

    gene_map <- data.frame(
      gene_id = as.character(rd$gene_id %||% rownames(counts)),
      gene_name = as.character(rd$gene_name %||% rownames(counts)),
      gene_type = as.character(rd$gene_type %||% rep(NA_character_, nrow(rd))),
      stringsAsFactors = FALSE
    )

    deduped <- dedupe_gene_rows(counts, gene_map, reference_mat = tpm)
    counts <- deduped$matrix
    gene_map <- deduped$gene_map
    tpm <- tpm[match(gene_map$gene_id, as.character(rd$gene_id)), , drop = FALSE]
    rownames(tpm) <- gene_map$feature_id

    clinical <- data.frame(
      project_id = project,
      sample_id = sample_ids,
      Sample = sample_ids,
      patient = as.character(cd$patient %||% cd$Patient %||% sample_ids),
      sample = as.character(cd$sample %||% sample_ids),
      barcode = as.character(cd$barcode %||% sample_ids),
      definition = as.character(cd$definition %||% cd$sample_type %||% ""),
      sample_type = as.character(cd$sample_type %||% ""),
      shortLetterCode = as.character(cd$shortLetterCode %||% ""),
      stringsAsFactors = FALSE
    )
    clinical$patient_barcode <- clinical$patient
    clinical$patient_short <- substr(clinical$patient_barcode, 1, 12)

    preferred_cols <- c(
      "gender", "race", "ethnicity", "vital_status", "days_to_death",
      "days_to_last_follow_up", "age_at_index", "age_at_diagnosis",
      "primary_diagnosis", "ajcc_pathologic_stage", "tumor_stage",
      "tumor_descriptor", "sample_submitter_id"
    )
    for (col in preferred_cols) {
      if (col %in% colnames(cd)) {
        clinical[[col]] <- cd[[col]]
      }
    }
    clinical$Group <- derive_group_from_clinical(cbind(clinical, cd))

    result <- list(
      project = project,
      counts = counts,
      tpm = tpm,
      clinical = clinical,
      gene_map = gene_map,
      source = basename(files$mrna)
    )
  }

  saveRDS(result, cache_path)
  result
}

gene_alias_map <- list(
  APITD1 = c("CENPS", "FAAP16"),
  BRE = c("BABAM2", "BRCC45"),
  STRA13 = c("CENPX", "FAAP10"),
  CENPS_CORT = c("CENPS"),
  CENPX = c("CENPX", "FAAP10"),
  FAAP20 = c("FAAP20", "FAAP20-MARCHF8"),
  MAD2L2 = c("MAD2L2", "REV7"),
  TELO2 = c("TELO2", "KIAA0683"),
  C1ORF86 = c("FAAP20"),
  C17ORF70 = c("FAAP100"),
  C19ORF40 = c("FAAP24"),
  KIAA1018 = c("FAN1")
)

resolve_feature_id <- function(gene, mat, gene_map = NULL) {
  if (is.null(gene) || is.na(gene) || !nzchar(gene)) {
    return(NULL)
  }
  if (gene %in% rownames(mat)) {
    return(gene)
  }
  if (!is.null(gene_map)) {
    if ("gene_id" %in% colnames(gene_map)) {
      idx <- match(gene, gene_map$gene_id)
      if (!is.na(idx)) {
        return(gene_map$feature_id[[idx]] %||% rownames(mat)[idx])
      }
    }
    if ("gene_name" %in% colnames(gene_map)) {
      idx <- match(gene, gene_map$gene_name)
      if (!is.na(idx)) {
        return(gene_map$feature_id[[idx]] %||% rownames(mat)[idx])
      }
    }
  }
  aliases <- gene_alias_map[[gene]]
  if (!is.null(aliases)) {
    for (alias in aliases) {
      if (alias %in% rownames(mat)) {
        return(alias)
      }
      if (!is.null(gene_map) && "gene_name" %in% colnames(gene_map)) {
        idx <- match(alias, gene_map$gene_name)
        if (!is.na(idx)) {
          return(gene_map$feature_id[[idx]] %||% rownames(mat)[idx])
        }
      }
    }
  }
  NULL
}

resolve_expression_vector <- function(data, gene, use_tpm = TRUE, log2_transform = TRUE) {
  mat <- if (use_tpm) data$tpm else data$counts
  gene_map <- data$gene_map

  feature_id <- resolve_feature_id(gene, mat, gene_map)
  if (is.null(feature_id)) {
    fail("Gene %s not found in project %s", gene, data$project)
  }

  values <- as.numeric(mat[feature_id, , drop = TRUE])
  names(values) <- colnames(mat)
  if (log2_transform && use_tpm) {
    values <- log2(values + 1)
  }
  data.frame(
    sample_id = names(values),
    feature = gene,
    feature_id = feature_id,
    value = values,
    stringsAsFactors = FALSE
  )
}

feature_label_from_config <- function(config) {
  if (!is.null(config$gene)) {
    return(as.character(config$gene))
  }
  if (!is.null(config$signature_name)) {
    return(as.character(config$signature_name))
  }
  if (!is.null(config$signature_file)) {
    return(tools::file_path_sans_ext(basename(config$signature_file)))
  }
  fail("Config requires gene or signature_file.")
}

compute_feature_scores_from_config <- function(data, config) {
  if (!is.null(config$gene)) {
    scores <- resolve_expression_vector(data, config$gene)
    names(scores)[names(scores) == "value"] <- "score"
    scores$n_genes <- 1L
    scores$feature_label <- as.character(config$gene)
    return(scores[, c("sample_id", "feature_label", "score", "n_genes"), drop = FALSE])
  }
  if (!is.null(config$signature_file)) {
    markers <- read_marker_table(config$signature_file)
    scores <- compute_signature_score(data, marker_table = markers, method = config$signature_method %||% "mean_z")
    scores$feature_label <- feature_label_from_config(config)
    return(scores[, c("sample_id", "feature_label", "score", "n_genes"), drop = FALSE])
  }
  fail("Config requires gene or signature_file.")
}

correlate_target_with_gene_set <- function(data, target_gene, genes, method = "pearson", sample_filter = c("Tumor")) {
  scores_target <- compute_feature_scores_from_config(data, list(gene = target_gene))
  clinical <- data$clinical
  keep <- rep(TRUE, nrow(clinical))
  if (!is.null(sample_filter) && length(sample_filter) && "Group" %in% colnames(clinical)) {
    keep <- as.character(clinical$Group) %in% sample_filter
  }
  clinical <- clinical[keep, , drop = FALSE]
  target <- scores_target[scores_target$sample_id %in% clinical$sample_id, c("sample_id", "score"), drop = FALSE]
  names(target)[2] <- "target_score"

  rows <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    gene <- genes[[i]]
    rows[[i]] <- tryCatch({
      gene_df <- resolve_expression_vector(data, gene)
      gene_df <- gene_df[gene_df$sample_id %in% clinical$sample_id, c("sample_id", "value"), drop = FALSE]
      names(gene_df)[2] <- "gene_score"
      merged <- merge(target, gene_df, by = "sample_id")
      merged <- merged[stats::complete.cases(merged), , drop = FALSE]
      if (nrow(merged) < 10 || stats::sd(merged$target_score) == 0 || stats::sd(merged$gene_score) == 0) {
        data.frame(gene = gene, n = nrow(merged), correlation = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
      } else {
        test <- suppressWarnings(stats::cor.test(merged$target_score, merged$gene_score, method = method))
        data.frame(
          gene = gene,
          n = nrow(merged),
          correlation = unname(test$estimate),
          p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      data.frame(gene = gene, n = 0, correlation = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
    })
  }
  do.call(rbind, rows)
}

read_marker_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  gene_col <- intersect(c("gene", "Gene", "symbol", "Symbol"), colnames(x))
  if (!length(gene_col)) {
    fail("Marker/signature file %s must contain a gene column.", path)
  }
  names(x)[match(gene_col[[1L]], names(x))] <- "gene"
  x
}

read_gene_list_file <- function(path, sheet = NULL, gene_column = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      fail("Reading Excel gene lists requires the readxl package.")
    }
    args <- list(path = path)
    if (!is.null(sheet)) {
      args$sheet <- sheet
    }
    x <- do.call(readxl::read_excel, args)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else if (ext == "csv") {
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (!nrow(x)) {
    fail("Gene list file %s is empty.", path)
  }

  if (is.null(gene_column)) {
    candidates <- c("Gene Name", "gene", "Gene", "symbol", "Symbol")
    gene_column <- intersect(candidates, colnames(x))[1]
    if (is.na(gene_column) || is.null(gene_column)) {
      gene_column <- colnames(x)[1]
    }
  }
  if (!(gene_column %in% colnames(x))) {
    fail("Gene column %s not found in %s", gene_column, path)
  }

  genes_raw <- as.character(x[[gene_column]])
  genes_raw <- genes_raw[!is.na(genes_raw) & nzchar(trimws(genes_raw))]
  genes_clean <- trimws(genes_raw)
  genes_clean <- sub("\\s*\\(.*\\)$", "", genes_clean)
  genes_clean <- sub("\\s+.*$", "", genes_clean)
  genes_clean <- unique(genes_clean[nzchar(genes_clean)])

  data.frame(
    gene = genes_clean,
    raw_label = genes_clean,
    stringsAsFactors = FALSE
  )
}

compute_signature_score <- function(data, marker_table, method = "mean_z") {
  genes <- unique(marker_table$gene)
  feature_ids <- vapply(genes, function(gene) {
    if (gene %in% rownames(data$tpm)) {
      return(gene)
    }
    idx <- match(gene, data$gene_map$gene_id)
    if (!is.na(idx)) {
      return(data$gene_map$feature_id[[idx]])
    }
    idx <- match(gene, data$gene_map$gene_name)
    if (!is.na(idx)) {
      return(data$gene_map$feature_id[[idx]])
    }
    NA_character_
  }, character(1))
  feature_ids <- feature_ids[!is.na(feature_ids)]
  if (!length(feature_ids)) {
    fail("None of the signature genes are available in %s.", data$project)
  }

  mat <- log2(data$tpm[feature_ids, , drop = FALSE] + 1)
  weights <- rep(1, nrow(mat))
  names(weights) <- rownames(mat)
  if ("weight" %in% colnames(marker_table)) {
    for (i in seq_len(nrow(marker_table))) {
      gene <- marker_table$gene[[i]]
      idx <- match(gene, data$gene_map$gene_name)
      if (!is.na(idx)) {
        fid <- data$gene_map$feature_id[[idx]]
        if (fid %in% names(weights)) {
          weights[[fid]] <- as.numeric(marker_table$weight[[i]] %||% 1)
        }
      }
    }
  }
  if ("direction" %in% colnames(marker_table)) {
    for (i in seq_len(nrow(marker_table))) {
      gene <- marker_table$gene[[i]]
      idx <- match(gene, data$gene_map$gene_name)
      if (!is.na(idx)) {
        fid <- data$gene_map$feature_id[[idx]]
        if (fid %in% names(weights)) {
          if (tolower(as.character(marker_table$direction[[i]])) %in% c("down", "-", "negative")) {
            weights[[fid]] <- -abs(weights[[fid]])
          }
        }
      }
    }
  }

  if (identical(method, "mean_z")) {
    z <- t(scale(t(mat)))
    z[is.na(z)] <- 0
    score <- colMeans(z * weights[rownames(z)], na.rm = TRUE)
  } else {
    score <- colMeans(mat * weights[rownames(mat)], na.rm = TRUE)
  }

  data.frame(
    sample_id = names(score),
    score = as.numeric(score),
    n_genes = nrow(mat),
    stringsAsFactors = FALSE
  )
}

select_group_samples <- function(data, group_column = "Group", numerator = NULL, denominator = NULL) {
  clinical <- data$clinical
  if (!(group_column %in% colnames(clinical))) {
    fail("Group column %s not found in clinical table.", group_column)
  }
  keep_levels <- Filter(function(x) !is.null(x) && nzchar(x), c(numerator, denominator))
  if (!length(keep_levels)) {
    keep <- !is.na(clinical[[group_column]]) & nzchar(as.character(clinical[[group_column]]))
  } else {
    keep <- as.character(clinical[[group_column]]) %in% keep_levels
  }
  clinical[keep, , drop = FALSE]
}

coerce_status <- function(x) {
  lower <- tolower(as.character(x))
  ifelse(
    lower %in% c("dead", "deceased", "1", "true"),
    1,
    ifelse(lower %in% c("alive", "0", "false"), 0, NA_real_)
  )
}

survival_columns <- function(clinical) {
  time_cols <- c("days_to_death", "days_to_last_follow_up")
  status_cols <- c("vital_status", "status", "event")
  list(
    time_death = time_cols[time_cols %in% colnames(clinical)][1] %||% NA_character_,
    time_followup = time_cols[time_cols %in% colnames(clinical)][2] %||% NA_character_,
    status = status_cols[status_cols %in% colnames(clinical)][1] %||% NA_character_
  )
}

run_survival_models <- function(score_df, clinical) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    return(NULL)
  }
  merged <- merge(score_df, clinical, by = "sample_id", all.x = TRUE)
  cols <- survival_columns(merged)
  if (all(is.na(unlist(cols)))) {
    return(NULL)
  }
  time <- suppressWarnings(as.numeric(merged[[cols$time_death]]))
  if (!is.na(cols$time_followup) && cols$time_followup %in% colnames(merged) && !all(is.na(merged[[cols$time_followup]]))) {
    followup <- suppressWarnings(as.numeric(merged[[cols$time_followup]]))
    time <- ifelse(is.na(time), followup, time)
  }
  status <- coerce_status(merged[[cols$status]])
  ok <- !is.na(time) & !is.na(status) & !is.na(merged$score)
  if (sum(ok) < 10) {
    return(NULL)
  }

  dat <- merged[ok, , drop = FALSE]
  dat$survival_time <- time[ok]
  dat$survival_event <- status[ok]
  fit <- survival::coxph(survival::Surv(time = survival_time, event = survival_event) ~ score, data = dat)
  cox <- summary(fit)

  cutoff <- stats::median(dat$score, na.rm = TRUE)
  dat$score_group <- ifelse(dat$score >= cutoff, "High", "Low")
  km <- survival::survfit(survival::Surv(time = survival_time, event = survival_event) ~ score_group, data = dat)

  list(
    cox = data.frame(
      term = rownames(cox$coefficients),
      hazard_ratio = unname(cox$coefficients[, "exp(coef)"]),
      conf_low = if ("lower .95" %in% colnames(cox$conf.int)) unname(cox$conf.int[, "lower .95"]) else NA_real_,
      conf_high = if ("upper .95" %in% colnames(cox$conf.int)) unname(cox$conf.int[, "upper .95"]) else NA_real_,
      p_value = unname(cox$coefficients[, "Pr(>|z|)"]),
      stringsAsFactors = FALSE
    ),
    km = list(
      median_cutoff = cutoff,
      n = nrow(dat),
      groups = table(dat$score_group)
    )
  )
}

clinical_associations <- function(score_df, clinical, columns = NULL) {
  merged <- merge(score_df, clinical, by = "sample_id", all.x = TRUE)
  columns <- columns %||% c("Group", "gender", "race", "ajcc_pathologic_stage", "tumor_stage", "age_at_index")
  columns <- columns[columns %in% colnames(merged)]
  results <- list()

  for (col in columns) {
    value <- merged[[col]]
    ok <- !is.na(merged$score) & !is.na(value) & nzchar(as.character(value))
    if (sum(ok) < 5) {
      next
    }
    dat <- merged[ok, c("sample_id", "score", col), drop = FALSE]
    numeric_value <- suppressWarnings(as.numeric(as.character(dat[[col]])))
    if (is.numeric(value) || sum(!is.na(numeric_value)) >= max(5, floor(0.8 * nrow(dat)))) {
      test <- suppressWarnings(stats::cor.test(dat$score, numeric_value, method = "spearman"))
      results[[length(results) + 1L]] <- data.frame(
        variable = col,
        test = "spearman",
        estimate = unname(test$estimate),
        p_value = test$p.value,
        n = nrow(dat),
        stringsAsFactors = FALSE
      )
    } else {
      lv <- unique(as.character(dat[[col]]))
      if (length(lv) == 2L) {
        test <- suppressWarnings(stats::wilcox.test(score ~ group, data = transform(dat, group = as.factor(dat[[col]]))))
        results[[length(results) + 1L]] <- data.frame(
          variable = col,
          test = "wilcox",
          estimate = NA_real_,
          p_value = test$p.value,
          n = nrow(dat),
          stringsAsFactors = FALSE
        )
      } else if (length(lv) > 2L) {
        test <- suppressWarnings(stats::kruskal.test(score ~ group, data = transform(dat, group = as.factor(dat[[col]]))))
        results[[length(results) + 1L]] <- data.frame(
          variable = col,
          test = "kruskal",
          estimate = NA_real_,
          p_value = test$p.value,
          n = nrow(dat),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(results)) {
    return(data.frame())
  }
  do.call(rbind, results)
}

tcga_project_code <- function(project) {
  gsub("^TCGA-", "", project)
}

normalize_stage_label <- function(x) {
  if (is.na(x) || !nzchar(as.character(x))) {
    return(NA_character_)
  }
  value <- toupper(trimws(as.character(x)))
  value <- gsub("[[:space:]]+", " ", value)
  if (!grepl("^STAGE", value)) {
    return(NA_character_)
  }
  if (grepl("STAGE IV", value, fixed = TRUE)) {
    return("Stage IV")
  }
  if (grepl("STAGE III", value, fixed = TRUE)) {
    return("Stage III")
  }
  if (grepl("STAGE II", value, fixed = TRUE)) {
    return("Stage II")
  }
  if (grepl("STAGE I", value, fixed = TRUE)) {
    return("Stage I")
  }
  NA_character_
}

collapse_stage_vector <- function(x) {
  factor(
    vapply(as.character(x), normalize_stage_label, character(1)),
    levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
    ordered = TRUE
  )
}

load_pancancer_subtypes <- function() {
  path <- file.path(project_root(), "0-Data", "PanCancer_subtypes.rda")
  if (!file.exists(path)) {
    fail("Missing pan-cancer subtype file: %s", path)
  }
  obj <- load_rdata_file(path)
  if (is.list(obj) && !is.data.frame(obj) && "subtypes" %in% names(obj)) {
    obj <- obj$subtypes
  }
  x <- as.data.frame(obj, stringsAsFactors = FALSE)
  if (!("pan.samplesID" %in% colnames(x))) {
    fail("PanCancer subtype table is missing pan.samplesID")
  }
  x$patient_short <- substr(as.character(x$pan.samplesID), 1, 12)
  x
}

join_project_subtypes <- function(data, subtype_column = "Subtype_Selected") {
  subtypes <- load_pancancer_subtypes()
  code <- tcga_project_code(data$project)
  sub <- subtypes[subtypes$cancer.type %in% c(code, data$project), , drop = FALSE]
  if (!nrow(sub)) {
    fail("No pan-cancer subtype records available for %s", data$project)
  }
  sub <- sub[!duplicated(sub$patient_short), c("patient_short", subtype_column), drop = FALSE]
  merged <- merge(
    data$clinical,
    sub,
    by = "patient_short",
    all.x = TRUE
  )
  names(merged)[names(merged) == subtype_column] <- "subtype"
  merged
}

resolve_gtex_tissues <- function(project, explicit_tissues = NULL) {
  if (!is.null(explicit_tissues) && length(explicit_tissues)) {
    return(explicit_tissues)
  }
  pair_file <- file.path(project_root(), "GTEX", "samplepair.txt")
  if (!file.exists(pair_file)) {
    return(character())
  }
  pair <- utils::read.delim(pair_file, stringsAsFactors = FALSE, check.names = FALSE)
  code <- tcga_project_code(project)
  hits <- pair[pair$TCGA %in% c(code, project), , drop = FALSE]
  unique(hits$GTEx[hits$GTEx != "-"])
}

read_gtex_samples <- function() {
  path <- file.path(project_root(), "GTEX", "tissue_gtex.csv")
  if (!file.exists(path)) {
    fail("Missing GTEx tissue annotation: %s", path)
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(x)[names(x) == "Barcode"] <- "sample_id"
  names(x)[names(x) == "Tissue"] <- "tissue"
  x
}

extract_gtex_gene_vector <- function(gene_id) {
  path <- file.path(project_root(), "GTEX", "gtex_RSEM_gene_tpm.gz")
  if (!file.exists(path)) {
    fail("Missing GTEx matrix: %s", path)
  }
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  sample_ids <- strsplit(header, "\t", fixed = TRUE)[[1L]][-1L]
  gene_core <- sub("\\..*$", "", gene_id)
  cmd <- sprintf("zgrep -E -m1 '^%s(\\.[0-9]+)?[[:space:]]' %s", gene_core, shQuote(path))
  line <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
  if (!length(line)) {
    return(NULL)
  }
  fields <- strsplit(line[[1L]], "\t", fixed = TRUE)[[1L]]
  values <- suppressWarnings(as.numeric(fields[-1L]))
  values <- 2^values - 0.001
  names(values) <- sample_ids
  values
}

init_run_context <- function(config, config_path) {
  task <- config$task %||% fail("Config must include a task field.")
  task_id <- normalize_task_id(task, config$task_id %||% config$project %||% config$label)
  run_dir <- file.path(output_root(), "tcga_runs", task_id)
  overwrite <- isTRUE(getOption("tcga_toolkit.overwrite", FALSE))
  if (overwrite && dir.exists(run_dir)) {
    unlink(run_dir, recursive = TRUE)
  }
  ensure_dir(file.path(run_dir, "results"))
  ensure_dir(file.path(run_dir, "plots"))
  ensure_dir(file.path(run_dir, "objects"))
  list(
    task = task,
    task_id = task_id,
    run_dir = run_dir,
    results_dir = file.path(run_dir, "results"),
    plots_dir = file.path(run_dir, "plots"),
    objects_dir = file.path(run_dir, "objects"),
    config_path = normalizePath(config_path, mustWork = TRUE),
    toolkit_version = toolkit_version()
  )
}

write_report <- function(ctx, title, bullets) {
  lines <- c(
    sprintf("# %s", title),
    "",
    sprintf("- Task: `%s`", ctx$task),
    sprintf("- Task ID: `%s`", ctx$task_id),
    sprintf("- Generated: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    bullets
  )
  writeLines(lines, con = file.path(ctx$run_dir, "report.md"))
}

resolve_pipeline_vars <- function(obj, step_contexts) {
  if (is.list(obj) && !is.data.frame(obj)) {
    out <- obj
    for (nm in names(out)) {
      out[[nm]] <- resolve_pipeline_vars(out[[nm]], step_contexts)
    }
    return(out)
  }
  if (!is.character(obj)) {
    return(obj)
  }

  pattern <- "\\{\\{step([0-9]+)\\.([^}]+)\\}\\}"
  matches <- gregexpr(pattern, obj, perl = TRUE)

  if (matches[[1]][[1]] == -1) {
    return(obj)
  }

  result <- obj
  capture_text <- regmatches(obj, matches)[[1]]
  capture_groups <- regmatches(obj, regexec(pattern, obj, perl = TRUE))

  for (i in seq_along(capture_text)) {
    token <- capture_text[i]
    groups <- capture_groups[[i]]
    step_idx <- as.integer(groups[2]) + 1L
    key <- groups[3]

    if (step_idx < 1 || step_idx > length(step_contexts)) {
      fail("Pipeline references unknown step index: %s", groups[2])
    }

    sc <- step_contexts[[step_idx]]
    value <- NULL

    if (grepl("^glob:", key)) {
      glob_pattern <- sub("^glob:", "", key)
      candidates <- sort(Sys.glob(file.path(sc$results_dir, glob_pattern)))
      if (!length(candidates)) {
        candidates <- sort(Sys.glob(file.path(sc$objects_dir, glob_pattern)))
      }
      if (!length(candidates)) {
        fail("Pipeline glob pattern '%s' matched no files in step %s results/objects.", glob_pattern, groups[2])
      }
      value <- candidates[1]
    } else if (key %in% c("run_dir", "results_dir", "plots_dir", "objects_dir", "task_id")) {
      value <- sc[[key]]
    } else {
      fail("Unknown pipeline variable key: %s", key)
    }

    result <- sub(token, value, result, fixed = TRUE)
  }

  result
}

validate_config <- function(task, config) {
  if (!nzchar(task)) {
    fail("Config must include a task field.")
  }
  required <- switch(
    task,
    run_deg = c("project", "numerator", "denominator"),
    run_enrichment = "input_deg",
    run_clinical_assoc = if (is.null(config$gene) && is.null(config$signature_file)) c("project", "gene_or_signature_file") else "project",
    run_tme = "project",
    run_gtex_compare = c("project", "gene"),
    validate_sc_signatures = c("project", "signature_file"),
    pan_cancer_expression = if (is.null(config$gene) && is.null(config$signature_file)) "gene_or_signature_file" else NULL,
    stage_analysis = if (is.null(config$gene) && is.null(config$signature_file)) c("project", "gene_or_signature_file") else "project",
    survival_map = if (is.null(config$gene) && is.null(config$signature_file)) "gene_or_signature_file" else NULL,
    subtype_analysis = if (is.null(config$gene) && is.null(config$signature_file)) c("project", "gene_or_signature_file") else "project",
    gene_correlation_heatmap = c("target_gene", "gene_list_file"),
    pipeline = "steps",
    render_report = "run_dirs",
    maf_summary = "project",
    mutation_survival = c("project", "gene"),
    prognostic_model = {
      missing_features <- is.null(config$feature_genes) && is.null(config$signature_file) && is.null(config$gene_list_file)
      if (missing_features) c("project", "feature_genes_or_signature_or_gene_list") else "project"
    },
    ssgsea_score = c("project", "gmt_paths"),
    tmb_analysis = "project",
    cohort_qc = "project",
    cnv_summary = "project",
    methylation_diff = "project",
    external_validate = c("cohort_name", "expression_file", "clinical_file"),
    immune_phenotype = "project",
    wgcna_modules = "project",
    drug_response = "project",
    NULL
  )
  if (!is.null(required)) {
    missing <- c()
    for (field in required) {
      if (field == "gene_or_signature_file") {
        if (is.null(config$gene) && is.null(config$signature_file)) {
          missing <- c(missing, "gene or signature_file")
        }
      } else if (field == "feature_genes_or_signature_or_gene_list") {
        if (is.null(config$feature_genes) && is.null(config$signature_file) && is.null(config$gene_list_file)) {
          missing <- c(missing, "feature_genes or signature_file or gene_list_file")
        }
      } else if (is.null(config[[field]]) || (is.character(config[[field]]) && !nzchar(config[[field]]))) {
        missing <- c(missing, field)
      }
    }
    if (length(missing)) {
      fail("Task '%s' is missing required config fields: %s", task, paste(missing, collapse = ", "))
    }
  }
  invisible(TRUE)
}

project_data_cache_path <- function(project) {
  cache_dir <- file.path(output_root(), ".tcga_cache")
  ensure_dir(cache_dir)
  file.path(cache_dir, sprintf("%s_prepared.rds", project))
}

save_plot <- function(plot, basename, plots_dir, width = 8, height = 6, dpi = 150) {
  pdf_path <- file.path(plots_dir, sprintf("%s.pdf", basename))
  png_path <- file.path(plots_dir, sprintf("%s.png", basename))
  is_ggsurv <- is.list(plot) && all(c("plot", "table") %in% names(plot))
  grDevices::pdf(pdf_path, width = width, height = height)
  if (is_ggsurv) {
    print(plot, newpage = FALSE)
  } else {
    print(plot)
  }
  grDevices::dev.off()
  grDevices::png(png_path, width = width * dpi, height = height * dpi, res = dpi)
  if (is_ggsurv) {
    print(plot, newpage = FALSE)
  } else {
    print(plot)
  }
  grDevices::dev.off()
  invisible(list(pdf = pdf_path, png = png_path))
}

theme_tcga <- function(base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    fail("Plotting requires the ggplot2 package.")
  }
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}

plot_volcano <- function(deg_df, title = "Volcano Plot", padj_cutoff = 0.05, log2fc_cutoff = 1) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  df <- deg_df[, c("log2FoldChange", "padj", "gene_name")]
  df$significance <- "NS"
  df$significance[abs(df$log2FoldChange) >= log2fc_cutoff] <- "FC"
  df$significance[!is.na(df$padj) & df$padj < padj_cutoff] <- "p-value"
  df$significance[abs(df$log2FoldChange) >= log2fc_cutoff & !is.na(df$padj) & df$padj < padj_cutoff] <- "Significant"
  df$significance <- factor(df$significance, levels = c("NS", "FC", "p-value", "Significant"))
  colors <- c("NS" = "grey70", "FC" = "#00BFC4", "p-value" = "#F8766D", "Significant" = "#C77CFF")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(title = title, x = expression(log[2] ~ FoldChange), y = expression(-log[10] ~ p[adj])) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey50") +
    theme_tcga() +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  p
}

plot_boxplot_by_group <- function(df, x_col, y_col, xlab = NULL, ylab = NULL, title = NULL, fill_col = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  aes_map <- ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
  if (!is.null(fill_col) && fill_col %in% colnames(df)) {
    aes_map <- ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])
  }
  p <- ggplot2::ggplot(df, aes_map) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    ggplot2::labs(title = title, x = xlab %||% x_col, y = ylab %||% y_col) +
    theme_tcga() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (!is.null(fill_col) && fill_col %in% colnames(df)) {
    p <- p + ggplot2::theme(legend.position = "bottom")
  }
  p
}

plot_km_curve <- function(score_df, clinical, title = "Survival Curve") {
  if (!requireNamespace("survival", quietly = TRUE) || !requireNamespace("survminer", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  merged <- merge(score_df, clinical, by = "sample_id", all.x = TRUE)
  cols <- survival_columns(merged)
  if (all(is.na(unlist(cols)))) {
    return(NULL)
  }
  time <- suppressWarnings(as.numeric(merged[[cols$time_death]]))
  if (!is.na(cols$time_followup) && cols$time_followup %in% colnames(merged) && !all(is.na(merged[[cols$time_followup]]))) {
    followup <- suppressWarnings(as.numeric(merged[[cols$time_followup]]))
    time <- ifelse(is.na(time), followup, time)
  }
  status <- coerce_status(merged[[cols$status]])
  ok <- !is.na(time) & !is.na(status) & !is.na(merged$score)
  if (sum(ok) < 10) {
    return(NULL)
  }
  dat <- merged[ok, , drop = FALSE]
  dat$survival_time <- time[ok]
  dat$survival_event <- status[ok]
  cutoff <- stats::median(dat$score, na.rm = TRUE)
  dat$score_group <- ifelse(dat$score >= cutoff, "High", "Low")
  fit <- survival::survfit(survival::Surv(time = survival_time, event = survival_event) ~ score_group, data = dat)
  n_strata <- length(fit$strata)
  leg_labs <- if (n_strata == 2) c("High", "Low") else if (n_strata == 1) c("High") else NULL
  palette_vals <- if (n_strata == 2) c("#F8766D", "#00BFC4") else if (n_strata == 1) c("#F8766D") else NULL
  p <- survminer::ggsurvplot(
    fit, data = dat, pval = TRUE, risk.table = TRUE,
    title = title, xlab = "Time (days)", ylab = "Survival probability",
    palette = palette_vals,
    legend.title = "Score group", legend.labs = leg_labs,
    risk.table.y.text.col = TRUE, risk.table.y.text = FALSE
  )
  p$plot <- p$plot + theme_tcga(base_size = 10) + ggplot2::theme(legend.position = "bottom")
  p
}

plot_forest <- function(df, title = "Forest Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(df) == 0) {
    return(NULL)
  }
  df$hazard_ratio_log <- log2(df$hazard_ratio)
  df$project <- factor(df$project, levels = rev(unique(df$project)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = hazard_ratio_log, y = project)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 3, color = "#2c3e50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = hazard_ratio_log * 0.8, xmax = hazard_ratio_log * 1.2), height = 0.2, color = "#2c3e50") +
    ggplot2::labs(title = title, x = expression(log[2] ~ Hazard ~ Ratio), y = NULL) +
    theme_tcga()
  p
}

plot_dotplot_enrichment <- function(enrich_df, title = "Enrichment") {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(enrich_df) == 0) {
    return(NULL)
  }
  top_n <- min(15, nrow(enrich_df))
  df <- head(enrich_df[order(enrich_df$pvalue), ], top_n)
  if (!("Count" %in% colnames(df)) && "setSize" %in% colnames(df)) {
    df$Count <- df$setSize
  }
  if (!("Count" %in% colnames(df))) {
    df$Count <- 1
  }
  df$Description <- factor(df$Description, levels = rev(df$Description))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Count, y = Description, color = -log10(pvalue), size = Count)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "#56B1F7", high = "#132B43") +
    ggplot2::labs(title = title, x = "Gene count", y = NULL) +
    theme_tcga() +
    ggplot2::theme(legend.position = "right")
  p
}
