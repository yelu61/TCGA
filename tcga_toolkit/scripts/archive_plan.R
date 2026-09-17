# Validate the read-only plan emitted by manage-research-project's lifecycle
# helper. Candidate paths authorize archiving; unreferenced rows are review only.
validate_archive_plan <- function(plan_path, project) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Archive-plan verification requires the digest R package.", call. = FALSE)
  plan <- jsonlite::read_json(plan_path, simplifyVector = FALSE)
  project <- normalizePath(project, mustWork = TRUE)
  if (is.null(plan$project_root) || !identical(plan$project_root, project)) {
    stop("Archive plan belongs to a different project root.", call. = FALSE)
  }
  relative_path <- function(path) {
    if (!is.character(path) || length(path) != 1L || !nzchar(path) || startsWith(path, "/") ||
        grepl("(^|/)(\\.\\.?)(/|$)|[\\\\\t\r\n]", path)) {
      stop("Unsafe project-relative path in archive plan.", call. = FALSE)
    }
    parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
    current <- project
    for (part in parts) {
      current <- file.path(current, part)
      link <- Sys.readlink(current)
      if (!is.na(link) && nzchar(link)) stop("Archive plan paths cannot traverse symbolic links.", call. = FALSE)
    }
    current
  }
  freshness <- plan$freshness
  hashes <- freshness$sha256
  scopes <- freshness$scope_files
  if (is.null(hashes) || is.null(names(hashes)) || is.null(scopes) || !length(scopes) ||
      !all(c("provenance/RUNS.tsv", "provenance/ARTIFACTS.tsv") %in% names(hashes))) {
    stop("Archive plan lacks registry and result-scope freshness evidence.", call. = FALSE)
  }
  for (relative in names(hashes)) {
    path <- relative_path(relative)
    if (!file.exists(path) || dir.exists(path) ||
        !identical(digest::digest(file = path, algo = "sha256", serialize = FALSE), hashes[[relative]])) {
      stop("Stale archive plan: changed or missing file ", relative, call. = FALSE)
    }
  }
  for (relative in names(scopes)) {
    scope <- relative_path(relative)
    if (file.exists(scope) && !dir.exists(scope)) stop("Stale archive plan: result scope is not a directory.", call. = FALSE)
    walk <- function(directory) {
      children <- list.files(directory, recursive = FALSE, full.names = TRUE,
                              all.files = TRUE, no.. = TRUE, include.dirs = TRUE)
      links <- Sys.readlink(children)
      if (any(!is.na(links) & nzchar(links))) stop("Result scope contains a symbolic link; archive plan is not verifiable.", call. = FALSE)
      c(children, unlist(lapply(children[file.info(children)$isdir], walk), use.names = FALSE))
    }
    entries <- if (dir.exists(scope)) walk(scope) else character()
    files <- entries[!file.info(entries)$isdir]
    actual <- sort(substring(files, nchar(project) + 2L))
    expected <- sort(as.character(unlist(scopes[[relative]], use.names = FALSE)))
    if (!identical(actual, expected) || !all(actual %in% names(hashes))) {
      stop("Stale archive plan: result-scope files changed in ", relative, call. = FALSE)
    }
  }
  candidates <- as.character(unlist(plan$candidate_native_directories %||% character(), use.names = FALSE))
  protected <- as.character(unlist(plan$protected_native_directories %||% character(), use.names = FALSE))
  if (any(candidates %in% protected)) stop("Archive plan protects a candidate directory.", call. = FALSE)
  for (path in unique(c(candidates, protected))) relative_path(path)
  list(plan = plan, candidates = unique(candidates), protected = unique(protected))
}
