task_pipeline <- function(config, ctx) {
  steps <- config$steps %||% fail("pipeline requires a steps array.")
  if (!length(steps)) {
    fail("pipeline steps array is empty.")
  }

  step_contexts <- list()
  step_results <- list()

  for (i in seq_along(steps)) {
    step_raw <- steps[[i]]
    if (is.null(step_raw$task)) {
      fail("Pipeline step %s is missing a task field.", i - 1)
    }

    step_config <- resolve_pipeline_vars(step_raw, step_contexts)
    step_config$task <- step_raw$task

    step_task_id <- step_config$task_id %||% sprintf("%s-step%s", ctx$task_id, i - 1)
    step_config$task_id <- step_task_id

    info("Pipeline step %s/%s: %s (task_id=%s)", i, length(steps), step_config$task, step_task_id)

    step_ctx <- init_run_context(step_config, ctx$config_path)

    handler <- dispatch[[step_config$task]]
    if (is.null(handler)) {
      fail("Unknown pipeline task: %s", step_config$task)
    }

    result <- NULL
    error_message <- NULL
    step_status <- "success"

    tryCatch({
      result <- handler(step_config, step_ctx)
    }, error = function(e) {
      error_message <<- conditionMessage(e)
      step_status <<- "failed"
      message(sprintf("Pipeline step %s failed: %s", i, error_message))
    })

    write_json(
      list(
        task = step_config$task,
        task_id = step_ctx$task_id,
        toolkit_version = step_ctx$toolkit_version,
        config_path = step_ctx$config_path,
        run_dir = step_ctx$run_dir,
        status = step_status,
        error = error_message,
        completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ),
      file.path(step_ctx$run_dir, "run_metadata.json")
    )

    step_contexts[[length(step_contexts) + 1L]] <- step_ctx
    step_results[[length(step_results) + 1L]] <- list(
      status = step_status,
      run_dir = step_ctx$run_dir,
      error = error_message
    )

    if (step_status == "failed") {
      fail("Pipeline aborted at step %s/%s (%s): %s", i, length(steps), step_config$task, error_message)
    }
  }

  summary_lines <- c(
    sprintf("- Steps executed: `%s`", length(steps)),
    "",
    "## Step Details"
  )
  for (i in seq_along(step_contexts)) {
    sc <- step_contexts[[i]]
    summary_lines <- c(
      summary_lines,
      sprintf("- Step %s: `%s` → `%s`", i, steps[[i]]$task, sc$run_dir)
    )
  }

  write_report(
    ctx,
    sprintf("Pipeline: %s", ctx$task_id),
    summary_lines
  )

  list(step_contexts = step_contexts, step_results = step_results)
}
