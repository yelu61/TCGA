# Extension Guide

How to add a new task to the toolkit without breaking existing runs.

A "task" is a single config-driven analysis dispatched by
[`run_task.R`](../scripts/run_task.R). Adding one touches four files plus
optional supporting material.

---

## Step 1 — Create `task_<name>.R`

Place a new file under `tcga_toolkit/scripts/task_<name>.R`. The function
must be named `task_<name>` and accept exactly two arguments:

```r
task_<name> <- function(config, ctx) {
  # config: parsed list from the JSON/YAML config file
  # ctx: list with task, task_id, run_dir, results_dir, plots_dir, objects_dir,
  #      config_path, toolkit_version

  project <- config$project %||% fail("<name> requires a project.")
  data <- standardize_project_data(project)

  # ... write CSVs to ctx$results_dir, plots via save_plot(...), rds to ctx$objects_dir

  write_report(
    ctx,
    title = sprintf("My task for %s", project),
    bullets = c(
      sprintf("- Samples: `%s`", ncol(data$counts)),
      "",
      "## Outputs",
      sprintf("- `results/%s_something.csv`", project)
    )
  )

  invisible(list(...))
}
```

Conventions to follow (consistent with existing tasks):

- Use `%||%` for default values and `fail(...)` for hard errors.
- Use `requireNamespace("pkg", quietly = TRUE)` before calling optional
  packages, and call `fail()` with a clear message if missing.
- Write all CSVs through `write_table_safe()`.
- Save figures with `save_plot()` so you get both PDF and PNG.
- Save R objects with `saveRDS()` into `ctx$objects_dir`.
- Always call `write_report(ctx, title, bullets)` so the run produces a
  `report.md` that `render_report` can pick up.

## Step 2 — Register the dispatcher

In [`scripts/run_task.R`](../scripts/run_task.R) add:

```r
source(file.path(script_dir, "task_<name>.R"))
```

and add an entry to the `dispatch` list:

```r
dispatch <- list(
  ...,
  <name> = task_<name>
)
```

## Step 3 — Update the schema

In [`scripts/common.R`](../scripts/common.R) extend `validate_config()`
with required fields:

```r
required <- switch(
  task,
  ...,
  <name> = c("project", "field_a"),
  NULL
)
```

If the requirement is "one of A or B", reuse the `gene_or_signature_file`
pattern.

## Step 4 — Document

- Add a row to [`tcga_toolkit/scripts/list_tasks.R`](../scripts/list_tasks.R)
  describing required/optional fields. Keep the summary one line.
- Add a section to [`task_index.md`](task_index.md) with the same
  schema, plus output filenames and gotchas.
- Append a row to [`tcga_toolkit/specs/README.md`](../specs/README.md).
- Add an example config under `tcga_toolkit/templates/`.
- Add a recipe to [`pipeline_recipes.md`](pipeline_recipes.md) if it
  composes naturally with existing tasks.
- Add a CHANGELOG entry in [`CHANGELOG.md`](../CHANGELOG.md).
- Update `check_deps.R` if your task introduces new R packages.

## Step 5 — Smoke test

Validate that the new task is wired in correctly without running the
analysis:

```bash
Rscript tcga_toolkit/scripts/list_tasks.R --task <name>
Rscript tcga_toolkit/scripts/check_deps.R --task <name>
Rscript tcga_toolkit/scripts/validate_config.R --config tcga_toolkit/templates/<name>_example.json
```

Then run a tiny end-to-end test against a small project (e.g.
`TCGA-CHOL`, `TARGET-OS`) to ensure outputs and the report markdown render
correctly.

## Style and stability

- Don't break existing config files. If you change defaults, bump the
  toolkit minor version in [`VERSION`](../VERSION) and note it in the
  CHANGELOG.
- Avoid hard-coding paths outside the helpers in `common.R`; everything
  should go through `project_root()`, `project_files()` or
  `output_root()`.
- Prefer extending `common.R` (e.g. a new helper such as
  `compute_signature_score()`) over duplicating logic across tasks.
- If your task needs new clinical fields, add them to the
  `preferred_cols` list inside `standardize_project_data()` so they are
  available to all tasks.
