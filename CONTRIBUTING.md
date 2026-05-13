# Contributing to tcga-toolkit

Thanks for considering a contribution. The toolkit is small and config-driven,
so most contributions take one of three shapes:

1. **A new task** — adds an analysis capability.
2. **A bug fix** — corrects a task's behaviour on existing data.
3. **Documentation** — improves `references/`, `specs/`, README, or templates.

## Before you start

- Run the smoke tests locally: `bash tcga_toolkit/tests/run_all.sh`. They
  should pass before and after your change.
- Open an issue first for substantive task additions or schema-breaking
  changes so we can agree on the design.
- Bump `tcga_toolkit/VERSION` only on coordinated releases.

## Adding a new task

The full procedure is in
[`tcga_toolkit/references/extension_guide.md`](tcga_toolkit/references/extension_guide.md).
The checklist:

- [ ] `tcga_toolkit/scripts/task_<name>.R` follows existing style:
      `function(config, ctx)`, uses `%||%`, `fail(...)`, `write_table_safe()`,
      `save_plot()`, ends with `write_report(ctx, ...)`.
- [ ] Registered in `run_task.R` (source + dispatch entry).
- [ ] Required fields in `common.R::validate_config()`.
- [ ] One-line summary in `scripts/list_tasks.R`.
- [ ] Dependency declaration in `scripts/check_deps.R`.
- [ ] Long-form documentation in `references/task_index.md`.
- [ ] Schema row in `specs/README.md`.
- [ ] Example template in `templates/`.
- [ ] CHANGELOG entry.

## Style

- 2-space indentation.
- `<-` for assignment.
- snake_case for functions and config fields.
- No comments that restate the obvious; do explain hidden constraints.
- Never introduce hard-coded paths outside the helpers in `common.R`.

## Testing

Smoke tests live under `tcga_toolkit/tests/`. They verify wiring (every task
sources, every template passes schema validation, Python adapter imports).
They are intentionally **not** end-to-end — that requires data the user
provides — but they should always be green.

If your change is non-trivial, please also run an end-to-end check against
a small project (e.g. `TCGA-CHOL`, `TARGET-OS`) and paste the
`run_metadata.json` snippet into the PR description.

## Data hygiene

Do not commit:

- `tcga_runs/` outputs, `*.rda`, `*.rds`, large CSVs, raw GDC downloads.
- Anything that identifies a specific patient cohort beyond the standard
  TCGA / TARGET / GTEx sample IDs.
- Personal API keys, login tokens, local paths in scripts.

The `.gitignore` already covers the common cases.

## License

By contributing you agree that your contributions are licensed under the
project's [Apache License 2.0](LICENSE).
