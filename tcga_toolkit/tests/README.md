# Smoke tests

Lightweight tests that verify the toolkit's wiring without running real
analyses. They take seconds, not minutes, and are safe to run any time.

## Run them all

```bash
bash tcga_toolkit/tests/run_all.sh
```

## Individually

```bash
Rscript tcga_toolkit/tests/test_introspection.R
Rscript tcga_toolkit/tests/test_templates.R
Rscript tcga_toolkit/tests/test_dispatch.R
python3 tcga_toolkit/tests/test_python_adapter.py
```

## What each test checks

| Test | Verifies |
| --- | --- |
| `test_introspection.R` | `list_tasks.R`, `inspect_project.R`, `validate_config.R`, `list_runs.R`, `check_deps.R` all execute and return the expected shape. |
| `test_templates.R` | Every JSON in `templates/` passes `validate_config.R` (schema-level; missing local data files are tolerated). |
| `test_dispatch.R` | `common.R` + every `task_*.R` source without error and every advertised task has a `task_<name>` function in scope. |
| `test_python_adapter.py` | `tcga_toolkit/adapters/python` imports cleanly and exports the public API. (Does not require sksurv / shap / lightgbm.) |

When adding a new task, also bump the assertions in
[`test_dispatch.R`](test_dispatch.R) — the test pulls the advertised list
from `list_tasks.R` automatically, so the only thing to do is make sure
your new task is added there.
