# TCGA toolkit — Python adapter

Lightweight helpers that read the output of the R `tcga_toolkit` tasks and
make them analysable with the Python stack (pandas, scikit-survival,
SHAP). It is not a Python re-implementation of the toolkit — only a
*reader* and a few downstream conveniences.

## Install dependencies

```bash
pip install pandas numpy matplotlib
pip install scikit-survival          # for survival.py
pip install lightgbm shap            # for shap_explainer.py
```

## Quickstart

```python
import sys, os
sys.path.insert(0, os.path.expanduser("tcga_toolkit/adapters"))
from python import (
    available_runs, load_run, load_run_table, load_risk_scores,
    fit_external_cox, plot_external_km, time_dependent_auc,
    explain_prognostic_model_lasso,
)

# 1. discover prior runs
runs = available_runs(".", task="prognostic_model", status="success")

# 2. open the latest one
run = load_run(runs[-1])
print(run.task_id, run.toolkit_version, run.status)

# 3. read its outputs as pandas DataFrames
risk = load_risk_scores(run)
coefs = load_run_table(run, "*prognostic_coefficients.csv")

# 4. validate the risk score externally
result = fit_external_cox(risk[risk["split"] == "test"])
print(result)

# 5. time-dependent AUC
auc = time_dependent_auc(risk[risk["split"] == "test"], times=[365, 1095])
print(auc)

# 6. lasso coefficients sorted by importance
top = explain_prognostic_model_lasso(coefs)
print(top.head(10))
```

## When to use which file

| Module | Purpose |
| --- | --- |
| `runs.py` | Discover and load run directories produced by `Rscript run_task.R`. |
| `survival.py` | scikit-survival wrappers — Cox, KM, time-dependent AUC. |
| `shap_explainer.py` | Surrogate-model SHAP for RSF risk scores; helper for Lasso. |

## Design notes

- The adapter avoids loading any R objects (`.rds`, `.rda`) — it reads the
  CSV / JSON / PNG artefacts the R tasks write to disk. This keeps the
  Python side dependency-free of `rpy2`.
- The SHAP module deliberately does not try to read an `RSF.rds` directly.
  Instead it trains a Python surrogate (LightGBM) on the exported risk
  scores so SHAP works on a native Python model.
- For Lasso-Cox there is no need for SHAP — the L1 coefficients *are*
  the explanation. `explain_prognostic_model_lasso` simply sorts them.

## File layout

```
tcga_toolkit/adapters/python/
├── __init__.py          # re-exports the public API as `python` package
├── runs.py              # tcga_runs/ reader
├── survival.py          # scikit-survival helpers
├── shap_explainer.py    # SHAP for RSF, sorter for Lasso
└── README.md            # this file
```
