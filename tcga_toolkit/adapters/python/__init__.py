"""TCGA toolkit Python adapter — read tcga_runs/ outputs in Python.

Lightweight helpers that load a run directory produced by the R toolkit so the
results can be analysed with the broader Python data-science stack
(pandas / scikit-survival / SHAP / matplotlib).
"""

from .runs import (
    available_runs,
    load_run,
    load_run_table,
    load_risk_scores,
)
from .survival import (
    fit_external_cox,
    plot_external_km,
    time_dependent_auc,
)
from .shap_explainer import (
    explain_prognostic_model_rsf,
    explain_prognostic_model_lasso,
)

__all__ = [
    "available_runs",
    "load_run",
    "load_run_table",
    "load_risk_scores",
    "fit_external_cox",
    "plot_external_km",
    "time_dependent_auc",
    "explain_prognostic_model_rsf",
    "explain_prognostic_model_lasso",
]
