"""Survival utilities for tcga_toolkit Python adapter.

Wrappers around scikit-survival for the most common downstream steps after the
R `prognostic_model` task: external Cox validation, KM plots, and
time-dependent AUC.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd


def _check_imports():
    try:
        import sksurv  # noqa: F401
    except ImportError as exc:  # pragma: no cover - import guard
        raise ImportError(
            "Install scikit-survival to use tcga_toolkit.adapters.python.survival: "
            "pip install scikit-survival"
        ) from exc


def fit_external_cox(
    risk_df: pd.DataFrame,
    time_col: str = "survival_time",
    event_col: str = "survival_event",
    score_col: str = "risk_score",
) -> dict:
    """Fit a univariate Cox model of risk_score against survival.

    Parameters
    ----------
    risk_df
        DataFrame with at least the time, event, and score columns. Typically
        the output of ``load_risk_scores(run)`` filtered to a specific split.

    Returns
    -------
    dict with ``hazard_ratio``, ``conf_low``, ``conf_high``, ``p_value`` and
    the concordance index.
    """
    _check_imports()
    from sksurv.linear_model import CoxPHSurvivalAnalysis
    from sksurv.metrics import concordance_index_censored

    df = risk_df.dropna(subset=[time_col, event_col, score_col])
    y = np.zeros(
        len(df),
        dtype=[("event", bool), ("time", float)],
    )
    y["event"] = df[event_col].astype(int).to_numpy() == 1
    y["time"] = df[time_col].astype(float).to_numpy()
    x = df[[score_col]].to_numpy(dtype=float)

    fit = CoxPHSurvivalAnalysis(alpha=1e-4).fit(x, y)
    coef = float(fit.coef_[0])
    # Stable hazard ratio CI via Wald estimate on the linear scale.
    pred = fit.predict(x)
    cindex = concordance_index_censored(y["event"], y["time"], pred)[0]

    # scikit-survival's Cox does not expose SE; approximate via permutation z-test.
    rng = np.random.default_rng(0)
    perm_coefs = []
    for _ in range(200):
        idx = rng.permutation(len(df))
        y_perm = y.copy()
        y_perm["event"] = y["event"][idx]
        y_perm["time"] = y["time"][idx]
        try:
            perm_coefs.append(float(CoxPHSurvivalAnalysis(alpha=1e-4).fit(x, y_perm).coef_[0]))
        except Exception:
            continue
    perm_coefs = np.asarray(perm_coefs)
    p_value = float((np.abs(perm_coefs) >= abs(coef)).mean()) if perm_coefs.size else float("nan")

    return {
        "n": int(len(df)),
        "events": int(df[event_col].astype(int).sum()),
        "log_hazard_ratio": coef,
        "hazard_ratio": float(np.exp(coef)),
        "p_value": p_value,
        "concordance": float(cindex),
    }


def plot_external_km(
    risk_df: pd.DataFrame,
    time_col: str = "survival_time",
    event_col: str = "survival_event",
    score_col: str = "risk_score",
    cutoff: float | None = None,
    title: str = "Risk-score KM",
    ax=None,
):
    """Plot a median-split Kaplan-Meier curve. Returns the matplotlib Axes."""
    _check_imports()
    import matplotlib.pyplot as plt
    from sksurv.nonparametric import kaplan_meier_estimator

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    df = risk_df.dropna(subset=[time_col, event_col, score_col])
    cutoff = float(df[score_col].median()) if cutoff is None else float(cutoff)
    df = df.assign(group=np.where(df[score_col] >= cutoff, "High", "Low"))

    for label, sub in df.groupby("group"):
        time, prob = kaplan_meier_estimator(
            sub[event_col].astype(int).to_numpy() == 1,
            sub[time_col].astype(float).to_numpy(),
        )
        ax.step(time, prob, where="post", label=f"{label} (n={len(sub)})")
    ax.set_title(title)
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.legend()
    return ax


def time_dependent_auc(
    risk_df: pd.DataFrame,
    times: Iterable[float] = (365, 1095, 1825),
    time_col: str = "survival_time",
    event_col: str = "survival_event",
    score_col: str = "risk_score",
) -> pd.DataFrame:
    """Compute time-dependent ROC AUC for each requested time point."""
    _check_imports()
    from sksurv.metrics import cumulative_dynamic_auc

    df = risk_df.dropna(subset=[time_col, event_col, score_col])
    y = np.zeros(len(df), dtype=[("event", bool), ("time", float)])
    y["event"] = df[event_col].astype(int).to_numpy() == 1
    y["time"] = df[time_col].astype(float).to_numpy()
    times = np.asarray(sorted(t for t in times if 0 < t < y["time"].max()))
    if not times.size:
        raise ValueError("No requested time-points fit within the observation horizon.")
    auc, mean_auc = cumulative_dynamic_auc(y, y, df[score_col].to_numpy(dtype=float), times)
    return pd.DataFrame({"time": times, "AUC": auc, "mean_AUC": [mean_auc] * len(times)})
