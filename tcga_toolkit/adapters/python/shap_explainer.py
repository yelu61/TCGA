"""SHAP explainability for prognostic models produced by the R toolkit.

The R `prognostic_model` task exports raw expression features and risk
scores so this side can re-fit a Python surrogate model (gradient boosted
trees) and explain it with SHAP. This sidesteps cross-language model object
loading entirely.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _require(name: str, install: str):
    try:
        return __import__(name)
    except ImportError as exc:  # pragma: no cover
        raise ImportError(f"Install '{name}' to use this function: {install}") from exc


def _prepare_surrogate(
    risk_df: pd.DataFrame,
    feature_matrix: pd.DataFrame,
    time_col: str,
    event_col: str,
    score_col: str,
):
    merged = feature_matrix.merge(
        risk_df[["sample_id", score_col, time_col, event_col]], on="sample_id"
    )
    merged = merged.dropna(subset=[score_col, time_col, event_col])
    y_target = merged[score_col].to_numpy(dtype=float)
    feature_cols = [c for c in feature_matrix.columns if c != "sample_id"]
    x = merged[feature_cols].to_numpy(dtype=float)
    return x, y_target, feature_cols, merged


def explain_prognostic_model_rsf(
    risk_df: pd.DataFrame,
    feature_matrix: pd.DataFrame,
    *,
    time_col: str = "survival_time",
    event_col: str = "survival_event",
    score_col: str = "risk_score",
    top_features: int = 20,
):
    """Explain an RSF risk score via a gradient-boosted surrogate + TreeSHAP.

    Parameters
    ----------
    risk_df
        DataFrame with sample_id, risk_score, survival_time, survival_event.
        Typically ``load_risk_scores(run)`` from the R toolkit.
    feature_matrix
        DataFrame with sample_id plus one column per feature. You produce this
        in R via ``data$tpm`` for the resolved features and write it out as
        a CSV; load with ``pd.read_csv``.

    Returns
    -------
    dict with ``shap_values`` (np.ndarray), ``feature_names``, ``top``
    (DataFrame of mean(|SHAP|) per feature), and ``model``.
    """
    shap = _require("shap", "pip install shap")
    lgb = _require("lightgbm", "pip install lightgbm")

    x, y, feature_cols, _ = _prepare_surrogate(
        risk_df, feature_matrix, time_col, event_col, score_col
    )
    model = lgb.LGBMRegressor(n_estimators=400, learning_rate=0.05, max_depth=6)
    model.fit(x, y)

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(x)

    importance = np.abs(shap_values).mean(axis=0)
    top_df = (
        pd.DataFrame({"feature": feature_cols, "mean_abs_shap": importance})
        .sort_values("mean_abs_shap", ascending=False)
        .head(top_features)
        .reset_index(drop=True)
    )

    return {
        "shap_values": shap_values,
        "feature_names": feature_cols,
        "top": top_df,
        "model": model,
        "explainer": explainer,
    }


def explain_prognostic_model_lasso(
    coefficients_csv: str | "pd.DataFrame",
    top_features: int = 20,
) -> pd.DataFrame:
    """For Lasso-Cox there is no need for SHAP — the coefficients *are* the
    explanation. This helper just sorts them by absolute magnitude and
    returns the top features. Accepts either the path to
    ``<project>_prognostic_coefficients.csv`` or an already-loaded DataFrame.
    """
    if isinstance(coefficients_csv, str):
        df = pd.read_csv(coefficients_csv)
    else:
        df = coefficients_csv.copy()
    df["abs_coef"] = df["coefficient"].abs()
    return df.sort_values("abs_coef", ascending=False).head(top_features).reset_index(drop=True)
