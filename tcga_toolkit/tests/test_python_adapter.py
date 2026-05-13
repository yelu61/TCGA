"""Smoke test for the Python adapter — only checks that the package imports
cleanly without requiring scikit-survival / SHAP / LightGBM. Run with::

    python tcga_toolkit/tests/test_python_adapter.py
"""

from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path


def main() -> int:
    here = Path(__file__).resolve().parent
    adapter_root = here.parent / "adapters"
    sys.path.insert(0, str(adapter_root))

    passes = 0
    failures: list[str] = []

    def check(name: str, ok: bool, detail: str = ""):
        nonlocal passes
        if ok:
            passes += 1
            print(f"  PASS  {name}")
        else:
            failures.append(name)
            print(f"  FAIL  {name}")
            if detail:
                print(f"        {detail}")

    # 1. `python` package exists (adapters/python/__init__.py).
    check("adapters/python directory exists", (adapter_root / "python").is_dir())

    # 2. importable.
    try:
        mod = importlib.import_module("python")
    except Exception as exc:  # noqa: BLE001
        check("import python adapter", False, repr(exc))
        return 1
    check("import python adapter", True)

    # 3. Public API.
    expected = {
        "available_runs",
        "load_run",
        "load_run_table",
        "load_risk_scores",
        "fit_external_cox",
        "plot_external_km",
        "time_dependent_auc",
        "explain_prognostic_model_rsf",
        "explain_prognostic_model_lasso",
    }
    missing = expected - set(dir(mod))
    check("public API surface", not missing, f"missing: {missing}")

    # 4. available_runs handles a non-existent root gracefully.
    res = mod.available_runs("/nonexistent/path")
    check("available_runs handles missing root", res == [])

    print(f"\nPython adapter: {passes} passed, {len(failures)} failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
