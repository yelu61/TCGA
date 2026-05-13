"""Read tcga_runs/<task_id>/ outputs in Python."""

from __future__ import annotations

import json
from dataclasses import dataclass
from glob import glob
from pathlib import Path
from typing import Iterable

import pandas as pd


@dataclass
class Run:
    """Lightweight handle to a tcga_runs/<task_id>/ directory."""

    run_dir: Path
    metadata: dict
    results: dict[str, Path]
    plots: dict[str, Path]
    objects: dict[str, Path]
    report_md: Path | None

    @property
    def task(self) -> str:
        return self.metadata.get("task", "unknown")

    @property
    def task_id(self) -> str:
        return self.metadata.get("task_id", self.run_dir.name)

    @property
    def toolkit_version(self) -> str:
        return self.metadata.get("toolkit_version", "0.0.0")

    @property
    def status(self) -> str:
        return self.metadata.get("status", "unknown")


def available_runs(
    output_root: str | Path = ".",
    task: str | None = None,
    status: str | None = None,
) -> list[Path]:
    """Return run directories under ``output_root/tcga_runs/``.

    Pass the *project* root (the directory that contains ``tcga_toolkit/``);
    this helper appends ``tcga_runs``.
    """
    root = Path(output_root) / "tcga_runs"
    if not root.is_dir():
        return []
    runs: list[Path] = []
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        meta_path = d / "run_metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if task is not None and meta.get("task") != task:
            continue
        if status is not None and meta.get("status") != status:
            continue
        runs.append(d)
    return runs


def load_run(run_dir: str | Path) -> Run:
    """Wrap a single run directory in a :class:`Run` dataclass."""
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    meta_path = run_dir / "run_metadata.json"
    metadata = json.loads(meta_path.read_text()) if meta_path.exists() else {}

    def collect(subdir: str, patterns: Iterable[str]) -> dict[str, Path]:
        out: dict[str, Path] = {}
        d = run_dir / subdir
        if not d.is_dir():
            return out
        for pat in patterns:
            for p in sorted(d.glob(pat)):
                out[p.name] = p
        return out

    return Run(
        run_dir=run_dir,
        metadata=metadata,
        results=collect("results", ("*.csv", "*.json")),
        plots=collect("plots", ("*.png", "*.pdf")),
        objects=collect("objects", ("*.rds", "*.rda")),
        report_md=(run_dir / "report.md") if (run_dir / "report.md").exists() else None,
    )


def load_run_table(run: Run, name: str) -> pd.DataFrame:
    """Load a CSV from ``results/`` by filename or glob (e.g. ``*deg.csv``)."""
    matches: list[Path] = []
    if name in run.results:
        matches = [run.results[name]]
    else:
        matches = [p for n, p in run.results.items() if Path(n).match(name)]
    if not matches:
        raise FileNotFoundError(f"No table matching {name!r} in {run.run_dir / 'results'}")
    return pd.read_csv(matches[0])


def load_risk_scores(run: Run) -> pd.DataFrame:
    """Load risk-score output written by ``prognostic_model``.

    Returns a DataFrame with columns at least ``sample_id``, ``split``,
    ``risk_score``, ``survival_time``, ``survival_event``.
    """
    candidates = [n for n in run.results if n.endswith("_prognostic_risk_scores.csv")]
    if not candidates:
        raise FileNotFoundError(
            f"No prognostic_risk_scores.csv found in {run.run_dir / 'results'}"
        )
    return pd.read_csv(run.results[candidates[0]])
