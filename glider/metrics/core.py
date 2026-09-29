"""Frozen metric definitions used by both prospective GLIDER campaigns."""

from __future__ import annotations

import numpy as np
import pandas as pd


def aggregate(
    molecule_rows: pd.DataFrame,
    *,
    esp_column: str = "esp_nrmse",
    dipole_mse_column: str = "dipole_mse",
) -> pd.DataFrame:
    """Aggregate with equal molecule weight, matching the frozen protocols."""
    result = molecule_rows.groupby("method", as_index=False).agg(
        esp_nrmse=(esp_column, "mean"), dipole_mse=(dipole_mse_column, "mean")
    )
    result["dipole_rmse_debye"] = np.sqrt(result.pop("dipole_mse"))
    return result


def paired_bootstrap(
    candidate: np.ndarray,
    comparator: np.ndarray,
    *,
    root_after_mean: bool,
    replicas: int = 100_000,
    seed: int = 20260814,
) -> tuple[float, float]:
    """Percentile CI for candidate-minus-comparator with molecule blocks."""
    candidate = np.asarray(candidate, dtype=float)
    comparator = np.asarray(comparator, dtype=float)
    if candidate.shape != comparator.shape or candidate.ndim != 1:
        raise ValueError("paired arrays must be one-dimensional and equal length")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(candidate), size=(replicas, len(candidate)))
    first = candidate[indices].mean(axis=1)
    second = comparator[indices].mean(axis=1)
    if root_after_mean:
        first, second = np.sqrt(first), np.sqrt(second)
    delta = first - second
    return float(np.quantile(delta, 0.025)), float(np.quantile(delta, 0.975))


def leave_one_molecule_out(
    candidate: np.ndarray, comparator: np.ndarray, *, root_after_mean: bool
) -> np.ndarray:
    """Diagnostic influence values; never used to replace the frozen aggregate."""
    values = []
    for held in range(len(candidate)):
        mask = np.arange(len(candidate)) != held
        first, second = candidate[mask].mean(), comparator[mask].mean()
        if root_after_mean:
            first, second = np.sqrt(first), np.sqrt(second)
        values.append(first - second)
    return np.asarray(values)
