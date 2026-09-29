"""Release-registry integrity checks."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def validate_registry(path: str | Path, *, molecules: int, configurations: int) -> None:
    table = pd.read_csv(path)
    required = {"config_id", "molecule_id", "regime"}
    if not required.issubset(table.columns):
        raise ValueError(f"registry lacks {sorted(required - set(table.columns))}")
    if len(table) != configurations or table.molecule_id.nunique() != molecules:
        raise ValueError("unexpected prospective panel dimensions")
    if table.config_id.duplicated().any():
        raise ValueError("duplicate configuration identifier")
    expected_regimes = {"equilibrium", "cooperative", "orientation", "compressed"}
    if set(table.regime) != expected_regimes:
        raise ValueError("unexpected regime set")
    counts = table.groupby("molecule_id").regime.nunique()
    if not (counts == 4).all():
        raise ValueError("each molecule must contain all four regimes")
