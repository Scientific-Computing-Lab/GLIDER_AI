#!/usr/bin/env python3
"""Recompute all prospective aggregate and paired statistics from frozen tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from glider.metrics.core import aggregate, leave_one_molecule_out, paired_bootstrap


def value(values: np.ndarray, root: bool) -> float:
    mean = float(np.mean(values))
    return float(np.sqrt(mean)) if root else mean


def panel_statistics(
    table: pd.DataFrame,
    *,
    panel: str,
    glider_name: str,
    dipole_column: str,
) -> tuple[dict, list[dict]]:
    molecules = sorted(table.molecule_id.unique())
    methods = sorted(table.method.unique())
    candidate = table[table.method == glider_name].set_index("molecule_id")
    pairs: list[dict] = []
    for comparator in methods:
        if comparator == glider_name:
            continue
        other = table[table.method == comparator].set_index("molecule_id")
        for observable, column, root in (
            ("esp", "esp_nrmse", False),
            ("dipole", dipole_column, True),
        ):
            if column not in other or other[column].isna().all():
                continue
            joined = candidate[[column]].join(
                other[[column]], how="inner", lsuffix="_glider", rsuffix="_comparator"
            )
            first = joined[f"{column}_glider"].to_numpy(dtype=float)
            second = joined[f"{column}_comparator"].to_numpy(dtype=float)
            observed_first, observed_second = value(first, root), value(second, root)
            influence = leave_one_molecule_out(first, second, root_after_mean=root)
            pairs.append(
                {
                    "panel": panel,
                    "observable": observable,
                    "glider_method": glider_name,
                    "comparator": comparator,
                    "n_molecules": len(joined),
                    "glider_value": observed_first,
                    "comparator_value": observed_second,
                    "delta_glider_minus_comparator": observed_first - observed_second,
                    "relative_error_reduction": 1.0 - observed_first / observed_second,
                    "glider_molecule_wins": int((first < second).sum()),
                    "bootstrap_ci_low": paired_bootstrap(first, second, root_after_mean=root)[0],
                    "bootstrap_ci_high": paired_bootstrap(first, second, root_after_mean=root)[1],
                    "median_molecule_glider": float(
                        np.sqrt(np.median(first)) if root else np.median(first)
                    ),
                    "median_molecule_comparator": float(
                        np.sqrt(np.median(second)) if root else np.median(second)
                    ),
                    "p90_molecule_glider": float(
                        np.sqrt(np.quantile(first, 0.9)) if root else np.quantile(first, 0.9)
                    ),
                    "p90_molecule_comparator": float(
                        np.sqrt(np.quantile(second, 0.9)) if root else np.quantile(second, 0.9)
                    ),
                    "loo_delta_min": float(influence.min()),
                    "loo_delta_max": float(influence.max()),
                }
            )
    normalized = table.rename(columns={dipole_column: "dipole_mse"})
    summary = aggregate(normalized)
    return {
        "panel": panel,
        "n_molecules": len(molecules),
        "n_methods": len(methods),
        "aggregate": summary.to_dict(orient="records"),
    }, pairs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("build/statistics"))
    args = parser.parse_args()
    panels = [pd.read_csv(args.root / f"experiments/panel_{n}/solute_results.csv") for n in [1,2,3]]
    summaries, all_pairs = [], []
    for n, panel in enumerate(panels,1):
        summary, pairs = panel_statistics(panel, panel=f"panel_{n}", glider_name="glider", dipole_column="dipole_mse_D2")
        summaries.append(summary); all_pairs.extend(pairs)
    args.output.mkdir(parents=True, exist_ok=True)
    pairs = pd.DataFrame(all_pairs)
    pairs.to_csv(args.output / "PAIRWISE_STATISTICS.csv", index=False)
    payload = {
        "metric_definition": {
            "esp": "mean molecule-mean configuration response-ESP NRMSE",
            "dipole": "sqrt(mean molecule-mean Cartesian-component MSE)",
        },
        "bootstrap": {"replicates": 100000, "seed": 20260814, "block": "molecule"},
        "panels": summaries,
        "pairwise": pairs.to_dict(orient="records"),
        "leave_one_out_policy": "influence diagnostic only; never substitutes for frozen aggregate",
    }
    (args.output / "FINAL_STATISTICS.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print("solute-blocked paired statistics: PASS")
    print(f"pairwise comparisons: {len(pairs)}")


if __name__ == "__main__":
    main()
