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


def verify_headlines(panel_1: pd.DataFrame, panel_2: pd.DataFrame, panel_3: pd.DataFrame) -> None:
    first = aggregate(panel_1).set_index("method")
    second = aggregate(panel_2.rename(columns={"dipole_mse_D2": "dipole_mse"})).set_index("method")
    third = aggregate(panel_3.rename(columns={"dipole_mse_D2": "dipole_mse"})).set_index("method")
    expected = {
        ("p1", "candidate", "esp_nrmse"): 0.32179935443942215,
        ("p1", "candidate", "dipole_rmse_debye"): 0.18935929367382542,
        ("p2", "previous_candidate", "esp_nrmse"): 0.34217487331895485,
        ("p2", "previous_candidate", "dipole_rmse_debye"): 0.11254123272430284,
        ("p2", "candidate", "esp_nrmse"): 0.363855548302602,
        ("p2", "candidate", "dipole_rmse_debye"): 0.25021615552684345,
        ("p3", "glider", "esp_nrmse"): 0.42658060581223267,
        ("p3", "glider", "dipole_rmse_debye"): 0.22669745000951386,
    }
    for (panel, method, metric), target in expected.items():
        observed = {"p1": first, "p2": second, "p3": third}[panel].loc[method, metric]
        if not np.isclose(observed, target, atol=1.0e-14, rtol=0):
            raise AssertionError((panel, method, metric, observed, target))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("evidence"))
    args = parser.parse_args()
    panel_1 = pd.read_csv(args.root / "results/prospective_1/FINAL_MOLECULE_COMPARISONS.csv")
    panel_2 = pd.read_csv(args.root / "results/prospective_2/NEXTGEN_FINAL_MOLECULE_RESULTS.csv")
    panel_3 = pd.read_csv(args.root / "canonical_panel_3/MOLECULE_LEVEL_METRICS.csv")
    panel_3["dipole_mse_D2"] = panel_3.dipole_rmse_D**2
    verify_headlines(panel_1, panel_2, panel_3)
    summary_1, pairs_1 = panel_statistics(
        panel_1, panel="prospective_1", glider_name="candidate", dipole_column="dipole_mse"
    )
    summary_2, pairs_2 = panel_statistics(
        panel_2,
        panel="prospective_2",
        glider_name="previous_candidate",
        dipole_column="dipole_mse_D2",
    )
    summary_3, pairs_3 = panel_statistics(
        panel_3,
        panel="prospective_3_canonical_freesolv_identities",
        glider_name="glider",
        dipole_column="dipole_mse_D2",
    )
    args.output.mkdir(parents=True, exist_ok=True)
    pairs = pd.DataFrame(pairs_1 + pairs_2 + pairs_3)
    pairs.to_csv(args.output / "PAIRWISE_STATISTICS.csv", index=False)
    payload = {
        "metric_definition": {
            "esp": "mean molecule-mean configuration response-ESP NRMSE",
            "dipole": "sqrt(mean molecule-mean Cartesian-component MSE)",
        },
        "bootstrap": {"replicates": 100000, "seed": 20260814, "block": "molecule"},
        "panels": [summary_1, summary_2, summary_3],
        "pairwise": pairs.to_dict(orient="records"),
        "leave_one_out_policy": "influence diagnostic only; never substitutes for frozen aggregate",
    }
    (args.output / "FINAL_STATISTICS.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print("headline verification: PASS")
    print(f"pairwise comparisons: {len(pairs)}")


if __name__ == "__main__":
    main()
