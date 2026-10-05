#!/usr/bin/env python3
"""Apply the immutable conjunctive prospective SOTA gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from response_learning import sha256


def load_prediction(path: Path):
    table = pd.read_csv(path / "prediction_registry.csv")
    return {
        row.config_id: path / row.prediction_file for row in table.itertuples()
    }, table


def aggregate(frame: pd.DataFrame):
    molecule = frame.groupby(["method", "molecule_id"], as_index=False).agg(
        esp_nrmse=("esp_nrmse", "mean"),
        dipole_mse=("dipole_mse", "mean"),
    )
    summary = molecule.groupby("method", as_index=False).agg(
        esp_nrmse=("esp_nrmse", "mean"), dipole_mse=("dipole_mse", "mean")
    )
    summary["dipole_vector_rmse_debye"] = np.sqrt(summary.dipole_mse)
    return molecule, summary.drop(columns="dipole_mse")


def bootstrap(molecule, candidate, comparator, metric, replicas=100000):
    pivot = molecule.pivot(index="molecule_id", columns="method", values=metric)
    if metric == "dipole_mse":

        def value(x):
            return np.sqrt(np.mean(x, axis=1))
    else:

        def value(x):
            return np.mean(x, axis=1)

    delta = np.column_stack([pivot[candidate].to_numpy(), pivot[comparator].to_numpy()])
    rng = np.random.default_rng(20260813)
    indices = rng.integers(0, len(delta), size=(replicas, len(delta)))
    difference = value(delta[indices, 0]) - value(delta[indices, 1])
    return float(np.quantile(difference, 0.025)), float(np.quantile(difference, 0.975))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--prediction", action="append", required=True, help="NAME=DIR")
    parser.add_argument("--candidate-runtime-seconds", type=float, required=True)
    parser.add_argument("--rotation-audit", type=Path, required=True)
    parser.add_argument("--freeze-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reference_table = pd.read_csv(args.reference / "observable_registry.csv")
    references = {
        row.config_id: args.reference / row.observable_file
        for row in reference_table.itertuples()
    }
    prediction_sets = {}
    registries = {}
    for item in args.prediction:
        name, raw = item.split("=", 1)
        prediction_sets[name], registries[name] = load_prediction(Path(raw))
    if "candidate" not in prediction_sets:
        raise RuntimeError("A prediction named candidate is required")
    expected = set(references)
    if any(set(values) != expected for values in prediction_sets.values()):
        raise RuntimeError("Prediction/reference configuration sets differ")
    rows = []
    for method, predictions in prediction_sets.items():
        for config_id in sorted(expected):
            ref = np.load(references[config_id])
            pred = np.load(predictions[config_id])
            if not np.allclose(
                ref["points_angstrom"], pred["points_angstrom"], atol=1e-10, rtol=0
            ):
                raise RuntimeError(f"ESP grid mismatch: {method}/{config_id}")
            target, estimate = (
                ref["delta_esp_hartree_per_e"],
                pred["predicted_esp_hartree_per_e"],
            )
            dip_error = pred["predicted_dipole_debye"] - ref["delta_dipole_debye"]
            meta = reference_table.loc[reference_table.config_id == config_id].iloc[0]
            rows.append(
                {
                    "method": method,
                    "config_id": config_id,
                    "molecule_id": meta.molecule_id,
                    "regime": meta.regime,
                    "esp_nrmse": float(
                        np.sqrt(np.mean((estimate - target) ** 2))
                        / max(np.sqrt(np.mean(target**2)), 1e-12)
                    ),
                    "dipole_mse": float(np.mean(dip_error**2)),
                    "dipole_vector_rmse_debye": float(np.sqrt(np.mean(dip_error**2))),
                    "net_charge_e": float(np.sum(pred["predicted_charges_e"])),
                }
            )
    frame = pd.DataFrame(rows)
    molecule, summary = aggregate(frame)
    baselines = [x for x in prediction_sets if x != "candidate"]
    esp_comparator = (
        summary[summary.method.isin(baselines)].sort_values("esp_nrmse").iloc[0].method
    )
    dip_comparator = (
        summary[summary.method.isin(baselines)]
        .sort_values("dipole_vector_rmse_debye")
        .iloc[0]
        .method
    )
    lookup = summary.set_index("method")
    candidate = lookup.loc["candidate"]
    esp_ci = bootstrap(molecule, "candidate", esp_comparator, "esp_nrmse")
    dip_ci = bootstrap(molecule, "candidate", dip_comparator, "dipole_mse")
    mol = molecule.pivot(index="molecule_id", columns="method")
    improved = int(
        np.sum(
            (mol["esp_nrmse"]["candidate"] < mol["esp_nrmse"][esp_comparator])
            & (mol["dipole_mse"]["candidate"] < mol["dipole_mse"][dip_comparator])
        )
    )

    def regime_ok(regime):
        part = frame[frame.regime == regime]
        c_esp = part[part.method == "candidate"].esp_nrmse.mean()
        b_esp = part[part.method == esp_comparator].esp_nrmse.mean()
        c_dip = np.sqrt(part[part.method == "candidate"].dipole_mse.mean())
        b_dip = np.sqrt(part[part.method == dip_comparator].dipole_mse.mean())
        return bool(c_esp <= b_esp and c_dip <= b_dip)

    rotation = json.loads(args.rotation_audit.read_text())
    qm_runtime = float(reference_table.runtime_seconds.sum())
    speedup = qm_runtime / args.candidate_runtime_seconds
    esp_base = lookup.loc[esp_comparator, "esp_nrmse"]
    dip_base = lookup.loc[dip_comparator, "dipole_vector_rmse_debye"]
    gates = {
        "lowest_esp": bool(candidate.esp_nrmse < lookup.loc[baselines].esp_nrmse.min()),
        "lowest_dipole": bool(
            candidate.dipole_vector_rmse_debye
            < lookup.loc[baselines].dipole_vector_rmse_debye.min()
        ),
        "esp_improvement_at_least_20pct": bool(candidate.esp_nrmse <= 0.8 * esp_base),
        "dipole_improvement_at_least_15pct": bool(
            candidate.dipole_vector_rmse_debye <= 0.85 * dip_base
        ),
        "esp_bootstrap_ci_favors": bool(esp_ci[1] < 0),
        "dipole_bootstrap_ci_favors": bool(dip_ci[1] < 0),
        "at_least_9_of_12_improve_both": bool(improved >= 9),
        "orientation_no_regression_both": regime_ok("orientation"),
        "compressed_no_regression_both": regime_ok("compressed"),
        "absolute_esp_nrmse_le_30pct": bool(candidate.esp_nrmse <= 0.30),
        "charge_conservation": bool(
            frame[frame.method == "candidate"].net_charge_e.abs().max() <= 1e-6
        ),
        "rotation_stability": bool(rotation["pass"]),
        "speedup_at_least_100x": bool(speedup >= 100),
    }
    passed = all(gates.values())
    args.output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output / "configuration_metrics.csv", index=False)
    molecule.to_csv(args.output / "molecule_metrics.csv", index=False)
    summary.to_csv(args.output / "comparator_table.csv", index=False)
    result = {
        "sota_gate_passed": passed,
        "gates": gates,
        "candidate": candidate.to_dict(),
        "esp_comparator": esp_comparator,
        "dipole_comparator": dip_comparator,
        "esp_relative_improvement": float(1 - candidate.esp_nrmse / esp_base),
        "dipole_relative_improvement": float(
            1 - candidate.dipole_vector_rmse_debye / dip_base
        ),
        "esp_paired_bootstrap_delta_95pct_ci": esp_ci,
        "dipole_paired_bootstrap_delta_95pct_ci": dip_ci,
        "molecules_improved_on_both": improved,
        "qm_runtime_seconds": qm_runtime,
        "candidate_runtime_seconds": args.candidate_runtime_seconds,
        "speedup": speedup,
        "freeze_manifest_sha256": sha256(args.freeze_manifest),
        "experimental_hydration_targets_accessed": False,
    }
    (args.output / "gate_result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
