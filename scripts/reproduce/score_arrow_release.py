#!/usr/bin/env python3
"""Rescore the released partial-coverage ARROW fields on shared QM probes."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np

from audit_arrow_package import bootstrap_paired, read_xyz

ROOT = Path(__file__).resolve().parents[2]
SETS = ("panel_1", "panel_2", "panel_3", "liquid", "shell_size")
METHODS = ("arrow_smeared", "arrow", "glider", "mace_polar_l")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def metric(prediction: Path, reference: Path) -> float:
    with np.load(prediction, allow_pickle=False) as p, np.load(reference, allow_pickle=False) as q:
        np.testing.assert_array_equal(p["points_angstrom"], q["points_angstrom"])
        delta = p["predicted_esp_hartree_per_e"] - q["delta_esp_hartree_per_e"]
        return float(np.sqrt(np.mean(delta * delta) / np.mean(q["delta_esp_hartree_per_e"] ** 2)))


def main() -> None:
    expected = json.loads((ROOT / "experiments/arrow_comparison/audit.json").read_text())
    scores = defaultdict(list)
    cases = set()
    for set_name in SETS:
        base = ROOT / "experiments" / set_name
        for method in ("arrow_smeared", "arrow"):
            prediction_dir = base / "predictions" / method
            report = {
                r["config_id"]: r for r in rows(
                    ROOT / "experiments/arrow_comparison/source_tables" / set_name /
                    f"{method}_results.csv"
                )
            }
            for row in rows(prediction_dir / "prediction_registry.csv"):
                config = row["config_id"]
                file = row["prediction_file"]
                prediction = prediction_dir / file
                if hashlib.sha256(prediction.read_bytes()).hexdigest() != row["prediction_sha256"]:
                    raise AssertionError(f"Hash mismatch: {prediction}")
                value = metric(prediction, base / "references" / file)
                if abs(value - float(report[config]["esp_nrmse"])) > 1e-10:
                    raise AssertionError(f"Archived score mismatch: {set_name}/{config}/{method}")
                scores[(set_name, row["molecule_id"], method)].append(value)
                if method == "arrow_smeared":
                    cases.add((set_name, config, row["molecule_id"], file))
    for set_name, config, molecule, file in cases:
        base = ROOT / "experiments" / set_name
        for method in ("glider", "mace_polar_l"):
            scores[(set_name, molecule, method)].append(
                metric(base / "predictions" / method / file, base / "references" / file)
            )
    means = {key: statistics.mean(values) for key, values in scores.items()}
    for set_name in SETS:
        solutes = sorted({mol for s, mol, _ in means if s == set_name})
        for method in METHODS:
            value = statistics.mean(means[(set_name, mol, method)] for mol in solutes)
            target = expected["sets"][set_name]["nrmse"][method]
            if abs(value - target) > 1e-10:
                raise AssertionError(f"Set summary mismatch: {set_name}/{method}")
    for label, selected in {
        "all_68": SETS,
        "without_shell_44": SETS[:-1],
        "without_panel1_or_shell_40": ("panel_2", "panel_3", "liquid"),
    }.items():
        groups = sorted({(s, mol) for s, mol, _ in means if s in selected})
        for method in METHODS:
            value = statistics.mean(means[(s, mol, method)] for s, mol in groups)
            target = expected["aggregates"][label]["nrmse"][method]
            if abs(value - target) > 1e-10:
                raise AssertionError(f"Aggregate mismatch: {label}/{method}")
        paired = [
            (s, mol, means[(s, mol, "arrow_smeared")] - means[(s, mol, "glider")])
            for s, mol in groups
        ]
        difference, lower, upper = bootstrap_paired(paired)
        comparison = expected["aggregates"][label]
        if abs(difference - comparison["arrow_minus_glider"]) > 1e-10:
            raise AssertionError(f"Paired mean mismatch: {label}")
        for observed, target in zip(
            (lower, upper), comparison["solute_cluster_bootstrap_95_percentile_interval"]
        ):
            if abs(observed - target) > 1e-10:
                raise AssertionError(f"Bootstrap interval mismatch: {label}")

    geometry = read_xyz(ROOT / "experiments/shell_size/geometries/configurations.extxyz")
    shell_files = {config: file for set_name, config, _, file in cases if set_name == "shell_size"}
    bias_count = 0
    for row in rows(ROOT / "experiments/arrow_comparison/shell_size_signed_bias_per_config.csv"):
        config = row["config_id"]
        method = row["method"]
        if config not in shell_files or method not in METHODS:
            continue
        file = shell_files[config]
        base = ROOT / "experiments/shell_size"
        positions, n_solute_atoms = geometry[config]
        with np.load(base / "references" / file, allow_pickle=False) as reference:
            points = reference["points_angstrom"]
            qm = reference["delta_esp_hartree_per_e"]
        with np.load(base / "predictions" / method / file, allow_pickle=False) as prediction:
            np.testing.assert_array_equal(prediction["points_angstrom"], points)
            predicted = prediction["predicted_esp_hartree_per_e"]
        nearest = np.argmin(
            np.linalg.norm(points[:, None, :] - positions[None, :, :], axis=2), axis=1
        )
        signed_error = float(np.mean((predicted - qm)[nearest < n_solute_atoms]) * 1000)
        if abs(signed_error - float(row["signed_error_solute_surface_mEh"])) > 1e-6:
            raise AssertionError(f"Shell-size signed-bias mismatch: {config}/{method}")
        bias_count += 1
    if bias_count != expected["signed_bias_rows_verified"]:
        raise AssertionError(f"Expected 96 signed-bias rows, found {bias_count}")
    print(
        f"PASS: {len(cases)} common configurations, {bias_count} signed-bias rows; "
        "ARROW hashes, probe alignment, NRMSE and paired intervals reproduced"
    )


if __name__ == "__main__":
    main()
