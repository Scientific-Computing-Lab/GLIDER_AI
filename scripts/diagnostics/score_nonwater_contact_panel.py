#!/usr/bin/env python3
"""Score the frozen contact-geometry follow-up after all 36 QM labels exist."""

from __future__ import annotations

import csv
import hashlib
import json
import argparse
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
METHODS = ("glider", "mace_polar_l")
BOOTSTRAPS = 100_000
SEED = 20261007


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, data: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def mean(data: list[float]) -> float:
    return float(np.mean(data))


def main(base: Path) -> None:
    base = base.resolve()
    if not base.is_dir():
        raise FileNotFoundError(base)
    freeze = json.loads((base / "prediction_freeze.json").read_text())
    if freeze["geometry_sha256"] != sha256(base / "configurations.extxyz"):
        raise ValueError("Frozen contact geometries changed")
    reference_path = base / "references"
    manifest = json.loads((reference_path / "REFERENCE_MANIFEST.json").read_text())
    if manifest["n_configurations"] != 36:
        raise ValueError("Expected 36 completed QM references")
    references = {row["config_id"]: row for row in rows(reference_path / "REFERENCE_REGISTRY.csv")}
    if len(references) != 36 or sha256(reference_path / "REFERENCE_REGISTRY.csv") != manifest["reference_registry_sha256"]:
        raise ValueError("Reference registry incomplete or changed")
    geometry = {row["config_id"]: row for row in rows(base / "configuration_registry.csv")}
    if set(geometry) != set(references):
        raise ValueError("Geometry/reference case mismatch")
    predictions: dict[str, dict[str, dict[str, str]]] = {}
    for method in METHODS:
        prediction_registry = base / "predictions" / method / "prediction_registry.csv"
        frozen = freeze[f"{method}_predictions"]
        if sha256(prediction_registry) != frozen["registry_sha256"]:
            raise ValueError(f"Frozen prediction registry changed for {method}")
        data = rows(prediction_registry)
        predictions[method] = {row["config_id"]: row for row in data}
        if len(data) != 36 or set(predictions[method]) != set(references):
            raise ValueError(f"Prediction coverage mismatch for {method}")
    scores: list[dict[str, object]] = []
    for config_id in sorted(references):
        ref_row = references[config_id]
        ref_file = reference_path / ref_row["observable_file"]
        if sha256(ref_file) != ref_row["observable_sha256"]:
            raise ValueError(f"Reference digest mismatch for {config_id}")
        with np.load(ref_file, allow_pickle=False) as archive:
            points = archive["points_angstrom"]
            ref_esp = archive["delta_esp_hartree_per_e"]
            ref_dipole = archive["delta_dipole_debye"]
        for method in METHODS:
            prediction_row = predictions[method][config_id]
            prediction_file = base / "predictions" / method / prediction_row["prediction_file"]
            if sha256(prediction_file) != prediction_row["prediction_sha256"]:
                raise ValueError(f"Prediction digest mismatch for {config_id}, {method}")
            if prediction_row["prediction_sha256"] != freeze[f"{method}_predictions"]["case_file_sha256"][config_id]:
                raise ValueError(f"Prediction changed after freeze: {config_id}, {method}")
            with np.load(prediction_file, allow_pickle=False) as archive:
                np.testing.assert_array_equal(points, archive["points_angstrom"])
                pred_esp = archive["predicted_esp_hartree_per_e"]
                pred_dipole = archive["predicted_dipole_debye"]
            if pred_esp.shape != ref_esp.shape:
                raise ValueError(f"ESP shape mismatch for {config_id}, {method}")
            difference = pred_esp - ref_esp
            dipole_difference = pred_dipole - ref_dipole
            scores.append({
                "config_id": config_id,
                "molecule_id": geometry[config_id]["molecule_id"],
                "perturbant": geometry[config_id]["perturbant"],
                "method": method,
                "esp_nrmse": float(np.sqrt(np.mean(difference**2) / np.mean(ref_esp**2))),
                "esp_rmse_mEh_per_e": float(np.sqrt(np.mean(difference**2)) * 1000),
                "response_rms_mEh_per_e": float(np.sqrt(np.mean(ref_esp**2)) * 1000),
                "dipole_mse_D2": float(np.mean(dipole_difference**2)),
                "cp_interaction_energy_kcal_mol": float(ref_row["cp_interaction_energy_kcal_mol"]),
            })
    write_rows(base / "results.csv", scores)
    molecules = sorted({str(row["molecule_id"]) for row in scores})
    methods = {method: [row for row in scores if row["method"] == method] for method in METHODS}
    per_solute: list[dict[str, object]] = []
    for method in METHODS:
        for molecule_id in molecules:
            subset = [row for row in methods[method] if row["molecule_id"] == molecule_id]
            per_solute.append({
                "molecule_id": molecule_id,
                "method": method,
                "esp_nrmse": mean([float(row["esp_nrmse"]) for row in subset]),
                "esp_rmse_mEh_per_e": mean([float(row["esp_rmse_mEh_per_e"]) for row in subset]),
                "dipole_mse_D2": mean([float(row["dipole_mse_D2"]) for row in subset]),
            })
    write_rows(base / "solute_results.csv", per_solute)
    average = {
        method: {
            "esp_nrmse": mean([float(row["esp_nrmse"]) for row in per_solute if row["method"] == method]),
            "esp_rmse_mEh_per_e": mean([float(row["esp_rmse_mEh_per_e"]) for row in per_solute if row["method"] == method]),
            "dipole_rmse_D": float(np.sqrt(mean([float(row["dipole_mse_D2"]) for row in per_solute if row["method"] == method]))),
        }
        for method in METHODS
    }
    groups = ("overall", "NH3", "CH3OH", "CH3CN")
    paired_by_group: dict[str, np.ndarray] = {}
    group_scores: dict[str, dict[str, float]] = {}
    for group in groups:
        glider_values = []
        baseline_values = []
        for molecule in molecules:
            for method, values in (("glider", glider_values), ("mace_polar_l", baseline_values)):
                subset = [
                    float(row["esp_nrmse"]) for row in methods[method]
                    if row["molecule_id"] == molecule
                    and (group == "overall" or row["perturbant"] == group)
                ]
                if len(subset) != (3 if group == "overall" else 1):
                    raise ValueError(f"Missing {group} case for {molecule}, {method}")
                values.append(mean(subset))
        glider_array = np.asarray(glider_values)
        baseline_array = np.asarray(baseline_values)
        paired_by_group[group] = glider_array - baseline_array
        group_scores[group] = {
            "glider": float(np.mean(glider_array)),
            "mace_polar_l": float(np.mean(baseline_array)),
        }
    paired = paired_by_group["overall"]
    generator = np.random.default_rng(SEED)
    indices = generator.integers(0, len(molecules), size=(BOOTSTRAPS, len(molecules)))
    resamples = paired[indices].mean(axis=1)
    group_intervals = {
        group: np.quantile(values[indices].mean(axis=1), [0.025, 0.975]).tolist()
        for group, values in paired_by_group.items()
    }
    effects = []
    for label, y in (("overall", 3), ("NH3", 2), ("CH3OH", 1), ("CH3CN", 0)):
        point = float(np.mean(paired_by_group[label]))
        low, high = group_intervals[label]
        effects.append({
            "label": label,
            "y": y,
            "esp_difference": point,
            "err_minus": point - low,
            "err_plus": high - point,
            "n_solutes": len(molecules),
            "n_configurations": 36 if label == "overall" else 12,
        })
    write_rows(base / "figure_effects.csv", effects)
    cp = np.array([float(row["cp_interaction_energy_kcal_mol"]) for row in methods["glider"]])
    summary = {
        "status": "new QM labels scored after geometry and prediction freeze; separate post-critique follow-up",
        "n_configurations": len(references),
        "n_solutes": len(molecules),
        "n_neighbour_species": 3,
        "equal_solute_results": average,
        "equal_solute_esp_nrmse_by_neighbour": group_scores,
        "glider_minus_mace_l_esp_nrmse": float(np.mean(paired)),
        "paired_95pct_solute_bootstrap_interval": np.quantile(resamples, [0.025, 0.975]).tolist(),
        "paired_95pct_solute_bootstrap_interval_by_neighbour": group_intervals,
        "paired_bootstrap_seed": SEED,
        "n_bootstrap": BOOTSTRAPS,
        "solute_wins_glider_vs_mace_l": int(np.sum(paired < 0)),
        "case_wins_glider_vs_mace_l": sum(
            float(next(row for row in methods["glider"] if row["config_id"] == config_id)["esp_nrmse"])
            < float(next(row for row in methods["mace_polar_l"] if row["config_id"] == config_id)["esp_nrmse"])
            for config_id in references
        ),
        "cp_interaction_energy_kcal_mol": {
            "median": float(np.median(cp)),
            "minimum": float(np.min(cp)),
            "maximum": float(np.max(cp)),
            "positive_count": int(np.sum(cp > 0)),
            "above_10_count": int(np.sum(cp > 10)),
        },
        "prediction_freeze_sha256": sha256(base / "prediction_freeze.json"),
        "reference_manifest_sha256": sha256(reference_path / "REFERENCE_MANIFEST.json"),
        "result_sha256": sha256(base / "results.csv"),
        "figure_effects_sha256": sha256(base / "figure_effects.csv"),
    }
    (base / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "build/nonwater_contact_panel")
    main(parser.parse_args().base)
