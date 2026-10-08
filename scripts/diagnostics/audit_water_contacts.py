#!/usr/bin/env python3
"""Audit short solute--water contacts in the three frozen water panels.

This is a post hoc geometry sensitivity analysis. It never edits the archived
geometries, predictions, reference observables, or primary panel scores.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "experiments" / "water_contact_audit"
THRESHOLD_ANGSTROM = 1.5
BOOTSTRAP_DRAWS = 100_000
PANELS = (
    ("I", "panel_1", "historical_results/final_configuration_metrics.csv", "candidate"),
    ("II", "panel_2", "results.csv", "glider"),
    ("III", "panel_3", "CONFIGURATION_LEVEL_METRICS.csv", "glider"),
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, data: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def closest_pairs(path: Path) -> dict[str, tuple[float, str, str]]:
    result = {}
    lines = path.read_text().splitlines()
    offset = 0
    while offset < len(lines):
        n_atoms = int(lines[offset])
        header = lines[offset + 1]
        config_id = header.split("config_id=", 1)[1].split()[0]
        n_solute = int(header.split("n_solute_atoms=", 1)[1].split()[0])
        atoms = []
        for line in lines[offset + 2 : offset + 2 + n_atoms]:
            fields = line.split()
            atoms.append((fields[0], tuple(map(float, fields[1:4]))))
        best = (float("inf"), "", "")
        for solute in atoms[:n_solute]:
            for water in atoms[n_solute:]:
                distance = math.dist(solute[1], water[1])
                if distance < best[0]:
                    best = (distance, solute[0], water[0])
        if config_id in result:
            raise ValueError(f"Repeated geometry: {config_id}")
        result[config_id] = best
        offset += n_atoms + 2
    return result


def closest_water_pairs(path: Path) -> dict[str, tuple[float, str, str, int, int]]:
    result = {}
    lines = path.read_text().splitlines()
    offset = 0
    while offset < len(lines):
        n_atoms = int(lines[offset])
        header = lines[offset + 1]
        config_id = header.split("config_id=", 1)[1].split()[0]
        atoms = []
        for line in lines[offset + 2 : offset + 2 + n_atoms]:
            fields = line.split()
            atoms.append((fields[0], tuple(map(float, fields[1:4])), int(fields[4])))
        best = (float("inf"), "", "", -1, -1)
        for first in atoms:
            for second in atoms:
                if first[2] < 1 or second[2] <= first[2]:
                    continue
                distance = math.dist(first[1], second[1])
                if distance < best[0]:
                    best = (distance, first[0], second[0], first[2], second[2])
        result[config_id] = best
        offset += n_atoms + 2
    return result


def summarize(
    panel: str,
    geometry: dict[str, dict[str, str]],
    metrics: dict[tuple[str, str], float],
    glider_method: str,
    exclude_short: bool,
    extra_exclusions: set[str] | None = None,
) -> dict[str, object]:
    by_solute: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {"glider": [], "baseline": []}
    )
    included = 0
    for config_id, geometry_row in geometry.items():
        short = (
            float(geometry_row["minimum_solute_water_distance_A"]) < THRESHOLD_ANGSTROM
            or config_id in (extra_exclusions or set())
        )
        if exclude_short and short:
            continue
        by_solute[geometry_row["molecule_id"]]["glider"].append(
            metrics[(config_id, glider_method)]
        )
        by_solute[geometry_row["molecule_id"]]["baseline"].append(
            metrics[(config_id, "mace_polar_l")]
        )
        included += 1
    if not by_solute or any(not values["glider"] for values in by_solute.values()):
        raise ValueError("Every solute must retain a configuration")
    solutes = sorted(by_solute)
    glider = np.array([np.mean(by_solute[s]["glider"]) for s in solutes])
    baseline = np.array([np.mean(by_solute[s]["baseline"]) for s in solutes])
    differences = glider - baseline
    rng = np.random.default_rng(20261008 + (1 if exclude_short else 0) + len(solutes))
    resampled = rng.integers(0, len(solutes), size=(BOOTSTRAP_DRAWS, len(solutes)))
    ci_low, ci_high = np.percentile(differences[resampled].mean(axis=1), [2.5, 97.5])
    glider_mean = float(glider.mean())
    baseline_mean = float(baseline.mean())
    return {
        "panel": panel,
        "analysis": (
            "exclude_any_interfragment_contact_below_1.5_A"
            if extra_exclusions else
            "exclude_solute_water_contacts_below_1.5_A"
            if exclude_short else "frozen_primary"
        ),
        "n_solutes": len(solutes),
        "n_configurations": included,
        "glider_equal_solute_nrmse": glider_mean,
        "mace_l_equal_solute_nrmse": baseline_mean,
        "reduction_pct": 100 * (1 - glider_mean / baseline_mean),
        "glider_minus_mace_l": float(differences.mean()),
        "paired_95pct_ci_low": float(ci_low),
        "paired_95pct_ci_high": float(ci_high),
        "solute_wins": int(np.sum(differences < 0)),
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    short_rows = []
    short_water_rows = []
    summaries = []
    total = 0
    input_hashes = {}
    for panel, directory, metric_name, glider_method in PANELS:
        folder = ROOT / "experiments" / directory
        for relative in (
            "geometries/configuration_registry.csv",
            "geometries/configurations.extxyz",
            "references/observable_registry.csv",
            metric_name,
        ):
            path = folder / relative
            input_hashes[str(path.relative_to(ROOT))] = sha256(path)
        geometry_rows = rows(folder / "geometries" / "configuration_registry.csv")
        geometry = {row["config_id"]: row for row in geometry_rows}
        if len(geometry) != len(geometry_rows):
            raise ValueError(f"Duplicate geometry IDs in Panel {panel}")
        pair = closest_pairs(folder / "geometries" / "configurations.extxyz")
        water_pair = closest_water_pairs(folder / "geometries" / "configurations.extxyz")
        if set(pair) != set(geometry):
            raise ValueError(f"Geometry registry mismatch in Panel {panel}")
        for config_id, row in geometry.items():
            reported = float(row["minimum_solute_water_distance_A"])
            if abs(reported - pair[config_id][0]) > 1e-6:
                raise ValueError(f"Distance mismatch: {config_id}")
        metric_rows = rows(folder / metric_name)
        reference_rows = {
            row["config_id"]: row
            for row in rows(folder / "references" / "observable_registry.csv")
        }
        if set(reference_rows) != set(geometry):
            raise ValueError(f"Reference registry mismatch in Panel {panel}")
        metrics = {
            (row["config_id"], row["method"]): float(row["esp_nrmse"])
            for row in metric_rows
            if row["esp_nrmse"]
        }
        if len(metrics) != sum(bool(row["esp_nrmse"]) for row in metric_rows):
            raise ValueError(f"Duplicate metrics in Panel {panel}")
        for row in geometry_rows:
            config_id = row["config_id"]
            distance, solute_element, water_element = pair[config_id]
            water_distance, water_element_a, water_element_b, water_a, water_b = water_pair[config_id]
            if water_distance < THRESHOLD_ANGSTROM:
                short_water_rows.append({
                    "panel": panel,
                    "config_id": config_id,
                    "solute_id": row["molecule_id"],
                    "regime": row["regime"],
                    "nearest_water_element_a": water_element_a,
                    "nearest_water_element_b": water_element_b,
                    "water_fragment_a": water_a,
                    "water_fragment_b": water_b,
                    "minimum_water_water_distance_A": water_distance,
                    "glider_esp_nrmse": metrics[(config_id, glider_method)],
                    "mace_l_esp_nrmse": metrics[(config_id, "mace_polar_l")],
                })
            if distance >= THRESHOLD_ANGSTROM:
                continue
            short_rows.append({
                "panel": panel,
                "config_id": config_id,
                "solute_id": row["molecule_id"],
                "regime": row["regime"],
                "nearest_solute_element": solute_element,
                "nearest_water_element": water_element,
                "minimum_distance_A": distance,
                "qm_response_rms_mEh_per_e":
                    1000 * float(reference_rows[config_id]["delta_esp_rms_hartree_per_e"]),
                "glider_esp_nrmse": metrics[(config_id, glider_method)],
                "mace_l_esp_nrmse": metrics[(config_id, "mace_polar_l")],
                "glider_esp_rmse_mEh_per_e":
                    1000 * float(reference_rows[config_id]["delta_esp_rms_hartree_per_e"])
                    * metrics[(config_id, glider_method)],
                "mace_l_esp_rmse_mEh_per_e":
                    1000 * float(reference_rows[config_id]["delta_esp_rms_hartree_per_e"])
                    * metrics[(config_id, "mace_polar_l")],
            })
        summaries.append(summarize(panel, geometry, metrics, glider_method, False))
        summaries.append(summarize(panel, geometry, metrics, glider_method, True))
        extra = {row["config_id"] for row in short_water_rows if row["panel"] == panel}
        if extra:
            summaries.append(summarize(panel, geometry, metrics, glider_method, True, extra))
        total += len(geometry)
    if total != 224 or len(short_rows) != 8 or len(short_water_rows) != 1:
        raise ValueError(
            f"Unexpected panel/contact counts: {total}, {len(short_rows)}, {len(short_water_rows)}"
        )
    write_csv(OUTPUT / "short_contacts.csv", short_rows)
    write_csv(OUTPUT / "short_water_water_contacts.csv", short_water_rows)
    write_csv(OUTPUT / "panel_sensitivity.csv", summaries)
    (OUTPUT / "analysis.json").write_text(json.dumps({
        "scope": "Three frozen water panels, post hoc geometry sensitivity",
        "short_contact_rule": "minimum interfragment atomic distance < 1.5 angstrom",
        "n_original_configurations": total,
        "n_short_solute_water_contacts": len(short_rows),
        "n_short_water_water_contacts": len(short_water_rows),
        "n_unique_configurations_with_a_short_interfragment_contact": len(
            {row["config_id"] for row in short_rows + short_water_rows}
        ),
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "primary_data_changed": False,
        "interaction_energies_available_in_released_water_references": False,
        "source_sha256": input_hashes,
    }, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
