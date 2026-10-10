#!/usr/bin/env python3
"""Independently rescore the collaborator's partial-coverage ARROW archive.

Usage: python scripts/reproduce/audit_arrow_package.py /path/to/ARROW_results.zip
The archive is read without extraction. No ARROW code or parameters are run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import random
import statistics
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PREFIX = "ARROW_results_for_GLIDER/for_gal_arrow_results/"
SETS = ("panel_1", "panel_2", "panel_3", "liquid", "shell_size")
METHODS = ("arrow_smeared", "arrow", "glider", "mace_polar_l")


def csv_rows(stream: io.BytesIO) -> list[dict[str, str]]:
    return list(csv.DictReader(io.TextIOWrapper(stream, encoding="utf-8")))


def archive_rows(archive: zipfile.ZipFile, name: str) -> list[dict[str, str]]:
    with archive.open(PREFIX + name) as stream:
        return csv_rows(stream)


def score(prediction: dict[str, np.ndarray], reference: dict[str, np.ndarray]) -> float:
    np.testing.assert_array_equal(prediction["points_angstrom"], reference["points_angstrom"])
    y = reference["delta_esp_hartree_per_e"]
    p = prediction["predicted_esp_hartree_per_e"]
    return float(np.sqrt(np.mean((p - y) ** 2) / np.mean(y**2)))


def read_xyz(path: Path) -> dict[str, tuple[np.ndarray, int]]:
    """Read the position and fragment fields of this release's simple extxyz."""
    lines = path.read_text().splitlines()
    frames = {}
    i = 0
    while i < len(lines):
        count = int(lines[i])
        header = lines[i + 1]
        fields = dict(token.split("=", 1) for token in header.split() if "=" in token)
        positions = np.asarray([
            [float(value) for value in lines[j].split()[1:4]]
            for j in range(i + 2, i + 2 + count)
        ])
        frames[fields["config_id"]] = positions, int(fields["n_solute_atoms"])
        i += count + 2
    return frames


def bootstrap_paired(
    pair_rows: list[tuple[str, str, float]], n_resamples: int = 100_000
) -> tuple[float, float, float]:
    by_solute = defaultdict(list)
    for _set, solute, difference in pair_rows:
        by_solute[solute].append(difference)
    solutes = sorted(by_solute)
    randomizer = random.Random(20261009)
    estimates = np.empty(n_resamples)
    for index in range(n_resamples):
        selected = [randomizer.choice(solutes) for _ in solutes]
        estimates[index] = statistics.mean(
            difference for solute in selected for difference in by_solute[solute]
        )
    return (
        statistics.mean(row[2] for row in pair_rows),
        float(np.quantile(estimates, 0.025)),
        float(np.quantile(estimates, 0.975)),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    with zipfile.ZipFile(args.archive) as archive:
        if any(name.startswith("/") or ".." in name.split("/") for name in archive.namelist()):
            raise ValueError("Unsafe archive path")
        scores: dict[tuple[str, str, str], float] = {}
        predictions: dict[tuple[str, str, str], dict[str, np.ndarray]] = {}
        metadata: dict[tuple[str, str], str] = {}
        file_count = 0
        maximum_csv_difference = 0.0

        for set_name in SETS:
            archive_base = f"experiments/{set_name}/"
            source_rows = archive_rows(archive, archive_base + "arrow_results.csv")
            point_report = {r["config_id"]: r for r in source_rows}
            smeared_report = {
                r["config_id"]: r
                for r in archive_rows(archive, archive_base + "arrow_smeared_results.csv")
            }
            native_report = {
                (r["config_id"], r["method"]): r
                for r in csv_rows((ROOT / "experiments" / set_name / "results.csv").open("rb"))
            }
            registries = {}
            for method in ("arrow_smeared", "arrow"):
                registry = archive_rows(
                    archive, archive_base + f"predictions/{method}/prediction_registry.csv"
                )
                registries[method] = registry
                for row in registry:
                    config_id = row["config_id"]
                    filename = row["prediction_file"]
                    data = archive.read(PREFIX + archive_base + f"predictions/{method}/" + filename)
                    if hashlib.sha256(data).hexdigest() != row["prediction_sha256"]:
                        raise ValueError(f"Prediction hash mismatch: {set_name}/{config_id}/{method}")
                    with np.load(io.BytesIO(data), allow_pickle=False) as npz:
                        prediction = {key: npz[key] for key in npz.files}
                    with np.load(ROOT / "experiments" / set_name / "references" / filename, allow_pickle=False) as npz:
                        reference = {key: npz[key] for key in npz.files}
                    value = score(prediction, reference)
                    report = (smeared_report if method == "arrow_smeared" else point_report)[config_id]
                    difference = abs(value - float(report["esp_nrmse"]))
                    maximum_csv_difference = max(maximum_csv_difference, difference)
                    if difference > 1e-10:
                        raise ValueError(f"Score mismatch: {set_name}/{config_id}/{method}: {difference}")
                    predictions[(set_name, config_id, method)] = prediction
                    scores[(set_name, config_id, method)] = value
                    metadata[(set_name, config_id)] = row["molecule_id"]
                    file_count += 1
            if {r["config_id"] for r in registries["arrow_smeared"]} != {
                r["config_id"] for r in registries["arrow"]
            }:
                raise ValueError(f"Point and smeared coverage differ in {set_name}")
            for row in registries["arrow_smeared"]:
                config_id = row["config_id"]
                filename = row["prediction_file"]
                with np.load(ROOT / "experiments" / set_name / "references" / filename, allow_pickle=False) as npz:
                    reference = {key: npz[key] for key in npz.files}
                for method in ("glider", "mace_polar_l"):
                    with np.load(ROOT / "experiments" / set_name / "predictions" / method / filename, allow_pickle=False) as npz:
                        prediction = {key: npz[key] for key in npz.files}
                    value = score(prediction, reference)
                    difference = abs(value - float(native_report[(config_id, method)]["esp_nrmse"]))
                    maximum_csv_difference = max(maximum_csv_difference, difference)
                    if difference > 1e-10:
                        raise ValueError(f"Native score mismatch: {set_name}/{config_id}/{method}")
                    scores[(set_name, config_id, method)] = value
                    predictions[(set_name, config_id, method)] = prediction

        grouped = defaultdict(list)
        for (set_name, config_id, method), value in scores.items():
            grouped[(set_name, metadata[(set_name, config_id)], method)].append(value)
        means = {key: statistics.mean(values) for key, values in grouped.items()}
        result = {
            "source_archive_sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest(),
            "prediction_files_verified": file_count,
            "configurations_verified": len(metadata),
            "unique_solutes": len(set(metadata.values())),
            "solute_set_pairs": len(set((set_name, solute) for (set_name, _), solute in metadata.items())),
            "maximum_absolute_score_csv_difference": maximum_csv_difference,
            "sets": {},
            "aggregates": {},
        }
        for set_name in SETS:
            pairs = sorted({(s, mol) for s, mol, method in means if s == set_name})
            result["sets"][set_name] = {
                "configurations": sum(s == set_name for s, _ in metadata),
                "solute_set_pairs": len(pairs),
                "nrmse": {
                    method: statistics.mean(means[(s, mol, method)] for s, mol in pairs)
                    for method in METHODS
                },
            }
        for label, selection in {
            "all_68": SETS,
            "without_shell_44": SETS[:-1],
            "without_panel1_or_shell_40": ("panel_2", "panel_3", "liquid"),
        }.items():
            pairs = sorted({(s, mol) for s, mol, method in means if s in selection})
            paired = [
                (s, mol, means[(s, mol, "arrow_smeared")] - means[(s, mol, "glider")])
                for s, mol in pairs
            ]
            difference, lo, hi = bootstrap_paired(paired)
            result["aggregates"][label] = {
                "solute_set_pairs": len(pairs),
                "unique_solutes": len({mol for _, mol in pairs}),
                "nrmse": {
                    method: statistics.mean(means[(s, mol, method)] for s, mol in pairs)
                    for method in METHODS
                },
                "arrow_minus_glider": difference,
                "solute_cluster_bootstrap_95_percentile_interval": [lo, hi],
            }

        geometry = read_xyz(ROOT / "experiments/shell_size/geometries/configurations.extxyz")
        bias_rows = archive_rows(archive, "shell_size_signed_bias_per_config.csv")
        maximum_bias_difference = 0.0
        checked_bias_rows = 0
        for row in bias_rows:
            key = ("shell_size", row["config_id"], row["method"])
            if key not in predictions:
                continue
            config_id = row["config_id"]
            position, n_solute = geometry[config_id]
            filename = hashlib.sha256(config_id.encode()).hexdigest()[:20] + ".npz"
            with np.load(ROOT / "experiments/shell_size/references" / filename, allow_pickle=False) as npz:
                points = npz["points_angstrom"]
                reference = npz["delta_esp_hartree_per_e"]
            nearest = np.argmin(np.linalg.norm(points[:, None, :] - position[None, :, :], axis=2), axis=1)
            mask = nearest < n_solute
            residual = (predictions[key]["predicted_esp_hartree_per_e"] - reference)[mask]
            value = float(np.mean(residual) * 1000)
            difference = abs(value - float(row["signed_error_solute_surface_mEh"]))
            checked_bias_rows += 1
            maximum_bias_difference = max(maximum_bias_difference, difference)
            if difference > 1e-6:
                raise ValueError(f"Signed-bias mismatch: {config_id}/{row['method']}: {difference}")
        result["maximum_absolute_shell_bias_csv_difference_mEh_per_e"] = maximum_bias_difference
        result["signed_bias_rows_verified"] = checked_bias_rows

    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
