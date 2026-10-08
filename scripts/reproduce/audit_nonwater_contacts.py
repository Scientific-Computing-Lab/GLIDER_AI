#!/usr/bin/env python3
"""Audit the frozen non-water panel's geometry and interaction energies.

The original 72-case panel is retained intact.  This is a post hoc diagnostic,
not a new prospective test or a filter for the published primary estimate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from ase import __version__ as ase_version
from ase.data import vdw_radii
from ase.io import read


ROOT = Path(__file__).resolve().parents[2]
HARTREE_TO_KCAL_MOL = 627.509474
COMPONENT_ORDER = [
    "full_complex",
    "solute_with_perturbant_ghosts",
    "perturbant_with_solute_ghosts",
]
METHODS = ("glider", "mace_polar_l", "mace_polar_m")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mean(values: list[float]) -> float:
    return float(np.mean(values))


def subset_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    per_solute: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        sid = str(row["molecule_id"])
        for method in METHODS:
            per_solute[sid][method].append(float(row[f"{method}_esp_nrmse"]))
    equal_solute = {
        method: mean([mean(values[method]) for values in per_solute.values()])
        for method in METHODS
    }
    return {
        "configurations": len(rows),
        "solutes": len(per_solute),
        "case_wins_glider_vs_mace_l": sum(
            float(row["glider_esp_nrmse"]) < float(row["mace_polar_l_esp_nrmse"])
            for row in rows
        ),
        "equal_solute_esp_nrmse": equal_solute,
        "relative_reduction_vs_mace_l": (
            1.0 - equal_solute["glider"] / equal_solute["mace_polar_l"]
        ),
    }


def audit(base: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    geometry_path = base / "geometries/configuration_registry.csv"
    reference_registry_path = base / "references/REFERENCE_REGISTRY.csv"
    scores_path = base / "results.csv"
    geometry = {row["config_id"]: row for row in read_csv(geometry_path)}
    frames = {str(atoms.info["config_id"]): atoms for atoms in read(base / "geometries/configurations.extxyz", index=":")}
    references = {row["config_id"]: row for row in read_csv(reference_registry_path)}
    scores: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in read_csv(scores_path):
        if row["method"] in METHODS:
            scores[row["config_id"]][row["method"]] = row
    ids = set(geometry)
    if not (len(ids) == 72 and ids == set(references) == set(scores) == set(frames)):
        raise ValueError("Geometry, QM reference and score registries do not match 72 cases")

    rows: list[dict[str, object]] = []
    for config_id in sorted(ids):
        geo = geometry[config_id]
        atoms = frames[config_id]
        n_solute = int(geo["n_solute_atoms"])
        positions = np.asarray(atoms.positions, dtype=float)
        numbers = np.asarray(atoms.numbers, dtype=int)
        pair_distances = np.linalg.norm(
            positions[:n_solute, None, :] - positions[None, n_solute:, :], axis=2
        )
        pair_vdw_radii = (
            vdw_radii[numbers[:n_solute, None]] + vdw_radii[numbers[None, n_solute:]]
        )
        if not np.all(np.isfinite(pair_vdw_radii)):
            raise ValueError(f"Undefined ASE van der Waals radius for {config_id}")
        vdw_ratio = pair_distances / pair_vdw_radii
        if abs(float(np.min(pair_distances)) - float(geo["minimum_interfragment_distance_A"])) > 1e-5:
            raise ValueError(f"Geometry registry distance mismatch for {config_id}")
        ref = references[config_id]
        if set(scores[config_id]) != set(METHODS):
            raise ValueError(f"Missing primary scores for {config_id}")
        array_path = base / "references" / ref["observable_file"]
        if sha256(array_path) != ref["observable_sha256"]:
            raise ValueError(f"Reference hash mismatch for {config_id}")
        record_path = array_path.with_suffix(".json")
        record = json.loads(record_path.read_text())
        if record["config_id"] != config_id or record["component_order"] != COMPONENT_ORDER:
            raise ValueError(f"Energy component order mismatch for {config_id}")
        with np.load(array_path, allow_pickle=False) as arrays:
            energies = np.asarray(arrays["component_energy_hartree"], dtype=float)
        if energies.shape != (3,) or not np.all(np.isfinite(energies)):
            raise ValueError(f"Bad component energies for {config_id}")
        interaction = float((energies[0] - energies[1] - energies[2]) * HARTREE_TO_KCAL_MOL)
        row: dict[str, object] = {
            "config_id": config_id,
            "molecule_id": geo["molecule_id"],
            "perturbant": geo["perturbant"],
            "regime": geo["regime"],
            "cp_interaction_energy_kcal_mol": interaction,
            "minimum_interfragment_distance_A": float(geo["minimum_interfragment_distance_A"]),
            "minimum_interfragment_covalent_ratio": float(
                geo["minimum_interfragment_covalent_ratio"]
            ),
            "minimum_interfragment_vdw_ratio_ase": float(np.min(vdw_ratio)),
            "clearance_threshold_met": geo["clearance_threshold_met"],
        }
        for method in METHODS:
            row[f"{method}_esp_nrmse"] = float(scores[config_id][method]["esp_nrmse"])
        row["glider_minus_mace_l_nrmse"] = (
            float(row["glider_esp_nrmse"]) - float(row["mace_polar_l_esp_nrmse"])
        )
        rows.append(row)

    energies = np.asarray([float(row["cp_interaction_energy_kcal_mol"]) for row in rows])
    distance = np.asarray([float(row["minimum_interfragment_distance_A"]) for row in rows])
    ratio = np.asarray([float(row["minimum_interfragment_covalent_ratio"]) for row in rows])
    vdw_ratio = np.asarray([float(row["minimum_interfragment_vdw_ratio_ase"]) for row in rows])
    summary = {
        "status": "post hoc geometry diagnostic; the frozen 72-case panel is unchanged",
        "energy_definition": "(E_complex - E_solute+ghosts - E_perturbant+ghosts) * 627.509474",
        "vdw_radius_source": f"ASE {ase_version} ase.data.vdw_radii; ratio is shortest interfragment distance divided by the corresponding sum of radii",
        "source_sha256": {
            "geometry_registry": sha256(geometry_path),
            "reference_registry": sha256(reference_registry_path),
            "scores": sha256(scores_path),
        },
        "all": subset_summary(rows),
        "cp_energy_kcal_mol": {
            "minimum": float(np.min(energies)),
            "median": float(np.median(energies)),
            "maximum": float(np.max(energies)),
            "positive_count": int(np.count_nonzero(energies > 0)),
            "above_10_count": int(np.count_nonzero(energies > 10)),
        },
        "minimum_interfragment_distance_A": {
            "minimum": float(np.min(distance)),
            "median": float(np.median(distance)),
        },
        "minimum_interfragment_covalent_ratio": {
            "minimum": float(np.min(ratio)),
            "median": float(np.median(ratio)),
        },
        "minimum_interfragment_vdw_ratio_ase": {
            "minimum": float(np.min(vdw_ratio)),
            "median": float(np.median(vdw_ratio)),
            "below_1_count": int(np.count_nonzero(vdw_ratio < 1)),
            "below_0_8_count": int(np.count_nonzero(vdw_ratio < 0.8)),
        },
        "descriptive_post_hoc_subsets": {
            "attractive_cp_energy_below_0": subset_summary(
                [row for row in rows if float(row["cp_interaction_energy_kcal_mol"]) < 0]
            ),
            "cp_energy_below_10": subset_summary(
                [row for row in rows if float(row["cp_interaction_energy_kcal_mol"]) < 10]
            ),
            "cp_energy_above_10": subset_summary(
                [row for row in rows if float(row["cp_interaction_energy_kcal_mol"]) > 10]
            ),
        },
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "experiments/nonwater")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.base
    rows, summary = audit(args.base)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "contact_audit.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "contact_audit_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
