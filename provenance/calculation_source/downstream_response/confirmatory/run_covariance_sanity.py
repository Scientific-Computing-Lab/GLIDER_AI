#!/usr/bin/env python3
"""Check common translation and rotation covariance on the smallest frozen case."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
COMMON = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(1, str(COMMON))
from common import molecule, sha256  # noqa: E402
from evaluate_couplings import METHODS, RESULT_ROOT, evaluate_geometry  # noqa: E402
from run_qm_densities import gpu_scf  # noqa: E402


def arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as values:
        return {key: np.asarray(values[key]) for key in values.files}


def registry_map(directory: Path, identifier: str):
    table = pd.read_csv(directory / "density_registry.csv")
    return {
        str(getattr(row, identifier)): directory / str(row.density_file)
        for row in table.itertuples()
    }


def prediction_map(directory: Path):
    table = pd.read_csv(directory / "prediction_registry.csv")
    return {
        str(row.config_id): directory / str(row.prediction_file)
        for row in table.itertuples()
    }


def transformed_scf_densities(base, probe_positions: np.ndarray):
    symbols = list(base.get_chemical_symbols())
    positions = np.asarray(base.positions)
    n_solute = int(base.info["n_solute_atoms"])
    solute = np.arange(len(base)) < n_solute
    masks = (np.ones(len(base), dtype=bool), solute, ~solute)
    components = [gpu_scf(molecule(symbols, positions, mask)) for mask in masks]
    probe = gpu_scf(molecule(["O", "H", "H"], probe_positions))
    component_density = np.stack([item["density"] for item in components])
    component_dipoles = np.stack([item["dipole_debye"] for item in components])
    return {
        "component_density_matrices": component_density,
        "response_density_matrix": component_density[0] - component_density[1] - component_density[2],
        "response_dipole_debye": component_dipoles[0] - component_dipoles[1] - component_dipoles[2],
        "probe_density_matrix": probe["density"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=RESULT_ROOT / "sanity/covariance")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    bases = read(RESULT_ROOT / "configurations/base_3water.extxyz", index=":")
    # Numerical-only, target-independent choice: the fewest-atom base minimizes the
    # cost of the independent rotated SCFs.
    base = min(bases, key=lambda value: (len(value), str(value.info["case_id"])))
    case_id = str(base.info["case_id"])
    probes = [
        probe
        for probe in read(RESULT_ROOT / "configurations/outer_w4_w12.extxyz", index=":")
        if str(probe.info["case_id"]) == case_id and int(probe.info["water_rank"]) == 4
    ]
    if len(probes) != 1:
        raise RuntimeError("expected exactly one W4 probe for covariance case")
    probe = probes[0]
    probe_id = str(probe.info["probe_id"])

    densities = registry_map(RESULT_ROOT / "qm/base_densities", "case_id")
    probe_densities = registry_map(RESULT_ROOT / "qm/probe_densities", "probe_id")
    glider = prediction_map(RESULT_ROOT / "predictions/glider")
    mace = prediction_map(RESULT_ROOT / "predictions/mace_polar_1_l")
    density = arrays(densities[case_id])
    probe_density = arrays(probe_densities[probe_id])
    glider_values = arrays(glider[case_id])
    mace_values = arrays(mace[case_id])
    probe_positions = probe_density["probe_positions_angstrom"]

    original = evaluate_geometry(
        base=base,
        probe_positions=probe_positions,
        probe_density=probe_density["probe_density_matrix"],
        density_values=density,
        glider_values=glider_values,
        mace_values=mace_values,
        include_components=True,
    )

    translation = np.asarray([1.234, -0.731, 0.419])
    translated_base = base.copy()
    translated_base.positions += translation
    translated = evaluate_geometry(
        base=translated_base,
        probe_positions=probe_positions + translation,
        probe_density=probe_density["probe_density_matrix"],
        density_values=density,
        glider_values=glider_values,
        mace_values=mace_values,
        include_components=True,
    )

    axis = np.asarray([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    matrix = Rotation.from_rotvec(axis * np.radians(37.0)).as_matrix()
    rotated_base = base.copy()
    rotated_base.positions = np.asarray(base.positions) @ matrix.T
    rotated_probe_positions = probe_positions @ matrix.T
    rotated_density = transformed_scf_densities(rotated_base, rotated_probe_positions)
    rotated_glider = dict(glider_values)
    rotated_mace = dict(mace_values)
    rotated_glider["predicted_dipoles_e_bohr"] = glider_values["predicted_dipoles_e_bohr"] @ matrix.T
    rotated_mace["predicted_dipoles_e_bohr"] = mace_values["predicted_dipoles_e_bohr"] @ matrix.T
    rotated = evaluate_geometry(
        base=rotated_base,
        probe_positions=rotated_probe_positions,
        probe_density=rotated_density["probe_density_matrix"],
        density_values=rotated_density,
        glider_values=rotated_glider,
        mace_values=rotated_mace,
        include_components=True,
    )
    expected_dipole = density["response_dipole_debye"] @ matrix.T
    dipole_error = float(np.max(np.abs(rotated_density["response_dipole_debye"] - expected_dipole)))

    array_path = args.output / "rotated_density_matrices.npz"
    np.savez_compressed(
        array_path,
        rotation_matrix=matrix,
        rotated_base_positions_angstrom=rotated_base.positions,
        rotated_probe_positions_angstrom=rotated_probe_positions,
        **rotated_density,
    )
    checks = {}
    for method in METHODS:
        original_value = float(original[f"{method}_energy_kcal_mol"])
        checks[method] = {
            "original_kcal_mol": original_value,
            "translated_kcal_mol": float(translated[f"{method}_energy_kcal_mol"]),
            "translation_abs_error_kcal_mol": abs(float(translated[f"{method}_energy_kcal_mol"]) - original_value),
            "rotated_kcal_mol": float(rotated[f"{method}_energy_kcal_mol"]),
            "rotation_abs_error_kcal_mol": abs(float(rotated[f"{method}_energy_kcal_mol"]) - original_value),
        }
    report = {
        "case_id": case_id,
        "probe_id": probe_id,
        "case_selection": "fewest base atoms, then lexical case ID; numerical-cost criterion only",
        "translation_angstrom": translation.tolist(),
        "rotation_axis": axis.tolist(),
        "rotation_degrees": 37.0,
        "response_dipole_rotation_max_abs_error_debye": dipole_error,
        "methods": checks,
        "rotated_density_file": array_path.name,
        "rotated_density_sha256": sha256(array_path),
        "prospective_claim": False,
    }
    (args.output / "covariance_check.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
