#!/usr/bin/env python3
"""Check scalar-energy translation and rotation covariance on the frozen first case."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from common import molecule, sha256
from evaluate_couplings import METHODS, evaluate_geometry
from run_qm_densities import gpu_scf
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response"


def load_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as values:
        return {key: np.asarray(values[key]) for key in values.files}


def registry_map(directory: Path, registry: str, identifier: str, filename: str):
    table = pd.read_csv(directory / registry)
    return {
        str(getattr(row, identifier)): directory / str(getattr(row, filename))
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
        "response_density_matrix": (
            component_density[0] - component_density[1] - component_density[2]
        ),
        "response_dipole_debye": (
            component_dipoles[0] - component_dipoles[1] - component_dipoles[2]
        ),
        "probe_density_matrix": probe["density"],
        "component_energy_hartree": np.asarray(
            [item["energy_hartree"] for item in components]
        ),
        "probe_energy_hartree": probe["energy_hartree"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=RESULT_ROOT / "sanity/covariance"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    bases = sorted(
        read(RESULT_ROOT / "configurations/base_3water.extxyz", index=":"),
        key=lambda item: str(item.info["case_id"]),
    )
    probes = {
        str(item.info["case_id"]): item
        for item in read(RESULT_ROOT / "configurations/heldout_w4.extxyz", index=":")
    }
    base = bases[0]
    case_id = str(base.info["case_id"])
    base_config_id = str(base.info["base_config_id"])
    density_files = registry_map(
        RESULT_ROOT / "qm_densities",
        "density_registry.csv",
        "case_id",
        "density_file",
    )
    glider_files = registry_map(
        ROOT / "results/liquid_bridge/predictions/glider",
        "prediction_registry.csv",
        "config_id",
        "prediction_file",
    )
    mace_files = registry_map(
        ROOT / "results/liquid_bridge/predictions/mace_polar_l",
        "prediction_registry.csv",
        "config_id",
        "prediction_file",
    )
    density = load_arrays(density_files[case_id])
    glider = load_arrays(glider_files[base_config_id])
    mace = load_arrays(mace_files[base_config_id])
    probe_positions = np.asarray(probes[case_id].positions)

    original = evaluate_geometry(
        base=base,
        probe_positions=probe_positions,
        probe_density=density["probe_density_matrix"],
        density_values=density,
        glider_values=glider,
        mace_values=mace,
        include_components=True,
    )

    translation = np.asarray([1.234, -0.731, 0.419])
    translated_base = base.copy()
    translated_base.positions += translation
    translated = evaluate_geometry(
        base=translated_base,
        probe_positions=probe_positions + translation,
        probe_density=density["probe_density_matrix"],
        density_values=density,
        glider_values=glider,
        mace_values=mace,
        include_components=True,
    )

    axis = np.asarray([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    matrix = Rotation.from_rotvec(axis * np.radians(37.0)).as_matrix()
    rotated_base = base.copy()
    rotated_base.positions = np.asarray(base.positions) @ matrix.T
    rotated_probe_positions = probe_positions @ matrix.T
    rotated_density = transformed_scf_densities(rotated_base, rotated_probe_positions)
    rotated_glider = dict(glider)
    rotated_mace = dict(mace)
    rotated_glider["predicted_dipoles_e_bohr"] = (
        glider["predicted_dipoles_e_bohr"] @ matrix.T
    )
    rotated_mace["predicted_dipoles_e_bohr"] = (
        mace["predicted_dipoles_e_bohr"] @ matrix.T
    )
    rotated = evaluate_geometry(
        base=rotated_base,
        probe_positions=rotated_probe_positions,
        probe_density=rotated_density["probe_density_matrix"],
        density_values=rotated_density,
        glider_values=rotated_glider,
        mace_values=rotated_mace,
        include_components=True,
    )
    response_dipole_expected = density["response_dipole_debye"] @ matrix.T
    response_dipole_rotation_error = float(
        np.max(
            np.abs(
                rotated_density["response_dipole_debye"]
                - response_dipole_expected
            )
        )
    )

    array_path = args.output / "rotated_density_matrices.npz"
    np.savez_compressed(
        array_path,
        rotation_matrix=matrix,
        rotated_base_positions_angstrom=rotated_base.positions,
        rotated_probe_positions_angstrom=rotated_probe_positions,
        **rotated_density,
    )
    method_checks = {}
    for method in METHODS:
        original_value = float(original[f"{method}_energy_kcal_mol"])
        method_checks[method] = {
            "original_kcal_mol": original_value,
            "translated_kcal_mol": float(
                translated[f"{method}_energy_kcal_mol"]
            ),
            "translation_abs_error_kcal_mol": abs(
                float(translated[f"{method}_energy_kcal_mol"]) - original_value
            ),
            "rotated_kcal_mol": float(rotated[f"{method}_energy_kcal_mol"]),
            "rotation_abs_error_kcal_mol": abs(
                float(rotated[f"{method}_energy_kcal_mol"]) - original_value
            ),
        }
    report = {
        "case_id": case_id,
        "translation_angstrom": translation.tolist(),
        "rotation_axis": axis.tolist(),
        "rotation_degrees": 37.0,
        "response_dipole_rotation_max_abs_error_debye": response_dipole_rotation_error,
        "methods": method_checks,
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
