#!/usr/bin/env python3
"""Verify confirmatory accounting, physics identities, and numerical stability."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read

ROOT = Path(__file__).resolve().parents[4]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response_confirmatory"


def max_step_vector_change(table: pd.DataFrame, step_column: str, primary: float) -> float:
    maximum = 0.0
    keys = ["probe_id", "method"]
    for _, group in table.groupby(keys, sort=False):
        primary_row = group[np.isclose(group[step_column], primary)]
        if len(primary_row) != 1:
            raise RuntimeError(f"missing primary derivative row: {group[keys].iloc[0].to_dict()}")
        primary_vector = primary_row[["x", "y", "z"]].to_numpy()[0]
        for row in group.itertuples():
            vector = np.asarray([row.x, row.y, row.z])
            maximum = max(maximum, float(np.linalg.norm(vector - primary_vector)))
    return maximum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=RESULT_ROOT / "sanity/numerical_checks.json")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    bases = read(RESULT_ROOT / "configurations/base_3water.extxyz", index=":")
    probes = read(RESULT_ROOT / "configurations/outer_w4_w12.extxyz", index=":")
    configuration = pd.read_csv(RESULT_ROOT / "configurations/configurations.csv")
    outer = pd.read_csv(RESULT_ROOT / "configurations/outer_probes.csv")
    base_density = pd.read_csv(RESULT_ROOT / "qm/base_densities/density_registry.csv")
    probe_density = pd.read_csv(RESULT_ROOT / "qm/probe_densities/density_registry.csv")
    variants = pd.read_csv(RESULT_ROOT / "qm/probe_variants/variant_registry.csv")
    energy = pd.read_csv(RESULT_ROOT / "evaluated/range_energy_per_probe.csv")
    force = pd.read_csv(RESULT_ROOT / "evaluated/force_step_convergence.csv")
    torque = pd.read_csv(RESULT_ROOT / "evaluated/torque_step_convergence.csv")
    direct = pd.read_csv(RESULT_ROOT / "evaluated/glider_direct_field_check.csv")
    energy_convergence = pd.read_csv(RESULT_ROOT / "evaluated/energy_numerical_convergence.csv")
    covariance = json.loads((RESULT_ROOT / "sanity/covariance/covariance_check.json").read_text())

    if not (
        len(bases) == len(configuration) == len(base_density) == 10
        and len(probes) == len(outer) == len(probe_density) == len(energy) == 90
    ):
        raise RuntimeError("confirmatory case accounting failed")
    if sorted(outer.groupby("case_id").water_rank.apply(list).tolist()) != [list(range(4, 13))] * 10:
        raise RuntimeError("W4-W12 rank accounting failed")
    if not all(int(base.info["n_waters"]) == 3 for base in bases):
        raise RuntimeError("an outer probe entered a base model/QM input")

    prediction_site_errors = []
    for directory in (RESULT_ROOT / "predictions/glider", RESULT_ROOT / "predictions/mace_polar_1_l"):
        registry = pd.read_csv(directory / "prediction_registry.csv")
        for row in registry.itertuples():
            base = next(value for value in bases if str(value.info["case_id"]) == row.config_id)
            with np.load(directory / row.prediction_file, allow_pickle=False) as values:
                prediction_site_errors.append(abs(len(values["predicted_charges_e"]) - len(base)))

    convergence_ranges = (
        energy_convergence.groupby(["probe_id", "check"]).energy_kcal_mol.agg(lambda x: float(np.ptp(x)))
    )
    glider_ranges = convergence_ranges[convergence_ranges.index.get_level_values("check") == "glider_gradient_step_bohr"]
    mace_ranges = convergence_ranges[convergence_ranges.index.get_level_values("check") == "mace_quadrature_level"]

    translation_error = max(
        float(values["translation_abs_error_kcal_mol"])
        for values in covariance["methods"].values()
    )
    rotation_error = max(
        float(values["rotation_abs_error_kcal_mol"])
        for values in covariance["methods"].values()
    )
    checks = {
        "accounting": {
            "all_10_bases_accounted": True,
            "all_90_outer_probes_accounted": True,
            "n_solutes": int(configuration.molecule_id.nunique()),
            "n_source_families": int(configuration.family.nunique()),
            "W4_W12_in_base_response_or_model_input": False,
            "water_ranks_complete": True,
            "water_indices_unique_within_case": bool(
                outer.groupby("case_id").source_water_index.nunique().eq(9).all()
            ),
            "water_distances_monotone": bool(
                outer.sort_values(["case_id", "water_rank"])
                .groupby("case_id")
                .oxygen_distance_A.apply(lambda x: np.all(np.diff(x) >= 0))
                .all()
            ),
        },
        "density_and_charge": {
            "all_base_scf_converged": bool(base_density.converged.all()),
            "all_probe_scf_converged": bool(probe_density.converged.all()),
            "all_variant_scf_converged": bool(variants.converged.all()),
            "maximum_abs_response_electron_number": float(base_density.response_electron_number.abs().max()),
            "maximum_response_dipole_reconstruction_error_debye": float(
                base_density.response_density_dipole_max_abs_discrepancy_debye.max()
            ),
        },
        "coupling": {
            "maximum_density_linearity_error_hartree": float(energy.linearity_error_hartree.abs().max()),
            "maximum_glider_direct_field_error_hartree_per_e": float(direct.max_abs_hartree_per_e.max()),
            "maximum_glider_site_gradient_energy_range_kcal_mol": float(glider_ranges.max()),
            "maximum_mace_quadrature_energy_range_kcal_mol": float(mace_ranges.max()),
        },
        "derivatives_and_rigid_geometry": {
            "maximum_force_vector_change_across_check_steps_kcal_mol_angstrom": max_step_vector_change(force, "step_angstrom", 0.002),
            "maximum_torque_vector_change_across_check_steps_kcal_mol": max_step_vector_change(torque, "step_degrees", 0.1),
            "maximum_rigid_bond_error_angstrom": float(variants.bond_length_max_discrepancy_angstrom.max()),
            "maximum_rigid_angle_error_degrees": float(variants.angle_discrepancy_degrees.max()),
        },
        "covariance": {
            "case_id": covariance["case_id"],
            "maximum_translation_scalar_error_kcal_mol": translation_error,
            "maximum_rotation_scalar_error_kcal_mol": rotation_error,
            "response_dipole_rotation_max_abs_error_debye": covariance[
                "response_dipole_rotation_max_abs_error_debye"
            ],
        },
        "model_inputs": {
            "maximum_prediction_site_count_error": int(max(prediction_site_errors)),
            "outer_probe_atoms_in_model_input": 0,
        },
        "units": {
            "hartree_to_kcal_mol": 627.5094740631,
            "force": "kcal mol^-1 Angstrom^-1",
            "torque": "kcal mol^-1 (derivative with respect to radians)",
        },
    }
    if checks["density_and_charge"]["maximum_abs_response_electron_number"] > 1.0e-6:
        raise RuntimeError("response charge tolerance exceeded")
    if checks["density_and_charge"]["maximum_response_dipole_reconstruction_error_debye"] > 1.0e-5:
        raise RuntimeError("response dipole reconstruction tolerance exceeded")
    if checks["model_inputs"]["maximum_prediction_site_count_error"] != 0:
        raise RuntimeError("prediction site accounting failed")
    args.output.write_text(json.dumps(checks, indent=2, sort_keys=True) + "\n")
    print(json.dumps(checks, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
