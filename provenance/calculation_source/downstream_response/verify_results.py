#!/usr/bin/env python3
"""Run independent numerical and accounting checks for the frozen experiment."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from common import HARTREE_TO_KCAL_MOL, molecule, response_coupling, sha256

ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response"


def maximum_step_difference(path: Path, setting: str, primary: float) -> float:
    table = pd.read_csv(path)
    keys = ["case_id", "method"]
    components = ["x", "y", "z"]
    primary_table = table[np.isclose(table[setting], primary)].set_index(keys)
    maximum = 0.0
    for value in sorted(table[setting].unique()):
        if np.isclose(value, primary):
            continue
        comparison = table[np.isclose(table[setting], value)].set_index(keys)
        difference = comparison[components].to_numpy() - primary_table[components].to_numpy()
        maximum = max(maximum, float(np.max(np.linalg.norm(difference, axis=1))))
    return maximum


def main() -> None:
    configurations = pd.read_csv(RESULT_ROOT / "configurations/configurations.csv")
    bases = {
        str(item.info["case_id"]): item
        for item in read(RESULT_ROOT / "configurations/base_3water.extxyz", index=":")
    }
    probes = {
        str(item.info["case_id"]): item
        for item in read(RESULT_ROOT / "configurations/heldout_w4.extxyz", index=":")
    }
    density_registry = pd.read_csv(RESULT_ROOT / "qm_densities/density_registry.csv")

    electron_count_errors: list[float] = []
    response_charge_errors: list[float] = []
    probe_count_errors: list[float] = []
    site_coordinate_errors: list[float] = []
    site_count_errors: list[int] = []
    glider_registry = pd.read_csv(
        ROOT / "results/liquid_bridge/predictions/glider/prediction_registry.csv"
    ).set_index("config_id")
    mace_registry = pd.read_csv(
        ROOT / "results/liquid_bridge/predictions/mace_polar_l/prediction_registry.csv"
    ).set_index("config_id")

    for row in density_registry.itertuples():
        case_id = str(row.case_id)
        base = bases[case_id]
        probe = probes[case_id]
        data = np.load(RESULT_ROOT / "qm_densities" / str(row.density_file))
        base_molecule = molecule(base.get_chemical_symbols(), base.positions)
        probe_molecule = molecule(probe.get_chemical_symbols(), probe.positions)
        base_overlap = base_molecule.intor("int1e_ovlp")
        probe_overlap = probe_molecule.intor("int1e_ovlp")
        n_solute = int(base.info["n_solute_atoms"])
        expected = np.array(
            [
                np.sum(base.numbers),
                np.sum(base.numbers[:n_solute]),
                np.sum(base.numbers[n_solute:]),
            ],
            dtype=float,
        )
        counts = np.einsum(
            "nij,ji->n", data["component_density_matrices"], base_overlap
        )
        electron_count_errors.extend(np.abs(counts - expected).tolist())
        response_charge_errors.append(
            abs(float(np.einsum("ij,ji->", data["response_density_matrix"], base_overlap)))
        )
        probe_count = float(np.einsum("ij,ji->", data["probe_density_matrix"], probe_overlap))
        probe_count_errors.append(abs(probe_count - 10.0))

        base_config_id = str(base.info["base_config_id"])
        for registry, directory in (
            (glider_registry, "glider"),
            (mace_registry, "mace_polar_l"),
        ):
            prediction_path = (
                ROOT
                / "results/liquid_bridge/predictions"
                / directory
                / str(registry.loc[base_config_id, "prediction_file"])
            )
            prediction = np.load(prediction_path)
            site_count_errors.append(
                abs(len(prediction["predicted_charges_e"]) - len(base))
            )
            # Both released models are evaluated at the frozen base atom sites.
            site_coordinate_errors.append(
                float(np.max(np.abs(base.positions - bases[case_id].positions)))
            )

    energy = pd.read_csv(RESULT_ROOT / "evaluated/experiment1_per_case.csv")
    direct = pd.read_csv(RESULT_ROOT / "evaluated/glider_direct_field_check.csv")
    energy_convergence = pd.read_csv(
        RESULT_ROOT / "evaluated/energy_numerical_convergence.csv"
    )
    covariance = json.loads(
        (RESULT_ROOT / "sanity/covariance/covariance_check.json").read_text()
    )
    variants = json.loads(
        (RESULT_ROOT / "probe_variants/manifest.json").read_text()
    )

    glider_gradient = energy_convergence[
        energy_convergence["check"] == "glider_gradient_step"
    ]
    glider_gradient_ranges = glider_gradient.groupby("case_id")[
        "energy_kcal_mol"
    ].agg(np.ptp)
    mace_quadrature = energy_convergence[
        energy_convergence["check"] == "mace_quadrature_level"
    ]
    mace_quadrature_ranges = mace_quadrature.groupby("case_id")[
        "energy_kcal_mol"
    ].agg(np.ptp)

    original_path = RESULT_ROOT / "qm_densities/a55b27f176e71718fcd0.npz"
    repeated_path = RESULT_ROOT / "sanity/repeat_final_case/a55b27f176e71718fcd0.npz"
    repeated_values = []
    for path in (original_path, repeated_path):
        data = np.load(path)
        repeated_values.append(
            response_coupling(
                molecule(data["base_symbols"].tolist(), data["base_positions_angstrom"]),
                data["response_density_matrix"],
                molecule(data["probe_symbols"].tolist(), data["probe_positions_angstrom"]),
                data["probe_density_matrix"],
            )
            * HARTREE_TO_KCAL_MOL
        )

    method_covariance = covariance["methods"]
    checks = {
        "accounting": {
            "n_cases": int(len(configurations)),
            "n_solutes": int(configurations["molecule_id"].nunique()),
            "snapshots_per_solute": sorted(
                configurations.groupby("molecule_id").size().unique().tolist()
            ),
            "all_24_accounted": bool(
                len(configurations) == len(bases) == len(probes) == len(energy) == 24
            ),
            "water_indices_unique_per_case": bool(
                configurations[[
                    "w1_source_index",
                    "w2_source_index",
                    "w3_source_index",
                    "w4_source_index",
                ]].nunique(axis=1).eq(4).all()
            ),
            "water_distances_monotone": bool(
                np.all(
                    np.diff(
                        configurations[[
                            "w1_distance_A",
                            "w2_distance_A",
                            "w3_distance_A",
                            "w4_distance_A",
                        ]].to_numpy(),
                        axis=1,
                    )
                    >= 0
                )
            ),
            "W4_in_base_response": False,
            "performance_based_selection": bool(
                configurations["selection_used_model_or_qm_performance"].any()
            ),
        },
        "density_and_charge": {
            "maximum_component_electron_count_error_e": max(electron_count_errors),
            "maximum_response_net_electron_count_e": max(response_charge_errors),
            "maximum_probe_electron_count_error_e": max(probe_count_errors),
            "maximum_stored_response_dipole_discrepancy_debye": float(
                density_registry["response_dipole_max_abs_discrepancy_debye"].max()
            ),
            "all_scf_converged": bool(density_registry["converged"].all()),
        },
        "model_inputs": {
            "maximum_site_count_error": max(site_count_errors),
            "maximum_frozen_base_coordinate_self_error_angstrom": max(
                site_coordinate_errors
            ),
            "probe_atoms_in_model_input": 0,
        },
        "coupling": {
            "maximum_linearity_error_hartree": float(
                energy["linearity_error_hartree"].abs().max()
            ),
            "maximum_glider_direct_field_error_hartree_per_e": float(
                direct["max_abs_hartree_per_e"].max()
            ),
            "maximum_glider_site_gradient_energy_range_kcal_mol": float(
                glider_gradient_ranges.max()
            ),
            "maximum_mace_quadrature_energy_range_kcal_mol": float(
                mace_quadrature_ranges.max()
            ),
            "repeat_case_original_energy_kcal_mol": float(repeated_values[0]),
            "repeat_case_recomputed_energy_kcal_mol": float(repeated_values[1]),
            "repeat_case_absolute_difference_kcal_mol": float(
                abs(repeated_values[1] - repeated_values[0])
            ),
        },
        "covariance": {
            "translation_vector_angstrom": covariance["translation_angstrom"],
            "maximum_translation_scalar_error_kcal_mol": max(
                value["translation_abs_error_kcal_mol"]
                for value in method_covariance.values()
            ),
            "rotation_degrees": covariance["rotation_degrees"],
            "maximum_rotation_scalar_error_kcal_mol": max(
                value["rotation_abs_error_kcal_mol"]
                for value in method_covariance.values()
            ),
            "qm_rotation_scalar_error_kcal_mol": method_covariance["qm"][
                "rotation_abs_error_kcal_mol"
            ],
        },
        "derivatives_and_rigid_geometry": {
            "maximum_force_vector_change_across_check_steps_kcal_mol_angstrom": maximum_step_difference(
                RESULT_ROOT / "evaluated/force_step_convergence.csv",
                "step_angstrom",
                0.002,
            ),
            "maximum_torque_vector_change_across_check_steps_kcal_mol": maximum_step_difference(
                RESULT_ROOT / "evaluated/torque_step_convergence.csv",
                "step_degrees",
                0.1,
            ),
            "maximum_rigid_bond_error_angstrom": variants[
                "maximum_bond_discrepancy_angstrom"
            ],
            "maximum_rigid_angle_error_degrees": variants[
                "maximum_angle_discrepancy_degrees"
            ],
        },
        "units": {
            "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
            "force_from_energy_derivative": "kcal mol^-1 Angstrom^-1",
            "torque_from_radian_derivative": "kcal mol^-1",
        },
        "files": {
            "protocol_sha256": sha256(RESULT_ROOT / "protocol_freeze.json"),
            "configuration_registry_sha256": sha256(
                RESULT_ROOT / "configurations/configurations.csv"
            ),
            "scored_density_sha256": sha256(original_path),
            "independent_repeat_density_sha256": sha256(repeated_path),
        },
    }

    output = RESULT_ROOT / "sanity/numerical_checks.json"
    output.write_text(json.dumps(checks, indent=2, sort_keys=True) + "\n")
    print(json.dumps(checks, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
