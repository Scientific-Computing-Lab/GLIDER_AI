#!/usr/bin/env python3
"""Validate the frozen exact-QM global-multipole hierarchy before scoring."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from multipole_control import (
    BOHR_TO_ANGSTROM,
    CONTROL_ROOT,
    DEBYE_PER_E_BOHR,
    HARTREE_TO_KCAL_MOL,
    ORIGIN_TO_EXISTING_DIPOLE,
    ORIGINS,
    MultipoleData,
    case_origins,
    hierarchy_couplings,
    irreducible_moments,
    molecule,
    point_multipole_coupling,
    raw_charge_moments,
    translate_raw_moments,
)


def tensor_norm(values: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(values).ravel()))


def moment_rows(data: MultipoleData) -> tuple[pd.DataFrame, dict]:
    rows = []
    maximum_charge = 0.0
    maximum_dipole_difference = 0.0
    maximum_trace_q = 0.0
    maximum_trace_o = 0.0
    maximum_translation = {"dipole": 0.0, "second": 0.0, "third": 0.0}
    for case_id, base in sorted(data.bases.items()):
        density_values = data.base_values(case_id)
        base_molecule = molecule(list(base.get_chemical_symbols()), np.asarray(base.positions))
        origins = case_origins(base)
        raw_by_origin = {
            name: raw_charge_moments(
                base_molecule, density_values["response_density_matrix"], origin
            )
            for name, origin in origins.items()
        }
        baseline = raw_by_origin["base_com"]
        for origin_name, raw in raw_by_origin.items():
            irreducible = irreducible_moments(raw)
            quadrupole = irreducible["quadrupole"]
            octupole = irreducible["octupole"]
            dipole_debye = np.asarray(raw["dipole"]) * DEBYE_PER_E_BOHR
            discrepancy = float(
                np.max(
                    np.abs(dipole_debye - density_values["response_dipole_debye"])
                )
            )
            maximum_charge = max(maximum_charge, abs(float(raw["charge"])))
            maximum_dipole_difference = max(maximum_dipole_difference, discrepancy)
            maximum_trace_q = max(
                maximum_trace_q, abs(float(np.trace(quadrupole)))
            )
            octupole_trace = np.einsum("iaa->i", octupole)
            maximum_trace_o = max(
                maximum_trace_o, float(np.max(np.abs(octupole_trace)))
            )
            if origin_name != "base_com":
                displacement = (
                    origins[origin_name] - origins["base_com"]
                ) / BOHR_TO_ANGSTROM
                translated = translate_raw_moments(baseline, displacement)
                for key in maximum_translation:
                    maximum_translation[key] = max(
                        maximum_translation[key],
                        float(
                            np.max(
                                np.abs(
                                    np.asarray(raw[key])
                                    - np.asarray(translated[key])
                                )
                            )
                        ),
                    )
            rows.append(
                {
                    "case_id": case_id,
                    "molecule_id": str(base.info["molecule_id"]),
                    "family": str(base.info["family"]),
                    "origin": origin_name,
                    "origin_x_angstrom": origins[origin_name][0],
                    "origin_y_angstrom": origins[origin_name][1],
                    "origin_z_angstrom": origins[origin_name][2],
                    "physical_response_charge_e": float(raw["charge"]),
                    "dipole_x_e_bohr": float(np.asarray(raw["dipole"])[0]),
                    "dipole_y_e_bohr": float(np.asarray(raw["dipole"])[1]),
                    "dipole_z_e_bohr": float(np.asarray(raw["dipole"])[2]),
                    "dipole_max_discrepancy_debye": discrepancy,
                    "quadrupole_frobenius_e_bohr2": tensor_norm(quadrupole),
                    "quadrupole_trace_e_bohr2": float(np.trace(quadrupole)),
                    "octupole_frobenius_e_bohr3": tensor_norm(octupole),
                    "octupole_max_trace_e_bohr3": float(
                        np.max(np.abs(octupole_trace))
                    ),
                    **{
                        f"Q_{i}{j}_e_bohr2": float(quadrupole[i, j])
                        for i in range(3)
                        for j in range(3)
                    },
                    **{
                        f"O_{i}{j}{k}_e_bohr3": float(octupole[i, j, k])
                        for i in range(3)
                        for j in range(3)
                        for k in range(3)
                    },
                }
            )
    checks = {
        "maximum_abs_physical_response_charge_e": maximum_charge,
        "maximum_density_vs_stored_dipole_discrepancy_debye": maximum_dipole_difference,
        "maximum_abs_traceless_quadrupole_trace_e_bohr2": maximum_trace_q,
        "maximum_abs_traceless_octupole_trace_e_bohr3": maximum_trace_o,
        "maximum_origin_translation_abs_errors": maximum_translation,
    }
    return pd.DataFrame(rows), checks


def exact_dipole_energy(
    base,
    density_values: dict[str, np.ndarray],
    positions: np.ndarray,
    density: np.ndarray,
) -> dict[str, float]:
    probe_molecule = molecule(["O", "H", "H"], positions)
    dipole = density_values["response_dipole_debye"] / DEBYE_PER_E_BOHR
    return {
        origin_name: point_multipole_coupling(
            origin[None],
            np.zeros(1),
            dipole[None],
            probe_molecule,
            density,
            gradient_step_bohr=1.0e-4,
        )
        * HARTREE_TO_KCAL_MOL
        for origin_name, origin in case_origins(base).items()
    }


def order1_reproduction(data: MultipoleData) -> dict:
    existing_energy = pd.read_csv(
        data.base_density_files[next(iter(data.base_density_files))].parents[2]
        / "evaluated/range_energy_per_probe.csv"
    ).set_index("probe_id")
    existing_torque = pd.read_csv(
        data.base_density_files[next(iter(data.base_density_files))].parents[2]
        / "evaluated/breadth_torque_per_solute.csv"
    )
    existing_orientation = pd.read_csv(
        data.base_density_files[next(iter(data.base_density_files))].parents[2]
        / "evaluated/orientation_profiles.csv"
    )
    maximum_energy = 0.0
    for probe_id in sorted(data.probes):
        probe = data.probes[probe_id]
        case_id = str(probe.info["case_id"])
        values = data.probe_values(probe_id)
        calculated = exact_dipole_energy(
            data.bases[case_id],
            data.base_values(case_id),
            values["probe_positions_angstrom"],
            values["probe_density_matrix"],
        )
        for origin_name, value in calculated.items():
            column = f"{ORIGIN_TO_EXISTING_DIPOLE[origin_name]}_energy_kcal_mol"
            maximum_energy = max(
                maximum_energy,
                abs(value - float(existing_energy.loc[probe_id, column])),
            )

    maximum_torque = 0.0
    for probe_id in data.w4_ids():
        case_id = str(data.probes[probe_id].info["case_id"])
        energies = {origin: np.empty((3, 2)) for origin in ORIGINS}
        for axis in range(3):
            for sign_index, label in enumerate(("minus", "plus")):
                variant = data.variant_values(
                    probe_id, f"torque_0.1deg_axis{axis}_{label}"
                )
                calculated = exact_dipole_energy(
                    data.bases[case_id],
                    data.base_values(case_id),
                    variant["positions_angstrom"],
                    variant["density_matrix"],
                )
                for origin_name, value in calculated.items():
                    energies[origin_name][axis, sign_index] = value
        for origin_name, values in energies.items():
            torque = -(values[:, 1] - values[:, 0]) / (2.0 * np.radians(0.1))
            method = ORIGIN_TO_EXISTING_DIPOLE[origin_name]
            existing = existing_torque[
                (existing_torque.probe_id == probe_id)
                & (existing_torque.method == method)
            ][["x", "y", "z"]].to_numpy()[0]
            maximum_torque = max(
                maximum_torque, float(np.max(np.abs(torque - existing)))
            )

    maximum_orientation = 0.0
    for probe_id in data.w4_ids():
        case_id = str(data.probes[probe_id].info["case_id"])
        for rotation_index in range(24):
            rotation_id = f"O_{rotation_index:02d}"
            variant = data.variant_values(
                probe_id, f"orientation_{rotation_id}"
            )
            calculated = exact_dipole_energy(
                data.bases[case_id],
                data.base_values(case_id),
                variant["positions_angstrom"],
                variant["density_matrix"],
            )
            for origin_name, value in calculated.items():
                method = ORIGIN_TO_EXISTING_DIPOLE[origin_name]
                existing = existing_orientation[
                    (existing_orientation.probe_id == probe_id)
                    & (existing_orientation.rotation_id == rotation_id)
                    & (existing_orientation.method == method)
                ].energy_kcal_mol.iloc[0]
                maximum_orientation = max(
                    maximum_orientation, abs(value - float(existing))
                )
    return {
        "maximum_initial_W4_W12_energy_difference_kcal_mol": maximum_energy,
        "maximum_W4_torque_component_difference_kcal_mol": maximum_torque,
        "maximum_W4_orientation_energy_difference_kcal_mol": maximum_orientation,
    }


def rotation_covariance(data: MultipoleData) -> dict:
    case_id = "confirmatory__methane__base3"
    base = data.bases[case_id]
    original_values = data.base_values(case_id)
    rotated_path = (
        data.base_density_files[case_id].parents[2]
        / "sanity/covariance/rotated_density_matrices.npz"
    )
    with np.load(rotated_path, allow_pickle=False) as stored:
        rotated_values = {key: np.asarray(stored[key]) for key in stored.files}
    rotation = rotated_values["rotation_matrix"]
    rotated_base = Atoms(
        symbols=base.get_chemical_symbols(),
        positions=rotated_values["rotated_base_positions_angstrom"],
    )
    rotated_base.info.update(base.info)
    original_molecule = molecule(
        list(base.get_chemical_symbols()), np.asarray(base.positions)
    )
    rotated_molecule = molecule(
        list(rotated_base.get_chemical_symbols()), np.asarray(rotated_base.positions)
    )
    rows = {}
    for origin_name in ORIGINS:
        original_origin = case_origins(base)[origin_name]
        rotated_origin = case_origins(rotated_base)[origin_name]
        original_raw = raw_charge_moments(
            original_molecule,
            original_values["response_density_matrix"],
            original_origin,
        )
        rotated_raw = raw_charge_moments(
            rotated_molecule,
            rotated_values["response_density_matrix"],
            rotated_origin,
        )
        original_irreducible = irreducible_moments(original_raw)
        rotated_irreducible = irreducible_moments(rotated_raw)
        expected_p = rotation @ np.asarray(original_raw["dipole"])
        expected_q = rotation @ original_irreducible["quadrupole"] @ rotation.T
        expected_o = np.einsum(
            "ia,jb,kc,abc->ijk",
            rotation,
            rotation,
            rotation,
            original_irreducible["octupole"],
        )
        values = {}
        for label, actual, expected in (
            ("dipole", np.asarray(rotated_raw["dipole"]), expected_p),
            ("quadrupole", rotated_irreducible["quadrupole"], expected_q),
            ("octupole", rotated_irreducible["octupole"], expected_o),
        ):
            error = tensor_norm(actual - expected)
            scale = tensor_norm(expected)
            values[f"{label}_frobenius_error"] = error
            values[f"{label}_relative_error"] = error / scale if scale else float("nan")
        rows[origin_name] = values
    return rows


def derivative_convergence(data: MultipoleData) -> list[dict]:
    rows = []
    for probe_id in (
        "confirmatory__methane__base3__W4",
        "confirmatory__2OJ9-original__base3__W4",
    ):
        probe = data.probes[probe_id]
        case_id = str(probe.info["case_id"])
        probe_values = data.probe_values(probe_id)
        for step in (0.0075, 0.01, 0.015):
            couplings = hierarchy_couplings(
                base=data.bases[case_id],
                density_values=data.base_values(case_id),
                probe_positions_angstrom=probe_values["probe_positions_angstrom"],
                probe_density=probe_values["probe_density_matrix"],
                step_bohr=step,
            )
            for (hierarchy, origin_name), values in couplings.items():
                rows.append(
                    {
                        "probe_id": probe_id,
                        "step_bohr": step,
                        "hierarchy": hierarchy,
                        "origin": origin_name,
                        "energy_kcal_mol": values["energy_kcal_mol"],
                    }
                )
    table = pd.DataFrame(rows)
    changes = []
    for (probe_id, hierarchy, origin), group in table.groupby(
        ["probe_id", "hierarchy", "origin"], sort=False
    ):
        primary = float(group[np.isclose(group.step_bohr, 0.01)].energy_kcal_mol.iloc[0])
        changes.append(
            {
                "probe_id": probe_id,
                "hierarchy": hierarchy,
                "origin": origin,
                "maximum_abs_change_from_primary_kcal_mol": float(
                    np.max(np.abs(group.energy_kcal_mol.to_numpy() - primary))
                ),
            }
        )
    return changes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=CONTROL_ROOT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = MultipoleData()
    moments, moment_checks = moment_rows(data)
    moments.to_csv(args.output / "multipole_moments.csv", index=False)
    validation = {
        "accounting": {
            "n_bases": len(data.bases),
            "n_initial_probes": len(data.probes),
            "n_W4": len(data.w4_ids()),
            "all_required_cases_present": len(data.bases) == 10
            and len(data.probes) == 90
            and len(data.w4_ids()) == 10,
        },
        "physical_charge_and_moments": moment_checks,
        "order1_reproduction": order1_reproduction(data),
        "rotation_covariance": rotation_covariance(data),
        "derivative_convergence": derivative_convergence(data),
        "comparative_performance_inspected": False,
    }
    (args.output / "validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
