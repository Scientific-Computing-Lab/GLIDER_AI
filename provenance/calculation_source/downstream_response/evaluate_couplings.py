#!/usr/bin/env python3
"""Evaluate frozen downstream energies, forces, torques, and orientations."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from common import (
    DEBYE_PER_E_BOHR,
    HARTREE_TO_KCAL_MOL,
    base_density_couplings,
    electron_density_at_nuclei_coupling,
    gaussian_multipole_coupling,
    molecule,
    point_multipole_coupling,
    sha256,
)
from scipy.spatial.transform import Rotation

from glider.inference import esp_from_sites

ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response"
METHODS = (
    "qm",
    "glider",
    "mace_polar_l",
    "exact_dipole_base_com",
    "exact_dipole_solute_com",
    "exact_dipole_nuclear_center",
    "zero",
)
FORCE_STEPS = (0.001, 0.002, 0.004)
TORQUE_STEPS = (0.05, 0.1, 0.2)


def load_file_map(directory: Path, registry_name: str, id_column: str, file_column: str):
    table = pd.read_csv(directory / registry_name)
    return {
        str(getattr(row, id_column)): directory / str(getattr(row, file_column))
        for row in table.itertuples()
    }


def load_predictions(directory: Path) -> dict[str, Path]:
    return load_file_map(
        directory, "prediction_registry.csv", "config_id", "prediction_file"
    )


def centre_of_mass(positions: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.average(np.asarray(positions), axis=0, weights=np.asarray(weights))


def case_origins(base, n_solute: int) -> dict[str, np.ndarray]:
    positions = np.asarray(base.positions)
    masses = np.asarray(base.get_masses())
    nuclear_charges = np.asarray(base.get_atomic_numbers(), dtype=float)
    return {
        "exact_dipole_base_com": centre_of_mass(positions, masses),
        "exact_dipole_solute_com": centre_of_mass(
            positions[:n_solute], masses[:n_solute]
        ),
        "exact_dipole_nuclear_center": centre_of_mass(positions, nuclear_charges),
    }


def nuclear_nuclear_coupling(base_molecule, probe_molecule) -> float:
    value = 0.0
    for base_atom in range(base_molecule.natm):
        for probe_atom in range(probe_molecule.natm):
            distance = np.linalg.norm(
                base_molecule.atom_coord(base_atom)
                - probe_molecule.atom_coord(probe_atom)
            )
            value += (
                base_molecule.atom_charge(base_atom)
                * probe_molecule.atom_charge(probe_atom)
                / distance
            )
    return float(value)


def evaluate_geometry(
    *,
    base,
    probe_positions: np.ndarray,
    probe_density: np.ndarray,
    density_values: dict[str, np.ndarray],
    glider_values: dict[str, np.ndarray],
    mace_values: dict[str, np.ndarray],
    glider_gradient_step_bohr: float = 1.0e-4,
    mace_grid_level: int = 5,
    include_components: bool = False,
) -> dict[str, float]:
    base_symbols = list(base.get_chemical_symbols())
    base_positions = np.asarray(base.positions)
    base_molecule = molecule(base_symbols, base_positions)
    probe_molecule = molecule(["O", "H", "H"], probe_positions)
    response_density = density_values["response_density_matrix"]
    if include_components:
        component_densities = density_values["component_density_matrices"]
        couplings = base_density_couplings(
            base_molecule,
            np.concatenate((component_densities, response_density[None]), axis=0),
            probe_molecule,
            probe_density,
        )
        qm_hartree = float(couplings[-1])
        linear_combination = float(couplings[0] - couplings[1] - couplings[2])
        linearity_error_hartree = qm_hartree - linear_combination
        full_frozen = float(couplings[0])
        full_frozen += electron_density_at_nuclei_coupling(
            probe_molecule, probe_density, base_molecule
        )
        full_frozen += nuclear_nuclear_coupling(base_molecule, probe_molecule)
    else:
        qm_hartree = float(
            base_density_couplings(
                base_molecule,
                response_density,
                probe_molecule,
                probe_density,
            )[0]
        )
        linearity_error_hartree = float("nan")
        full_frozen = float("nan")

    glider_hartree = point_multipole_coupling(
        base_positions,
        glider_values["predicted_charges_e"],
        glider_values["predicted_dipoles_e_bohr"],
        probe_molecule,
        probe_density,
        gradient_step_bohr=glider_gradient_step_bohr,
    )
    mace_hartree = gaussian_multipole_coupling(
        base_positions,
        mace_values["predicted_charges_e"],
        mace_values["predicted_dipoles_e_bohr"],
        probe_molecule,
        probe_density,
        grid_level=mace_grid_level,
    )
    response_dipole_e_bohr = (
        np.asarray(density_values["response_dipole_debye"])
        / DEBYE_PER_E_BOHR
    )
    origins = case_origins(base, int(base.info["n_solute_atoms"]))
    values_hartree = {
        "qm": qm_hartree,
        "glider": glider_hartree,
        "mace_polar_l": mace_hartree,
        "zero": 0.0,
    }
    for method, origin in origins.items():
        values_hartree[method] = point_multipole_coupling(
            origin[None],
            np.zeros(1),
            response_dipole_e_bohr[None],
            probe_molecule,
            probe_density,
            gradient_step_bohr=glider_gradient_step_bohr,
        )
    result = {
        f"{method}_energy_hartree": float(values_hartree[method])
        for method in METHODS
    }
    result.update(
        {
            f"{method}_energy_kcal_mol": float(
                values_hartree[method] * HARTREE_TO_KCAL_MOL
            )
            for method in METHODS
        }
    )
    result["full_frozen_electrostatic_hartree"] = full_frozen
    result["full_frozen_electrostatic_kcal_mol"] = (
        full_frozen * HARTREE_TO_KCAL_MOL
    )
    result["linearity_error_hartree"] = linearity_error_hartree
    return result


class Evaluator:
    def __init__(self, output: Path):
        self.output = output
        self.cache = output / "raw_evaluations"
        self.cache.mkdir(parents=True, exist_ok=True)
        self.bases = {
            str(base.info["case_id"]): base
            for base in read(
                RESULT_ROOT / "configurations/base_3water.extxyz", index=":"
            )
        }
        self.probes = {
            str(probe.info["case_id"]): probe
            for probe in read(
                RESULT_ROOT / "configurations/heldout_w4.extxyz", index=":"
            )
        }
        self.densities = load_file_map(
            RESULT_ROOT / "qm_densities",
            "density_registry.csv",
            "case_id",
            "density_file",
        )
        self.glider = load_predictions(
            ROOT / "results/liquid_bridge/predictions/glider"
        )
        self.mace = load_predictions(
            ROOT / "results/liquid_bridge/predictions/mace_polar_l"
        )
        variants = pd.read_csv(RESULT_ROOT / "probe_variants/variant_registry.csv")
        self.variants = {
            (str(row.case_id), str(row.variant_id)): RESULT_ROOT
            / "probe_variants"
            / str(row.density_file)
            for row in variants.itertuples()
        }
        if not (
            len(self.bases)
            == len(self.probes)
            == len(self.densities)
            == 24
        ):
            raise RuntimeError("the evaluator requires all 24 frozen cases")

    @staticmethod
    def arrays(path: Path) -> dict[str, np.ndarray]:
        with np.load(path, allow_pickle=False) as values:
            return {key: np.asarray(values[key]) for key in values.files}

    def case_values(self, case_id: str):
        base = self.bases[case_id]
        base_config_id = str(base.info["base_config_id"])
        return (
            base,
            self.arrays(self.densities[case_id]),
            self.arrays(self.glider[base_config_id]),
            self.arrays(self.mace[base_config_id]),
        )

    def evaluate(
        self,
        case_id: str,
        evaluation_id: str,
        positions: np.ndarray,
        density: np.ndarray,
        *,
        glider_step: float = 1.0e-4,
        mace_level: int = 5,
        include_components: bool = False,
    ) -> dict[str, object]:
        specification = (
            f"{case_id}::{evaluation_id}::g{glider_step:.8g}::m{mace_level}::"
            f"c{int(include_components)}"
        )
        tag = hashlib.sha256(specification.encode()).hexdigest()[:24]
        path = self.cache / f"{tag}.json"
        if path.exists():
            return json.loads(path.read_text())
        base, density_values, glider_values, mace_values = self.case_values(case_id)
        values = evaluate_geometry(
            base=base,
            probe_positions=positions,
            probe_density=density,
            density_values=density_values,
            glider_values=glider_values,
            mace_values=mace_values,
            glider_gradient_step_bohr=glider_step,
            mace_grid_level=mace_level,
            include_components=include_components,
        )
        record: dict[str, object] = {
            "case_id": case_id,
            "evaluation_id": evaluation_id,
            "glider_gradient_step_bohr": glider_step,
            "mace_grid_level": mace_level,
            **values,
        }
        path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        return record

    def initial_probe(self, case_id: str) -> tuple[np.ndarray, np.ndarray]:
        values = self.arrays(self.densities[case_id])
        return values["probe_positions_angstrom"], values["probe_density_matrix"]

    def variant_probe(self, case_id: str, variant_id: str):
        values = self.arrays(self.variants[(case_id, variant_id)])
        return values["positions_angstrom"], values["density_matrix"]


def surface_reconstruction_check(evaluator: Evaluator) -> pd.DataFrame:
    rows = []
    for case_id, base in sorted(evaluator.bases.items()):
        base_config_id = str(base.info["base_config_id"])
        prediction = evaluator.arrays(evaluator.glider[base_config_id])
        reconstructed = esp_from_sites(
            base.positions,
            prediction["points_angstrom"],
            prediction["predicted_charges_e"],
            prediction["predicted_dipoles_e_bohr"],
        )
        discrepancy = float(
            np.max(
                np.abs(reconstructed - prediction["predicted_esp_hartree_per_e"])
            )
        )
        if discrepancy > 1.0e-8:
            raise RuntimeError(f"GLIDER direct-field check failed: {case_id}")
        rows.append({"case_id": case_id, "max_abs_hartree_per_e": discrepancy})
    return pd.DataFrame(rows)


def run_energy(evaluator: Evaluator) -> None:
    rows = []
    checks = []
    for index, case_id in enumerate(sorted(evaluator.bases), 1):
        positions, density = evaluator.initial_probe(case_id)
        record = evaluator.evaluate(
            case_id,
            "initial",
            positions,
            density,
            include_components=True,
        )
        base = evaluator.bases[case_id]
        row = {
            "case_id": case_id,
            "base_config_id": str(base.info["base_config_id"]),
            "molecule_id": str(base.info["molecule_id"]),
            "source_frame_index": int(base.info["source_frame_index"]),
            **{key: value for key, value in record.items() if key.endswith("kcal_mol")},
            "linearity_error_hartree": record["linearity_error_hartree"],
        }
        rows.append(row)
        if "__liquid_00__" in case_id:
            for glider_step in (5.0e-5, 1.0e-4, 2.0e-4):
                check = evaluator.evaluate(
                    case_id,
                    f"initial_glider_step_{glider_step:g}",
                    positions,
                    density,
                    glider_step=glider_step,
                )
                checks.append(
                    {
                        "case_id": case_id,
                        "check": "glider_gradient_step",
                        "setting": glider_step,
                        "energy_kcal_mol": check["glider_energy_kcal_mol"],
                    }
                )
            for mace_level in (4, 5, 6):
                check = evaluator.evaluate(
                    case_id,
                    f"initial_mace_level_{mace_level}",
                    positions,
                    density,
                    mace_level=mace_level,
                )
                checks.append(
                    {
                        "case_id": case_id,
                        "check": "mace_quadrature_level",
                        "setting": mace_level,
                        "energy_kcal_mol": check["mace_polar_l_energy_kcal_mol"],
                    }
                )
        print(f"energy [{index}/24] {case_id}", flush=True)
    pd.DataFrame(rows).to_csv(evaluator.output / "experiment1_per_case.csv", index=False)
    pd.DataFrame(checks).to_csv(
        evaluator.output / "energy_numerical_convergence.csv", index=False
    )
    surface_reconstruction_check(evaluator).to_csv(
        evaluator.output / "glider_direct_field_check.csv", index=False
    )


def force_for_step(
    evaluator: Evaluator, case_id: str, step: float
) -> dict[str, np.ndarray]:
    initial_positions, density = evaluator.initial_probe(case_id)
    energies = {method: np.empty((3, 2)) for method in METHODS}
    for axis in range(3):
        for sign_index, sign in enumerate((-1.0, 1.0)):
            positions = initial_positions.copy()
            positions[:, axis] += sign * step
            record = evaluator.evaluate(
                case_id,
                f"force_{step:g}A_axis{axis}_{'minus' if sign < 0 else 'plus'}",
                positions,
                density,
            )
            for method in METHODS:
                energies[method][axis, sign_index] = record[
                    f"{method}_energy_kcal_mol"
                ]
    return {
        method: -(values[:, 1] - values[:, 0]) / (2.0 * step)
        for method, values in energies.items()
    }


def vector_rows(
    case_id: str, molecule_id: str, vectors: dict[str, np.ndarray]
) -> list[dict[str, object]]:
    reference = np.asarray(vectors["qm"])
    rows = []
    for method in METHODS:
        value = np.asarray(vectors[method])
        error = value - reference
        ref_norm = float(np.linalg.norm(reference))
        value_norm = float(np.linalg.norm(value))
        if ref_norm > 0 and value_norm > 0:
            cosine = float(np.clip(np.dot(value, reference) / (value_norm * ref_norm), -1, 1))
            angle = float(np.degrees(np.arccos(cosine)))
        else:
            cosine = float("nan")
            angle = float("nan")
        rows.append(
            {
                "case_id": case_id,
                "molecule_id": molecule_id,
                "method": method,
                "x": value[0],
                "y": value[1],
                "z": value[2],
                "magnitude": value_norm,
                "reference_magnitude": ref_norm,
                "vector_error_norm": float(np.linalg.norm(error)),
                "magnitude_error": value_norm - ref_norm,
                "cosine_similarity": cosine,
                "angular_error_degrees": angle,
            }
        )
    return rows


def run_force(evaluator: Evaluator) -> None:
    primary_rows = []
    convergence_rows = []
    for index, case_id in enumerate(sorted(evaluator.bases), 1):
        base = evaluator.bases[case_id]
        molecule_id = str(base.info["molecule_id"])
        primary = force_for_step(evaluator, case_id, 0.002)
        primary_rows.extend(vector_rows(case_id, molecule_id, primary))
        if "__liquid_00__" in case_id:
            for step in FORCE_STEPS:
                values = primary if step == 0.002 else force_for_step(
                    evaluator, case_id, step
                )
                for row in vector_rows(case_id, molecule_id, values):
                    convergence_rows.append({"step_angstrom": step, **row})
        print(f"force [{index}/24] {case_id}", flush=True)
    pd.DataFrame(primary_rows).to_csv(
        evaluator.output / "force_per_case.csv", index=False
    )
    pd.DataFrame(convergence_rows).to_csv(
        evaluator.output / "force_step_convergence.csv", index=False
    )


def torque_for_step(
    evaluator: Evaluator, case_id: str, step_degrees: float
) -> dict[str, np.ndarray]:
    energies = {method: np.empty((3, 2)) for method in METHODS}
    for axis in range(3):
        for sign_index, label in enumerate(("minus", "plus")):
            variant_id = f"torque_{step_degrees:g}deg_axis{axis}_{label}"
            positions, density = evaluator.variant_probe(case_id, variant_id)
            record = evaluator.evaluate(
                case_id, variant_id, positions, density
            )
            for method in METHODS:
                energies[method][axis, sign_index] = record[
                    f"{method}_energy_kcal_mol"
                ]
    step_radians = np.radians(step_degrees)
    return {
        method: -(values[:, 1] - values[:, 0]) / (2.0 * step_radians)
        for method, values in energies.items()
    }


def run_torque(evaluator: Evaluator) -> None:
    primary_rows = []
    convergence_rows = []
    for index, case_id in enumerate(sorted(evaluator.bases), 1):
        base = evaluator.bases[case_id]
        molecule_id = str(base.info["molecule_id"])
        primary = torque_for_step(evaluator, case_id, 0.1)
        primary_rows.extend(vector_rows(case_id, molecule_id, primary))
        if "__liquid_00__" in case_id:
            for step in TORQUE_STEPS:
                values = primary if step == 0.1 else torque_for_step(
                    evaluator, case_id, step
                )
                for row in vector_rows(case_id, molecule_id, values):
                    convergence_rows.append({"step_degrees": step, **row})
        print(f"torque [{index}/24] {case_id}", flush=True)
    pd.DataFrame(primary_rows).to_csv(
        evaluator.output / "torque_per_case.csv", index=False
    )
    pd.DataFrame(convergence_rows).to_csv(
        evaluator.output / "torque_step_convergence.csv", index=False
    )


def run_orientation(evaluator: Evaluator) -> None:
    rows = []
    case_ids = [
        case_id for case_id in sorted(evaluator.bases) if "__liquid_00__" in case_id
    ]
    for case_index, case_id in enumerate(case_ids, 1):
        base = evaluator.bases[case_id]
        molecule_id = str(base.info["molecule_id"])
        for rotation_index in range(24):
            rotation_id = f"O_{rotation_index:02d}"
            variant_id = f"orientation_{rotation_id}"
            positions, density = evaluator.variant_probe(case_id, variant_id)
            variant_values = evaluator.arrays(evaluator.variants[(case_id, variant_id)])
            record = evaluator.evaluate(case_id, variant_id, positions, density)
            matrix = variant_values["rotation_matrix"]
            quaternion = Rotation.from_matrix(matrix).as_quat()
            for method in METHODS:
                rows.append(
                    {
                        "case_id": case_id,
                        "molecule_id": molecule_id,
                        "rotation_id": rotation_id,
                        "rotation_index": rotation_index,
                        "quaternion_x": quaternion[0],
                        "quaternion_y": quaternion[1],
                        "quaternion_z": quaternion[2],
                        "quaternion_w": quaternion[3],
                        "method": method,
                        "energy_kcal_mol": record[
                            f"{method}_energy_kcal_mol"
                        ],
                    }
                )
        print(f"orientation [{case_index}/6] {case_id}", flush=True)
    pd.DataFrame(rows).to_csv(
        evaluator.output / "orientation_profiles.csv", index=False
    )


def write_manifest(output: Path) -> None:
    files = [
        path
        for path in output.iterdir()
        if path.is_file() and path.name != "evaluation_manifest.json"
    ]
    manifest = {
        "experiment": "post hoc downstream response coupling",
        "protocol": "results/posthoc_downstream_response/protocol_freeze.json",
        "outputs": {path.name: sha256(path) for path in sorted(files)},
        "methods": list(METHODS),
        "manuscript_modified": False,
        "prospective_claim": False,
    }
    (output / "evaluation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=RESULT_ROOT / "evaluated"
    )
    parser.add_argument(
        "--stage",
        choices=("energy", "force", "torque", "orientation", "all"),
        default="all",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    evaluator = Evaluator(args.output)
    if args.stage in ("energy", "all"):
        run_energy(evaluator)
    if args.stage in ("force", "all"):
        run_force(evaluator)
    if args.stage in ("torque", "all"):
        run_torque(evaluator)
    if args.stage in ("orientation", "all"):
        run_orientation(evaluator)
    write_manifest(args.output)


if __name__ == "__main__":
    main()
