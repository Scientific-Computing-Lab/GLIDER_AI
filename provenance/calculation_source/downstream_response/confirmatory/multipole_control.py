"""Exact global response-multipole control for the confirmatory water probes."""

from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from pyscf import dft

HERE = Path(__file__).resolve().parent
COMMON = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(COMMON))
from common import (  # noqa: E402
    BOHR_TO_ANGSTROM,
    DEBYE_PER_E_BOHR,
    HARTREE_TO_KCAL_MOL,
    molecular_electrostatic_potential,
    molecule,
    point_multipole_coupling,
)

ROOT = Path(__file__).resolve().parents[4]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response_confirmatory"
CONTROL_ROOT = RESULT_ROOT / "multipole_control"

ORIGINS = (
    "base_com",
    "solute_com",
    "nuclear_center",
)
ORIGIN_TO_EXISTING_DIPOLE = {
    "base_com": "exact_dipole_base_com",
    "solute_com": "exact_dipole_solute_com",
    "nuclear_center": "exact_dipole_nuclear_center",
}
HIERARCHIES = (
    "exact_dipole",
    "exact_dipole_quadrupole",
    "exact_dipole_quadrupole_octupole",
)


def arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as values:
        return {key: np.asarray(values[key]) for key in values.files}


def load_file_map(
    directory: Path, registry: str, identifier: str, filename: str
) -> dict[str, Path]:
    table = pd.read_csv(directory / registry)
    return {
        str(getattr(row, identifier)): directory / str(getattr(row, filename))
        for row in table.itertuples()
    }


def centre(positions: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.average(np.asarray(positions), axis=0, weights=np.asarray(weights))


def case_origins(base) -> dict[str, np.ndarray]:
    positions = np.asarray(base.positions)
    masses = np.asarray(base.get_masses())
    charges = np.asarray(base.get_atomic_numbers(), dtype=float)
    n_solute = int(base.info["n_solute_atoms"])
    return {
        "base_com": centre(positions, masses),
        "solute_com": centre(positions[:n_solute], masses[:n_solute]),
        "nuclear_center": centre(positions, charges),
    }


def symmetrize_rank3(values: np.ndarray) -> np.ndarray:
    permutations = set(itertools.permutations((0, 1, 2)))
    return sum(np.transpose(values, axes) for axes in permutations) / len(permutations)


def raw_charge_moments(
    base_molecule, response_density: np.ndarray, origin_angstrom: np.ndarray
) -> dict[str, np.ndarray | float]:
    """Physical electronic response-charge moments through rank three."""
    origin_bohr = np.asarray(origin_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    density = np.asarray(response_density, dtype=float)
    overlap = base_molecule.intor("int1e_ovlp")
    with base_molecule.with_common_origin(origin_bohr):
        r = base_molecule.intor("int1e_r", comp=3)
        rr = base_molecule.intor("int1e_rr", comp=9).reshape(
            3, 3, base_molecule.nao_nr(), base_molecule.nao_nr()
        )
        rrr = base_molecule.intor("int1e_rrr", comp=27).reshape(
            3, 3, 3, base_molecule.nao_nr(), base_molecule.nao_nr()
        )
    # The stored matrix is a positive electron-number density difference.
    # Nuclear terms cancel, so physical response-charge moments carry a minus.
    charge = -float(np.einsum("ij,ji->", density, overlap))
    dipole = -np.einsum("aij,ji->a", r, density)
    second = -np.einsum("abij,ji->ab", rr, density)
    third = -np.einsum("abcij,ji->abc", rrr, density)
    second = 0.5 * (second + second.T)
    third = symmetrize_rank3(third)
    return {
        "charge": charge,
        "dipole": dipole,
        "second": second,
        "third": third,
    }


def irreducible_moments(raw: dict[str, np.ndarray | float]) -> dict[str, np.ndarray]:
    second = np.asarray(raw["second"])
    third = np.asarray(raw["third"])
    quadrupole = 3.0 * second - np.eye(3) * np.trace(second)
    trace_vector = np.einsum("iaa->i", third)
    octupole = 5.0 * third.copy()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                octupole[i, j, k] -= (
                    (float(i == j) * trace_vector[k])
                    + (float(i == k) * trace_vector[j])
                    + (float(j == k) * trace_vector[i])
                )
    return {"quadrupole": quadrupole, "octupole": symmetrize_rank3(octupole)}


def translate_raw_moments(
    raw: dict[str, np.ndarray | float], displacement_bohr: np.ndarray
) -> dict[str, np.ndarray | float]:
    """Translate moments from O to O'=O+a, with x'=x-a."""
    a = np.asarray(displacement_bohr, dtype=float)
    q = float(raw["charge"])
    p = np.asarray(raw["dipole"])
    second = np.asarray(raw["second"])
    third = np.asarray(raw["third"])
    translated_p = p - q * a
    translated_second = (
        second
        - np.einsum("i,j->ij", a, p)
        - np.einsum("j,i->ij", a, p)
        + q * np.einsum("i,j->ij", a, a)
    )
    translated_third = third.copy()
    translated_third -= np.einsum("i,jk->ijk", a, second)
    translated_third -= np.einsum("j,ik->ijk", a, second)
    translated_third -= np.einsum("k,ij->ijk", a, second)
    translated_third += np.einsum("i,j,k->ijk", a, a, p)
    translated_third += np.einsum("i,k,j->ijk", a, a, p)
    translated_third += np.einsum("j,k,i->ijk", a, a, p)
    translated_third -= q * np.einsum("i,j,k->ijk", a, a, a)
    return {
        "charge": q,
        "dipole": translated_p,
        "second": 0.5 * (translated_second + translated_second.T),
        "third": symmetrize_rank3(translated_third),
    }


def finite_difference_weights(derivative: int) -> np.ndarray:
    nodes = np.arange(-3, 4, dtype=float)
    matrix = np.vstack([nodes**power for power in range(7)])
    target = np.zeros(7)
    target[derivative] = float(math.factorial(derivative))
    return np.linalg.solve(matrix, target)


FD_WEIGHTS = {order: finite_difference_weights(order) for order in range(4)}


def probe_potential_derivatives(
    probe_molecule,
    probe_density: np.ndarray,
    origins_angstrom: dict[str, np.ndarray],
    step_bohr: float,
) -> dict[str, dict[str, np.ndarray]]:
    """Full frozen-water ESP derivatives through rank three at each origin."""
    origin_names = list(origins_angstrom)
    offsets = np.arange(-3, 4, dtype=float)
    lattice = np.stack(np.meshgrid(offsets, offsets, offsets, indexing="ij"), axis=-1)
    lattice_angstrom = lattice * step_bohr * BOHR_TO_ANGSTROM
    points = np.concatenate(
        [np.asarray(origins_angstrom[name])[None, None, None, :] + lattice_angstrom for name in origin_names],
        axis=0,
    ).reshape(-1, 3)
    potential = molecular_electrostatic_potential(
        probe_molecule, probe_density, points
    ).reshape(len(origin_names), 7, 7, 7)
    output: dict[str, dict[str, np.ndarray]] = {}
    for origin_index, name in enumerate(origin_names):
        grid = potential[origin_index]
        hessian = np.empty((3, 3))
        third = np.empty((3, 3, 3))
        for i in range(3):
            for j in range(3):
                derivative_counts = [0, 0, 0]
                derivative_counts[i] += 1
                derivative_counts[j] += 1
                hessian[i, j] = np.einsum(
                    "abc,a,b,c->",
                    grid,
                    FD_WEIGHTS[derivative_counts[0]],
                    FD_WEIGHTS[derivative_counts[1]],
                    FD_WEIGHTS[derivative_counts[2]],
                ) / step_bohr**2
                for k in range(3):
                    counts = derivative_counts.copy()
                    counts[k] += 1
                    third[i, j, k] = np.einsum(
                        "abc,a,b,c->",
                        grid,
                        FD_WEIGHTS[counts[0]],
                        FD_WEIGHTS[counts[1]],
                        FD_WEIGHTS[counts[2]],
                    ) / step_bohr**3
        output[name] = {
            "hessian": 0.5 * (hessian + hessian.T),
            "third": symmetrize_rank3(third),
        }
    return output


def hierarchy_couplings(
    *,
    base,
    density_values: dict[str, np.ndarray],
    probe_positions_angstrom: np.ndarray,
    probe_density: np.ndarray,
    step_bohr: float = 0.01,
) -> dict[tuple[str, str], dict[str, float]]:
    origins = case_origins(base)
    base_molecule = molecule(list(base.get_chemical_symbols()), np.asarray(base.positions))
    probe_molecule = molecule(["O", "H", "H"], probe_positions_angstrom)
    derivatives = probe_potential_derivatives(
        probe_molecule, probe_density, origins, step_bohr
    )
    response_density = density_values["response_density_matrix"]
    # Retain the published exact QM dipole for exact order-1 reproduction.
    dipole = np.asarray(density_values["response_dipole_debye"]) / DEBYE_PER_E_BOHR
    output: dict[tuple[str, str], dict[str, float]] = {}
    for origin_name, origin in origins.items():
        raw = raw_charge_moments(base_molecule, response_density, origin)
        irreducible = irreducible_moments(raw)
        dipole_energy = point_multipole_coupling(
            origin[None],
            np.zeros(1),
            dipole[None],
            probe_molecule,
            probe_density,
            gradient_step_bohr=1.0e-4,
        )
        quadrupole_energy = float(
            np.einsum(
                "ij,ij->",
                irreducible["quadrupole"],
                derivatives[origin_name]["hessian"],
            )
            / 6.0
        )
        octupole_energy = float(
            np.einsum(
                "ijk,ijk->",
                irreducible["octupole"],
                derivatives[origin_name]["third"],
            )
            / 30.0
        )
        cumulative = {
            "exact_dipole": dipole_energy,
            "exact_dipole_quadrupole": dipole_energy + quadrupole_energy,
            "exact_dipole_quadrupole_octupole": dipole_energy
            + quadrupole_energy
            + octupole_energy,
        }
        for hierarchy, value in cumulative.items():
            output[(hierarchy, origin_name)] = {
                "energy_hartree": float(value),
                "energy_kcal_mol": float(value * HARTREE_TO_KCAL_MOL),
                "dipole_term_hartree": float(dipole_energy),
                "quadrupole_term_hartree": float(quadrupole_energy),
                "octupole_term_hartree": float(octupole_energy),
            }
    return output


class MultipoleData:
    def __init__(self) -> None:
        self.bases = {
            str(base.info["case_id"]): base
            for base in read(RESULT_ROOT / "configurations/base_3water.extxyz", index=":")
        }
        self.probes = {
            str(probe.info["probe_id"]): probe
            for probe in read(RESULT_ROOT / "configurations/outer_w4_w12.extxyz", index=":")
        }
        self.base_density_files = load_file_map(
            RESULT_ROOT / "qm/base_densities",
            "density_registry.csv",
            "case_id",
            "density_file",
        )
        self.probe_density_files = load_file_map(
            RESULT_ROOT / "qm/probe_densities",
            "density_registry.csv",
            "probe_id",
            "density_file",
        )
        variant_table = pd.read_csv(RESULT_ROOT / "qm/probe_variants/variant_registry.csv")
        self.variant_files = {
            (str(row.probe_id), str(row.variant_id)): RESULT_ROOT
            / "qm/probe_variants"
            / str(row.density_file)
            for row in variant_table.itertuples()
        }
        self._base_values: dict[str, dict[str, np.ndarray]] = {}
        self._probe_values: dict[str, dict[str, np.ndarray]] = {}
        self._variant_values: dict[tuple[str, str], dict[str, np.ndarray]] = {}
        if not (
            len(self.bases) == len(self.base_density_files) == 10
            and len(self.probes) == len(self.probe_density_files) == 90
        ):
            raise RuntimeError("multipole control requires 10 bases and 90 probes")

    def base_values(self, case_id: str) -> dict[str, np.ndarray]:
        if case_id not in self._base_values:
            self._base_values[case_id] = arrays(self.base_density_files[case_id])
        return self._base_values[case_id]

    def probe_values(self, probe_id: str) -> dict[str, np.ndarray]:
        if probe_id not in self._probe_values:
            self._probe_values[probe_id] = arrays(self.probe_density_files[probe_id])
        return self._probe_values[probe_id]

    def variant_values(self, probe_id: str, variant_id: str) -> dict[str, np.ndarray]:
        key = (probe_id, variant_id)
        if key not in self._variant_values:
            self._variant_values[key] = arrays(self.variant_files[key])
        return self._variant_values[key]

    def metadata(self, probe_id: str) -> dict[str, object]:
        probe = self.probes[probe_id]
        case_id = str(probe.info["case_id"])
        base = self.bases[case_id]
        return {
            "probe_id": probe_id,
            "case_id": case_id,
            "molecule_id": str(probe.info["molecule_id"]),
            "family": str(base.info["family"]),
            "water_rank": int(probe.info["water_rank"]),
            "source_water_index": int(probe.info["source_water_index"]),
            "oxygen_distance_A": float(probe.info["oxygen_distance_A"]),
        }

    def w4_ids(self) -> list[str]:
        return sorted(
            probe_id
            for probe_id, probe in self.probes.items()
            if int(probe.info["water_rank"]) == 4
        )


def source_extent_for_case(base, density_values: dict[str, np.ndarray]) -> tuple[dict[str, dict], dict]:
    base_molecule = molecule(list(base.get_chemical_symbols()), np.asarray(base.positions))
    grids = dft.gen_grid.Grids(base_molecule)
    grids.level = 4
    grids.build()
    numerical_integrator = dft.numint.NumInt()
    coordinate_blocks = []
    signed_weight_blocks = []
    absolute_weight_blocks = []
    density = density_values["response_density_matrix"]
    for ao, _, weights, coordinates in numerical_integrator.block_loop(
        base_molecule,
        grids,
        base_molecule.nao_nr(),
        deriv=0,
        max_memory=4000,
    ):
        rho_electron = numerical_integrator.eval_rho(base_molecule, ao, density)
        coordinate_blocks.append(np.asarray(coordinates))
        signed_weight_blocks.append(np.asarray(rho_electron) * np.asarray(weights))
        absolute_weight_blocks.append(np.abs(rho_electron) * np.asarray(weights))
    coordinates = np.concatenate(coordinate_blocks)
    signed_weights = np.concatenate(signed_weight_blocks)
    absolute_weights = np.concatenate(absolute_weight_blocks)
    signed_electron = float(np.sum(signed_weights))
    absolute_integral = float(np.sum(absolute_weights))
    origins = case_origins(base)
    values: dict[str, dict] = {}
    for name, origin_angstrom in origins.items():
        origin_bohr = origin_angstrom / BOHR_TO_ANGSTROM
        radii_bohr = np.linalg.norm(coordinates - origin_bohr, axis=1)
        order = np.argsort(radii_bohr)
        cumulative = np.cumsum(absolute_weights[order]) / absolute_integral
        quantiles = {}
        for fraction in (0.90, 0.95, 0.99):
            index = min(int(np.searchsorted(cumulative, fraction)), len(order) - 1)
            quantiles[f"r{int(fraction * 100)}_angstrom"] = float(
                radii_bohr[order[index]] * BOHR_TO_ANGSTROM
            )
        values[name] = {
            **quantiles,
            "max_base_nuclear_radius_angstrom": float(
                np.max(np.linalg.norm(np.asarray(base.positions) - origin_angstrom, axis=1))
            ),
        }
    diagnostics = {
        "grid_level": 4,
        "n_grid_points": int(len(coordinates)),
        "integrated_response_electron_number": signed_electron,
        "integrated_absolute_response_charge_e": absolute_integral,
    }
    return values, diagnostics
