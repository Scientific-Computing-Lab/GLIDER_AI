"""Shared physical and numerical routines for the held-out-water experiment."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
from pyscf import df, dft, gto
from scipy.special import erf

BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_TO_KCAL_MOL = 627.5094740631
DEBYE_PER_E_BOHR = 2.541746473
MACE_GAUSSIAN_SIGMA_ANGSTROM = 1.5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def molecule(
    symbols: list[str],
    positions_angstrom: np.ndarray,
    real_atoms: np.ndarray | None = None,
) -> gto.Mole:
    positions = np.asarray(positions_angstrom, dtype=float)
    if real_atoms is None:
        real_atoms = np.ones(len(symbols), dtype=bool)
    atoms = [
        (symbol if bool(real_atoms[index]) else f"ghost-{symbol}", tuple(position))
        for index, (symbol, position) in enumerate(zip(symbols, positions, strict=True))
    ]
    return gto.M(
        atom=atoms,
        basis="def2-tzvpd",
        unit="Angstrom",
        charge=0,
        spin=0,
        verbose=0,
        max_memory=22000,
    )


def combined_molecule(
    base_symbols: list[str],
    base_positions_angstrom: np.ndarray,
    probe_symbols: list[str],
    probe_positions_angstrom: np.ndarray,
) -> gto.Mole:
    symbols = [*base_symbols, *probe_symbols]
    positions = np.vstack((base_positions_angstrom, probe_positions_angstrom))
    return molecule(symbols, positions)


def cross_electron_coulomb_df(
    base_molecule: gto.Mole,
    base_density: np.ndarray,
    probe_molecule: gto.Mole,
    probe_density: np.ndarray,
) -> float:
    """Density-fitted electron-electron Coulomb coupling in Hartree."""
    combo = combined_molecule(
        [base_molecule.atom_symbol(i) for i in range(base_molecule.natm)],
        base_molecule.atom_coords(unit="Angstrom"),
        [probe_molecule.atom_symbol(i) for i in range(probe_molecule.natm)],
        probe_molecule.atom_coords(unit="Angstrom"),
    )
    n_base_ao = base_molecule.nao_nr()
    n_probe_ao = probe_molecule.nao_nr()
    if combo.nao_nr() != n_base_ao + n_probe_ao:
        raise RuntimeError("combined AO ordering is not a base/probe concatenation")
    density = np.zeros((combo.nao_nr(), combo.nao_nr()))
    density[n_base_ao:, n_base_ao:] = np.asarray(probe_density)
    coulomb = df.DF(combo).get_jk(density, with_k=False)[0]
    return float(
        np.einsum(
            "ij,ji->",
            np.asarray(base_density),
            coulomb[:n_base_ao, :n_base_ao],
        )
    )


def probe_interaction_operator_on_base(
    base_molecule: gto.Mole,
    probe_molecule: gto.Mole,
    probe_density: np.ndarray,
) -> np.ndarray:
    """AO operator coupling a base electron density to the complete probe."""
    combo = combined_molecule(
        [base_molecule.atom_symbol(i) for i in range(base_molecule.natm)],
        base_molecule.atom_coords(unit="Angstrom"),
        [probe_molecule.atom_symbol(i) for i in range(probe_molecule.natm)],
        probe_molecule.atom_coords(unit="Angstrom"),
    )
    n_base_ao = base_molecule.nao_nr()
    n_probe_ao = probe_molecule.nao_nr()
    if combo.nao_nr() != n_base_ao + n_probe_ao:
        raise RuntimeError("combined AO ordering is not a base/probe concatenation")
    density = np.zeros((combo.nao_nr(), combo.nao_nr()))
    density[n_base_ao:, n_base_ao:] = np.asarray(probe_density)
    coulomb = df.DF(combo).get_jk(density, with_k=False)[0]
    operator = np.asarray(coulomb[:n_base_ao, :n_base_ao]).copy()
    for atom in range(probe_molecule.natm):
        with base_molecule.with_rinv_origin(probe_molecule.atom_coord(atom)):
            inverse_distance = base_molecule.intor("int1e_rinv")
        operator -= probe_molecule.atom_charge(atom) * inverse_distance
    return operator


def base_density_couplings(
    base_molecule: gto.Mole,
    base_densities: np.ndarray,
    probe_molecule: gto.Mole,
    probe_density: np.ndarray,
) -> np.ndarray:
    """Couple one or more base electron densities to the complete probe."""
    densities = np.asarray(base_densities)
    if densities.ndim == 2:
        densities = densities[None, :, :]
    operator = probe_interaction_operator_on_base(
        base_molecule, probe_molecule, probe_density
    )
    return np.einsum("nij,ji->n", densities, operator)


def electron_density_at_nuclei_coupling(
    electron_molecule: gto.Mole,
    electron_density: np.ndarray,
    nuclear_molecule: gto.Mole,
) -> float:
    """Attraction of positive nuclei to an electron number density."""
    value = 0.0
    for atom in range(nuclear_molecule.natm):
        with electron_molecule.with_rinv_origin(nuclear_molecule.atom_coord(atom)):
            inverse_distance = electron_molecule.intor("int1e_rinv")
        value -= nuclear_molecule.atom_charge(atom) * np.einsum(
            "ij,ji->", electron_density, inverse_distance
        )
    return float(value)


def response_coupling(
    base_molecule: gto.Mole,
    response_density: np.ndarray,
    probe_molecule: gto.Mole,
    probe_density: np.ndarray,
) -> float:
    return float(
        base_density_couplings(
            base_molecule, response_density, probe_molecule, probe_density
        )[0]
    )


def full_frozen_electrostatic_coupling(
    base_molecule: gto.Mole,
    base_density: np.ndarray,
    probe_molecule: gto.Mole,
    probe_density: np.ndarray,
) -> float:
    value = float(
        base_density_couplings(
            base_molecule, base_density, probe_molecule, probe_density
        )[0]
    )
    value += electron_density_at_nuclei_coupling(
        probe_molecule, probe_density, base_molecule
    )
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


def molecular_electrostatic_potential(
    molecule_value: gto.Mole,
    density: np.ndarray,
    points_angstrom: np.ndarray,
    chunk_size: int = 64,
) -> np.ndarray:
    """Frozen molecular ESP in Hartree/e at arbitrary points."""
    coordinates = np.asarray(points_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    nuclear = np.zeros(len(coordinates))
    for atom in range(molecule_value.natm):
        displacement = molecule_value.atom_coord(atom) - coordinates
        nuclear += molecule_value.atom_charge(atom) / np.linalg.norm(
            displacement, axis=1
        )
    electronic = np.empty(len(coordinates))
    for begin in range(0, len(coordinates), chunk_size):
        end = min(begin + chunk_size, len(coordinates))
        fake = gto.fakemol_for_charges(coordinates[begin:end])
        integrals = df.incore.aux_e2(molecule_value, fake)
        electronic[begin:end] = np.einsum(
            "ijp,ij->p", integrals, np.asarray(density)
        )
    return nuclear - electronic


def point_multipole_coupling(
    site_positions_angstrom: np.ndarray,
    charges_e: np.ndarray,
    dipoles_e_bohr: np.ndarray,
    probe_molecule: gto.Mole,
    probe_density: np.ndarray,
    gradient_step_bohr: float = 1.0e-4,
) -> float:
    """Couple point q/mu response sites to a frozen molecular probe."""
    positions = np.asarray(site_positions_angstrom, dtype=float)
    charges = np.asarray(charges_e, dtype=float)
    dipoles = np.asarray(dipoles_e_bohr, dtype=float)
    if len(positions) != len(charges) or dipoles.shape != (len(positions), 3):
        raise ValueError("site positions, charges, and dipoles differ")
    potential = molecular_electrostatic_potential(
        probe_molecule, probe_density, positions
    )
    gradient = np.empty((len(positions), 3))
    displacement = gradient_step_bohr * BOHR_TO_ANGSTROM
    for axis in range(3):
        plus = positions.copy()
        minus = positions.copy()
        plus[:, axis] += displacement
        minus[:, axis] -= displacement
        gradient[:, axis] = (
            molecular_electrostatic_potential(
                probe_molecule, probe_density, plus
            )
            - molecular_electrostatic_potential(
                probe_molecule, probe_density, minus
            )
        ) / (2.0 * gradient_step_bohr)
    return float(charges @ potential + np.einsum("ij,ij->", dipoles, gradient))


def gaussian_multipole_potential(
    site_positions_angstrom: np.ndarray,
    points_angstrom: np.ndarray,
    charges_e: np.ndarray,
    dipoles_e_bohr: np.ndarray,
    sigma_angstrom: float = MACE_GAUSSIAN_SIGMA_ANGSTROM,
) -> np.ndarray:
    """Official MACE-POLAR open-boundary normalized Gaussian l<=1 potential."""
    positions = np.asarray(site_positions_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    points = np.asarray(points_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    sigma = sigma_angstrom / BOHR_TO_ANGSTROM
    displacement = points[:, None, :] - positions[None, :, :]
    distance = np.linalg.norm(displacement, axis=-1)
    if np.any(distance < 1.0e-10):
        raise ValueError("Gaussian potential point coincides with a site")
    argument = distance / (math.sqrt(2.0) * sigma)
    gaussian = np.exp(-0.5 * (distance / sigma) ** 2)
    monopole_kernel = erf(argument) / distance
    radial_dipole_kernel = erf(argument) / distance**3 - math.sqrt(
        2.0 / math.pi
    ) * gaussian / (sigma * distance**2)
    monopoles = monopole_kernel @ np.asarray(charges_e, dtype=float)
    dot = np.einsum(
        "pnc,nc->pn", displacement, np.asarray(dipoles_e_bohr, dtype=float)
    )
    return monopoles + np.sum(dot * radial_dipole_kernel, axis=1)


def gaussian_multipole_coupling(
    site_positions_angstrom: np.ndarray,
    charges_e: np.ndarray,
    dipoles_e_bohr: np.ndarray,
    probe_molecule: gto.Mole,
    probe_density: np.ndarray,
    grid_level: int = 5,
) -> float:
    """Integrate the official MACE Gaussian response field over frozen W4."""
    grids = dft.gen_grid.Grids(probe_molecule)
    grids.level = int(grid_level)
    grids.build()
    ao = dft.numint.eval_ao(probe_molecule, grids.coords)
    rho = dft.numint.eval_rho(probe_molecule, ao, np.asarray(probe_density))
    electronic_potential = gaussian_multipole_potential(
        site_positions_angstrom,
        grids.coords * BOHR_TO_ANGSTROM,
        charges_e,
        dipoles_e_bohr,
    )
    nuclear_potential = gaussian_multipole_potential(
        site_positions_angstrom,
        probe_molecule.atom_coords(unit="Angstrom"),
        charges_e,
        dipoles_e_bohr,
    )
    nuclear = sum(
        probe_molecule.atom_charge(atom) * nuclear_potential[atom]
        for atom in range(probe_molecule.natm)
    )
    electronic = np.dot(grids.weights * rho, electronic_potential)
    return float(nuclear - electronic)
