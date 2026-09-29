"""Physical reconstruction used by the frozen GLIDER response representation."""

from __future__ import annotations

import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903
DEBYE_PER_E_BOHR = 2.541746473


def esp_from_sites(
    positions_angstrom: np.ndarray,
    points_angstrom: np.ndarray,
    charges_e: np.ndarray,
    dipoles_e_bohr: np.ndarray,
) -> np.ndarray:
    """Open-boundary point monopole/dipole response ESP in Hartree/e."""
    positions = np.asarray(positions_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    points = np.asarray(points_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    charges = np.asarray(charges_e, dtype=float)
    dipoles = np.asarray(dipoles_e_bohr, dtype=float)
    displacement = points[:, None, :] - positions[None, :, :]
    distance = np.linalg.norm(displacement, axis=2)
    if np.any(distance == 0):
        raise ValueError("ESP probe coincides with an atomic centre")
    monopoles = (charges[None, :] / distance).sum(axis=1)
    dipole_term = (dipoles[None, :, :] * displacement / distance[:, :, None] ** 3).sum(axis=(1, 2))
    return monopoles + dipole_term


def dipole_from_sites(
    positions_angstrom: np.ndarray,
    charges_e: np.ndarray,
    dipoles_e_bohr: np.ndarray,
) -> np.ndarray:
    """Molecular response dipole in Debye."""
    positions_bohr = np.asarray(positions_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    charges = np.asarray(charges_e, dtype=float)
    dipoles = np.asarray(dipoles_e_bohr, dtype=float)
    return (
        (charges[:, None] * positions_bohr).sum(axis=0) + dipoles.sum(axis=0)
    ) * DEBYE_PER_E_BOHR
