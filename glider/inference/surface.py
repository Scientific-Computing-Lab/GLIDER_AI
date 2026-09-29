"""Deterministic ESP probe surface used by the frozen response campaign."""

from __future__ import annotations

import math

import numpy as np

VDW_ANGSTROM = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}


def fibonacci_sphere(n: int = 26) -> np.ndarray:
    """Return the frozen deterministic unit-sphere directions."""
    index = np.arange(n, dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    z = 1.0 - 2.0 * (index + 0.5) / n
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = golden * index
    return np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))


def molecular_surface(symbols: list[str], positions_angstrom: np.ndarray) -> np.ndarray:
    """Build the two-shell, neighbour-filtered probe surface exactly."""
    positions = np.asarray(positions_angstrom, dtype=float)
    directions = fibonacci_sphere(26)
    radii = np.asarray([VDW_ANGSTROM[symbol] for symbol in symbols])
    points: list[np.ndarray] = []
    for atom, (center, radius) in enumerate(zip(positions, radii)):
        for scale in (1.40, 1.80):
            candidates = center[None, :] + scale * radius * directions
            for point in candidates:
                distances = np.linalg.norm(positions - point[None, :], axis=1)
                mask = np.arange(len(symbols)) != atom
                if np.any(distances[mask] < 1.05 * radii[mask]):
                    continue
                points.append(point)
    if not points:
        raise RuntimeError("empty molecular surface")
    array = np.asarray(points)
    _, keep = np.unique(np.round(array, 5), axis=0, return_index=True)
    return array[np.sort(keep)]
