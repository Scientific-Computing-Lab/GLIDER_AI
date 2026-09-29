#!/usr/bin/env python3
"""Save the frozen MINAO solute guess for same-Hamiltonian SCF continuation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from ase.io import read

WORKSPACE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(WORKSPACE))
from direct_response_discovery.src.acquire_response_qm import make_scf  # noqa: E402
from prospective_response_sota.src.acquire_cp_response import molecule_with_ghosts  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    atoms = read(args.configurations, args.index)
    nsolute = int(atoms.info["n_solute_atoms"])
    mask = np.arange(len(atoms)) < nsolute
    molecule = molecule_with_ghosts(atoms.get_chemical_symbols(), atoms.positions, mask)
    mean_field = make_scf(molecule)
    density = mean_field.get_init_guess(molecule, key="minao")
    if hasattr(density, "get"):
        density = density.get()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.output, np.asarray(density))
    print(atoms.info["config_id"], args.output)


if __name__ == "__main__":
    main()
