#!/usr/bin/env python3
"""Score the one short inter-water contact as an isolated frozen water dimer.

This is a geometry diagnostic, not a response reference or a replacement for
the frozen panel score. It uses the same CP energy convention and QM level as
the short solute--water energy audit.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pyscf
from ase.io import read

from acquire_nonwater_contact_qm_cpu import COMPONENTS, solve


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "experiments" / "water_contact_audit"
SOURCE = ROOT / "experiments" / "panel_1" / "geometries" / "configurations.extxyz"
HARTREE_TO_KCAL_MOL = 627.5094740631


def main() -> None:
    pyscf.lib.num_threads(4)
    with (AUDIT / "short_water_water_contacts.csv").open(newline="") as stream:
        flagged = list(csv.DictReader(stream))
    if len(flagged) != 1:
        raise ValueError(f"Expected one flagged inter-water pair, found {len(flagged)}")
    row = flagged[0]
    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    output = AUDIT / "short_water_pair_cp.json"
    if output.exists():
        old = json.loads(output.read_text())
        if old["source_sha256"] == source_hash and old["config_id"] == row["config_id"]:
            print(f"Verified existing result: {output}")
            return
        raise ValueError(f"Stale output: {output}")
    atoms = next(a for a in read(SOURCE, index=":") if a.info["config_id"] == row["config_id"])
    n_solute = int(atoms.info["n_solute_atoms"])
    water_a, water_b = int(row["water_fragment_a"]), int(row["water_fragment_b"])
    indices = [
        n_solute + 3 * (fragment - 1) + offset
        for fragment in (water_a, water_b)
        for offset in range(3)
    ]
    dimer = atoms[indices]
    if dimer.get_chemical_symbols() != ["O", "H", "H", "O", "H", "H"]:
        raise ValueError("Expected two whole water molecules in OHH order")
    masks = (
        np.ones(6, dtype=bool),
        np.array([True] * 3 + [False] * 3),
        np.array([False] * 3 + [True] * 3),
    )
    energies = {}
    attempts = {}
    for name, mask in zip(COMPONENTS, masks):
        mean_field, _, record = solve(dimer, mask, name)
        energies[name] = float(mean_field.e_tot)
        attempts[name] = record
    cp = energies[COMPONENTS[0]] - energies[COMPONENTS[1]] - energies[COMPONENTS[2]]
    result = {
        "config_id": row["config_id"],
        "water_fragments": [water_a, water_b],
        "minimum_interwater_distance_A": float(row["minimum_water_water_distance_A"]),
        "cp_pair_energy_kcal_mol": cp * HARTREE_TO_KCAL_MOL,
        "component_energies_hartree": energies,
        "attempts": attempts,
        "source_sha256": source_hash,
        "method": "DF-RKS omegaB97X-D3(BJ)/def2-TZVPD grid4 SCF 1e-10",
        "scope": "Isolated frozen pair of waters; the third water and solute are absent",
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"CP two-water pair energy: {result['cp_pair_energy_kcal_mol']:.3f} kcal/mol")


if __name__ == "__main__":
    main()
