#!/usr/bin/env python3
"""Build the frozen fresh sulfur/polar development stress set."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import write


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "manybody_completion/src"))
from build_clusters import cooperative_chain, reorient  # noqa: E402
sys.path.insert(0, str(WORKSPACE / "complete_interaction_hamiltonian/src"))
from build_configurations import (  # noqa: E402
    build_cluster,
    conformers,
    direction_outward,
    make_atoms,
    place_at_minimum_distance,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    molecules = pd.read_csv(ROOT / "config/development_molecules.csv")
    seed = 2026081307
    frames = []
    rows = []
    for molecule_index, row in molecules.iterrows():
        rng = np.random.default_rng(seed + 104729 * molecule_index)
        molecule, coordinates = conformers(row.smiles, seed + molecule_index)
        solute = coordinates[0]
        heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
        hetero = [
            atom.GetIdx()
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() in (7, 8, 9, 15, 16, 17)
        ]
        anchors = hetero + [index for index in heavy if index not in hetero]
        radial = direction_outward(solute, anchors[0], rng)
        one = [
            place_at_minimum_distance(
                solute, anchors[0], radial, 2.15, "water_donor", rng
            )
        ]
        two = build_cluster(solute, anchors, 2, rng, "equilibrium")
        equilibrium = build_cluster(solute, anchors, 3, rng, "equilibrium")
        specifications = [
            ("control_1water", "zero", one, "pair_zero"),
            ("control_2water", "zero", two, "threebody_zero"),
            ("cluster_3water", "equilibrium", equilibrium, "equilibrium"),
            (
                "cluster_3water",
                "cooperative",
                cooperative_chain(solute, anchors, 3, rng),
                "cooperative",
            ),
            (
                "cluster_3water",
                "orientation",
                reorient(equilibrium, rng),
                "orientation",
            ),
            (
                "cluster_3water",
                "compressed",
                cooperative_chain(solute, anchors, 3, rng, compressed=True),
                "compressed",
            ),
        ]
        for category, motif, waters, regime in specifications:
            config_id = f"development__{row.molecule_id}__{category}__{motif}"
            metadata = {
                "config_id": config_id,
                "molecule_id": row.molecule_id,
                "smiles": row.smiles,
                "split": "development",
                "chemical_role": row.chemical_role,
                "category": category,
                "motif": motif,
                "regime": regime,
                "conformer_index": 0,
            }
            atoms = make_atoms(molecule, solute, waters, metadata)
            frames.append(atoms)
            rows.append(
                {
                    **metadata,
                    "n_atoms": len(atoms),
                    "n_solute_atoms": len(solute),
                    "n_waters": len(waters),
                    "minimum_solute_water_distance_A": float(
                        np.linalg.norm(
                            atoms.positions[: len(solute), None]
                            - atoms.positions[None, len(solute) :],
                            axis=2,
                        ).min()
                    ),
                }
            )
    order = np.argsort([str(frame.info["config_id"]) for frame in frames])
    frames = [frames[index] for index in order]
    table = pd.DataFrame(rows).sort_values("config_id").reset_index(drop=True)
    output = ROOT / "data/development_stress"
    output.mkdir(parents=True, exist_ok=True)
    xyz = output / "configurations.extxyz"
    csv = output / "configuration_registry.csv"
    write(xyz, frames, format="extxyz")
    table.to_csv(csv, index=False)
    manifest = {
        "seed": seed,
        "n_molecules": int(table.molecule_id.nunique()),
        "n_configurations": len(table),
        "counts_by_n_waters": {
            str(key): int(value) for key, value in table.groupby("n_waters").size().items()
        },
        "counts_by_regime": {
            str(key): int(value) for key, value in table.groupby("regime").size().items()
        },
        "configuration_sha256": sha256(xyz),
        "registry_sha256": sha256(csv),
        "experimental_hydration_targets_accessed": False,
        "prospective_panel_accessed": False,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
