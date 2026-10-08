#!/usr/bin/env python3
"""Make a label-blind contact panel from the original non-water starting poses.

Each of the 12 original solutes and three neighbour species has two starting
poses.  Optimize each with MMFF94s, fixing every solute atom and allowing the
neighbour to move.  Keep the converged pose with the lower MMFF interaction
energy.  No GLIDER prediction, response label, or QM energy is read here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import ase
import numpy as np
import rdkit
from ase.io import read, write
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D


ROOT = Path(__file__).resolve().parents[2]
NEIGHBOUR_SMILES = {"NH3": "N", "CH3OH": "CO", "CH3CN": "CC#N"}
MIN_HEAVY_ATOM_DISTANCE_A = 2.4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def optimize(atoms: ase.Atoms) -> tuple[ase.Atoms, dict[str, object]]:
    n_solute = int(atoms.info["n_solute_atoms"])
    species = str(atoms.info["perturbant"])
    solute = Chem.AddHs(Chem.MolFromSmiles(str(atoms.info["smiles"])))
    neighbour = Chem.AddHs(Chem.MolFromSmiles(NEIGHBOUR_SMILES[species]))
    molecule = Chem.CombineMols(solute, neighbour)
    Chem.SanitizeMol(molecule)
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    if symbols != atoms.get_chemical_symbols() or n_solute != solute.GetNumAtoms():
        raise ValueError(f"RDKit and archived atom orders differ for {atoms.info['config_id']}")
    conformer = Chem.Conformer(len(atoms))
    for index, xyz in enumerate(atoms.positions):
        conformer.SetAtomPosition(index, Point3D(*map(float, xyz)))
    molecule.AddConformer(conformer)
    parameters = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94s")
    if parameters is None:
        raise ValueError(f"MMFF94s parameters unavailable for {atoms.info['config_id']}")
    force_field = AllChem.MMFFGetMoleculeForceField(
        molecule, parameters, ignoreInterfragInteractions=False
    )
    for index in range(n_solute):
        force_field.AddFixedPoint(index)
    status = int(force_field.Minimize(maxIts=1000))
    positions = np.asarray(force_field.Positions(), dtype=float).reshape(-1, 3)
    for index, xyz in enumerate(positions):
        molecule.GetConformer().SetAtomPosition(index, Point3D(*map(float, xyz)))
    internal_only = AllChem.MMFFGetMoleculeForceField(
        molecule, parameters, ignoreInterfragInteractions=True
    )
    interaction = float(force_field.CalcEnergy() - internal_only.CalcEnergy())
    heavy_solute = np.flatnonzero(np.asarray(atoms.numbers[:n_solute]) != 1)
    heavy_neighbour = np.flatnonzero(np.asarray(atoms.numbers[n_solute:]) != 1) + n_solute
    heavy_distances = np.linalg.norm(
        positions[heavy_solute, None, :] - positions[None, heavy_neighbour, :], axis=2
    )
    min_heavy = float(np.min(heavy_distances))
    min_any = float(
        np.min(
            np.linalg.norm(
                positions[:n_solute, None, :] - positions[None, n_solute:, :], axis=2
            )
        )
    )
    result = atoms.copy()
    result.positions = positions
    metadata: dict[str, object] = {
        "source_config_id": str(atoms.info["config_id"]),
        "molecule_id": str(atoms.info["molecule_id"]),
        "perturbant": species,
        "source_regime": str(atoms.info["regime"]),
        "mmff_status": status,
        "mmff_interaction_kcal_mol": interaction,
        "minimum_heavy_atom_distance_A": min_heavy,
        "minimum_any_atom_distance_A": min_any,
    }
    return result, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path,
        default=ROOT / "experiments/nonwater/geometries/configurations.extxyz",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "build/nonwater_contact_panel",
    )
    args = parser.parse_args()
    frames = read(args.source, index=":")
    if len(frames) != 72:
        raise ValueError("Expected all 72 original source configurations")
    candidates: dict[tuple[str, str], list[tuple[ase.Atoms, dict[str, object]]]] = defaultdict(list)
    for atoms in frames:
        result, metadata = optimize(atoms)
        key = (str(metadata["molecule_id"]), str(metadata["perturbant"]))
        candidates[key].append((result, metadata))
    if len(candidates) != 36 or any(len(values) != 2 for values in candidates.values()):
        raise ValueError("Expected 12 solutes × 3 species × 2 starting poses")

    selected: list[ase.Atoms] = []
    selected_rows: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []
    for key in sorted(candidates):
        choices = candidates[key]
        converged = [value for value in choices if value[1]["mmff_status"] == 0]
        if not converged:
            raise ValueError(f"No converged MMFF pose for {key}")
        winner = min(
            converged,
            key=lambda value: (
                float(value[1]["mmff_interaction_kcal_mol"]),
                str(value[1]["source_config_id"]),
            ),
        )
        if (
            float(winner[1]["mmff_interaction_kcal_mol"]) >= 0
            or float(winner[1]["minimum_heavy_atom_distance_A"])
            < MIN_HEAVY_ATOM_DISTANCE_A
        ):
            raise ValueError(f"Pre-QM contact criterion failed for {key}")
        selected_id = f"nonwater_contact__{key[0]}__{key[1].lower()}"
        atoms = winner[0].copy()
        atoms.info = {
            "config_id": selected_id,
            "molecule_id": key[0],
            "smiles": str(winner[0].info["smiles"]),
            "perturbant": key[1],
            "n_solute_atoms": int(winner[0].info["n_solute_atoms"]),
            "n_environment_atoms": len(winner[0]) - int(winner[0].info["n_solute_atoms"]),
            "environment_molecule_count": 1,
            "regime": "mmff94s_contact",
            "source_config_id": str(winner[1]["source_config_id"]),
            "geometry_protocol": "fixed_solute_MMFF94s_neighbour_optimization",
        }
        selected.append(atoms)
        row = {"config_id": selected_id, **winner[1]}
        selected_rows.append(row)
        for _, metadata in choices:
            all_rows.append({**metadata, "selected": metadata is winner[1]})

    args.output.mkdir(parents=True, exist_ok=True)
    geometry_file = args.output / "configurations.extxyz"
    write(geometry_file, selected)
    for name, rows in (
        ("configuration_registry.csv", selected_rows),
        ("candidate_registry.csv", all_rows),
    ):
        with (args.output / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    manifest = {
        "status": "geometry-only candidate panel; no new response labels or scores",
        "source_geometry_sha256": sha256(args.source),
        "generator_sha256": sha256(Path(__file__)),
        "geometry_sha256": sha256(geometry_file),
        "ase_version": ase.__version__,
        "rdkit_version": rdkit.__version__,
        "optimization": "MMFF94s, interfragment interactions enabled, all solute atoms fixed, max 1000 iterations",
        "selection": "lower converged MMFF interfragment energy of two archived starting poses per solute × species",
        "selection_requirements": {
            "mmff_interaction_kcal_mol_below": 0,
            "minimum_heavy_atom_distance_A_at_least": MIN_HEAVY_ATOM_DISTANCE_A,
        },
        "n_source_poses": len(frames),
        "n_selected": len(selected),
        "n_solutes": len({str(row["molecule_id"]) for row in selected_rows}),
        "n_neighbour_species": len({str(row["perturbant"]) for row in selected_rows}),
    }
    (args.output / "design_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
