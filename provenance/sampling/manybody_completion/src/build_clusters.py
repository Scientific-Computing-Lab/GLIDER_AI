#!/usr/bin/env python3
"""Create deterministic 1-/2-/3-/4-water panels with physical motifs."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import write


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "complete_interaction_hamiltonian/src"))
from build_configurations import (  # noqa: E402
    build_cluster,
    conformers,
    direction_outward,
    make_atoms,
    orient_water,
    place_at_minimum_distance,
    random_rotation,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_water(candidate, solute, waters, min_solute=1.30, min_inter=1.15):
    if np.linalg.norm(solute[:, None] - candidate[None, :], axis=2).min() < min_solute:
        return False
    return all(
        np.linalg.norm(candidate[:, None] - water[None, :], axis=2).min() >= min_inter
        for water in waters
    )


def cooperative_chain(solute, anchors, nwater, rng, compressed=False):
    """First-shell water plus an oriented water-water H-bond chain."""
    anchor = anchors[int(rng.integers(len(anchors)))]
    radial = direction_outward(solute, anchor, rng)
    first_distance = float(rng.uniform(1.55, 1.75) if compressed else rng.uniform(1.85, 2.35))
    first = place_at_minimum_distance(solute, anchor, radial, first_distance, "water_donor", rng)
    waters = [first]
    for _ in range(1, nwater):
        accepted = None
        for _attempt in range(300):
            direction = waters[-1][0] - solute.mean(0) + 0.5 * rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            oo = float(rng.uniform(2.65, 2.95))
            candidate = orient_water(direction, "water_donor", rng) + waters[-1][0] + oo * direction
            if valid_water(candidate, solute, waters):
                accepted = candidate
                break
        if accepted is None:
            # Deterministic fallback remains a valid cluster but is marked by
            # the caller's motif rather than silently dropping the geometry.
            independent = build_cluster(solute, anchors, 1, rng, "equilibrium")[0]
            if not valid_water(independent, solute, waters):
                raise RuntimeError("Failed cooperative-chain placement")
            accepted = independent
        waters.append(accepted)
    return waters


def outer_shell(solute, anchors, nwater, rng):
    """Two first-shell waters and a connected but solute-distant outer shell."""
    first_shell = build_cluster(solute, anchors, max(1, nwater - 1), rng, "equilibrium")
    parent = first_shell[int(rng.integers(len(first_shell)))]
    accepted = None
    for _ in range(400):
        direction = parent[0] - solute.mean(0) + 0.35 * rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        candidate = orient_water(direction, "water_donor", rng) + parent[0] + 2.8 * direction
        min_solute = np.linalg.norm(solute[:, None] - candidate[None, :], axis=2).min()
        if min_solute >= 3.2 and valid_water(candidate, solute, first_shell):
            accepted = candidate
            break
    if accepted is None:
        raise RuntimeError("Failed outer-shell placement")
    return first_shell + [accepted]


def reorient(waters, rng):
    result = []
    for water in waters:
        result.append((water - water[0]) @ random_rotation(rng).T + water[0])
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", choices=("development", "prospective"), required=True)
    args = parser.parse_args()
    seed = 20260821 if args.panel == "development" else 20260829
    molecules = pd.read_csv(ROOT / f"config/molecules_{args.panel}.csv")
    output = ROOT / f"data/configurations_{args.panel}"
    output.mkdir(parents=True, exist_ok=True)
    frames, rows = [], []
    for molecule_index, row in molecules.iterrows():
        rng = np.random.default_rng(seed + molecule_index * 1009)
        molecule, conformer_coordinates = conformers(row.smiles, seed + molecule_index)
        solute = conformer_coordinates[0]
        heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
        hetero = [
            atom.GetIdx()
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() in (7, 8, 9, 16, 17)
        ]
        anchors = hetero + [index for index in heavy if index not in hetero]
        first_anchor = anchors[0]
        radial = direction_outward(solute, first_anchor, rng)
        one = [place_at_minimum_distance(solute, first_anchor, radial, 2.15, "water_donor", rng)]
        specifications = [
            ("control_1water", "zero", one, "pair_zero"),
            ("cluster_2water", "independent", build_cluster(solute, anchors, 2, rng, "equilibrium"), "equilibrium"),
            ("cluster_2water", "cooperative", cooperative_chain(solute, anchors, 2, rng), "cooperative"),
            ("cluster_2water", "compressed", cooperative_chain(solute, anchors, 2, rng, compressed=True), "compressed"),
            ("cluster_3water", "independent", build_cluster(solute, anchors, 3, rng, "equilibrium"), "equilibrium"),
            ("cluster_3water", "cooperative", cooperative_chain(solute, anchors, 3, rng), "cooperative"),
            ("cluster_3water", "outer", outer_shell(solute, anchors, 3, rng), "outer_shell"),
            ("cluster_4water", "compact", build_cluster(solute, anchors, 4, rng, "equilibrium"), "equilibrium"),
            ("cluster_4water", "cooperative", cooperative_chain(solute, anchors, 4, rng), "cooperative"),
            ("cluster_4water", "orientation", reorient(build_cluster(solute, anchors, 4, rng, "equilibrium"), rng), "orientation"),
        ]
        for category, suffix, waters, regime in specifications:
            config_id = f"{row.split}__{row.molecule_id}__{category}__{suffix}"
            metadata = {
                "config_id": config_id,
                "molecule_id": row.molecule_id,
                "smiles": row.smiles,
                "split": row.split,
                "chemical_role": row.chemical_role,
                "category": category,
                "motif": suffix,
                "regime": regime,
                "conformer_index": 0,
            }
            atoms = make_atoms(molecule, solute, waters, metadata)
            min_sw = float(
                np.linalg.norm(
                    atoms.positions[: len(solute), None]
                    - atoms.positions[None, len(solute) :],
                    axis=2,
                ).min()
            )
            oxygen_indices = [len(solute) + 3 * i for i in range(len(waters))]
            min_oo = (
                min(
                    np.linalg.norm(atoms.positions[i] - atoms.positions[j])
                    for a, i in enumerate(oxygen_indices)
                    for j in oxygen_indices[a + 1 :]
                )
                if len(oxygen_indices) > 1
                else float("nan")
            )
            frames.append(atoms)
            rows.append(
                {
                    **metadata,
                    "n_atoms": len(atoms),
                    "n_solute_atoms": len(solute),
                    "n_waters": len(waters),
                    "minimum_solute_water_distance_A": min_sw,
                    "minimum_water_oxygen_distance_A": min_oo,
                }
            )
    order = np.argsort([atoms.info["config_id"] for atoms in frames])
    frames = [frames[index] for index in order]
    registry = pd.DataFrame(rows).sort_values("config_id").reset_index(drop=True)
    if registry.config_id.tolist() != [atoms.info["config_id"] for atoms in frames]:
        raise AssertionError("Registry/frame order mismatch")
    xyz = output / "configurations.extxyz"
    csv_path = output / "configuration_registry.csv"
    write(xyz, frames, format="extxyz")
    registry.to_csv(csv_path, index=False)
    manifest = {
        "panel": args.panel,
        "seed": seed,
        "n_molecules": int(registry.molecule_id.nunique()),
        "n_configurations": len(registry),
        "counts_by_body_order": registry.groupby("n_waters").size().astype(int).to_dict(),
        "counts_by_regime": registry.groupby("regime").size().astype(int).to_dict(),
        "molecules": sorted(registry.molecule_id.unique()),
        "configuration_sha256": sha256(xyz),
        "registry_sha256": sha256(csv_path),
        "experimental_hydration_targets_accessed": False,
    }
    (output / "configuration_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
