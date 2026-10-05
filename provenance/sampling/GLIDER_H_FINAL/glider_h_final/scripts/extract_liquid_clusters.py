#!/usr/bin/env python3
"""Extract deterministic nearest-water 3/5/8-water clusters from H0 liquid frames."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
PRIOR = Path("/home/galoren/freesolv_agent_benchmark/hamiltonian_lift")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    fractional = np.linalg.solve(box.T, np.asarray(delta).T).T
    fractional -= np.rint(fractional)
    return fractional @ box


def local_cluster(
    positions_a: np.ndarray,
    box_a: np.ndarray,
    symbols: list[str],
    n_solute: int,
    n_waters: int,
) -> Atoms:
    if (len(symbols) - n_solute) % 3:
        raise RuntimeError("TIP3P topology is not solute followed by OHH waters")
    solute = positions_a[:n_solute].copy()
    waters = positions_a[n_solute:].reshape(-1, 3, 3)
    ranked = []
    for index, water in enumerate(waters):
        displacements = minimum_image(water[0] - solute, box_a)
        nearest_atom = int(np.argmin(np.linalg.norm(displacements, axis=1)))
        oxygen = solute[nearest_atom] + displacements[nearest_atom]
        hydrogen_1 = oxygen + minimum_image(water[1] - water[0], box_a)
        hydrogen_2 = oxygen + minimum_image(water[2] - water[0], box_a)
        distance = float(np.linalg.norm(oxygen - solute[nearest_atom]))
        ranked.append((distance, index, np.vstack([oxygen, hydrogen_1, hydrogen_2])))
    chosen = sorted(ranked, key=lambda value: (value[0], value[1]))[:n_waters]
    coordinates = np.vstack([solute, *[value[2] for value in chosen]])
    cluster_symbols = symbols[:n_solute] + [item for _ in chosen for item in ("O", "H", "H")]
    atoms = Atoms(cluster_symbols, positions=coordinates)
    atoms.info["nearest_water_distances_A"] = ";".join(f"{value[0]:.6f}" for value in chosen)
    return atoms


def main() -> None:
    selected = pd.read_csv(ROOT / "liquid/DEVELOPMENT_SELECTION.csv")
    output = ROOT / "liquid/clusters"
    output.mkdir(parents=True, exist_ok=True)
    by_shell: dict[int, list[Atoms]] = {3: [], 5: [], 8: []}
    rows = []
    for row in selected.itertuples():
        if row.trajectory_source == "prior_H0_target_independent_trajectory":
            source = PRIOR / "data/endpoints" / row.compound_id / "aqueous.npz"
            topology = PRIOR / "data/h0" / row.compound_id / "initial_solvated.pdb"
            manifest = json.loads(
                (PRIOR / "data/h0" / row.compound_id / "system_manifest.json").read_text()
            )
            values = np.load(source)
            indices = np.linspace(0, len(values["positions_nm"]) - 1, 4, dtype=int)
        else:
            source = ROOT / "liquid/raw" / row.compound_id / "trajectory.npz"
            topology = ROOT / "liquid/raw" / row.compound_id / "initial_solvated.pdb"
            manifest = json.loads(
                (ROOT / "liquid/raw" / row.compound_id / "system_manifest.json").read_text()
            )
            values = np.load(source)
            indices = np.arange(4, dtype=int)
        topology_atoms = read(topology)
        symbols = topology_atoms.get_chemical_symbols()
        n_solute = int(manifest["n_solute_atoms"])
        for local_index, source_index in enumerate(indices):
            for shell in (3, 5, 8):
                if shell == 5 and local_index not in (1, 3):
                    continue
                if shell == 8 and local_index != 2:
                    continue
                atoms = local_cluster(
                    np.asarray(values["positions_nm"][source_index]) * 10.0,
                    np.asarray(values["box_vectors_nm"][source_index]) * 10.0,
                    symbols,
                    n_solute,
                    shell,
                )
                config_id = f"glider_h__{row.compound_id}__liquid_{local_index:02d}__{shell}w"
                atoms.info.update(
                    {
                        "config_id": config_id,
                        "molecule_id": row.compound_id,
                        "regime": "liquid",
                        "shell_waters": shell,
                        "n_solute_atoms": n_solute,
                        "source_frame_index": int(source_index),
                    }
                )
                by_shell[shell].append(atoms)
                rows.append(
                    {
                        "config_id": config_id,
                        "molecule_id": row.compound_id,
                        "shell_waters": shell,
                        "source_frame_index": int(source_index),
                        "source_trajectory": str(source),
                        "source_trajectory_sha256": sha256(source),
                        "n_solute_atoms": n_solute,
                        "n_atoms": len(atoms),
                        "nearest_water_distances_A": atoms.info[
                            "nearest_water_distances_A"
                        ],
                        "experimental_hfe_accessed": False,
                    }
                )
    for shell, frames in by_shell.items():
        write(output / f"liquid_{shell}water.extxyz", frames, format="extxyz")
    registry = ROOT / "liquid/LIQUID_CLUSTER_REGISTRY.csv"
    pd.DataFrame(rows).sort_values(["shell_waters", "config_id"]).to_csv(registry, index=False)
    manifest = {
        "selection_sha256": sha256(ROOT / "liquid/DEVELOPMENT_SELECTION.csv"),
        "registry_sha256": sha256(registry),
        "n_3water": len(by_shell[3]),
        "n_5water": len(by_shell[5]),
        "n_8water": len(by_shell[8]),
        "geometry_hashes": {
            f"liquid_{shell}water.extxyz": sha256(output / f"liquid_{shell}water.extxyz")
            for shell in (3, 5, 8)
        },
        "experimental_hfe_accessed": False,
    }
    (ROOT / "liquid/LIQUID_CLUSTER_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
