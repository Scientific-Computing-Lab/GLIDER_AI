#!/usr/bin/env python3
"""Materialize the four preregistered geometries for 20 selected identities."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import write

from geometry import build_molecule_environments


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    selected = pd.read_csv(args.selected)
    if len(selected) != 20 or selected.freesolv_id.nunique() != 20:
        raise RuntimeError("Panel III requires exactly 20 selected identities")
    protocol_hash = sha256(args.protocol)
    frames, rows = [], []
    for row in selected.sort_values("freesolv_id").itertuples():
        molecule_frames = build_molecule_environments(str(row.freesolv_id), str(row.canonical_smiles), f"stratum_{row.stratum}", protocol_hash)
        for atoms in molecule_frames:
            nsolute = int(atoms.info["n_solute_atoms"])
            frames.append(atoms)
            rows.append({
                **dict(atoms.info),
                "source_dataset": "FreeSolv_v0.52",
                "source_commit": "6c7d19b4b565537365ffd22006aa2cd4643200c6",
                "n_atoms": len(atoms),
                "minimum_solute_water_distance_A": float(np.linalg.norm(atoms.positions[:nsolute, None] - atoms.positions[None, nsolute:], axis=2).min()),
            })
    frames.sort(key=lambda atoms: str(atoms.info["config_id"]))
    registry = pd.DataFrame(rows).sort_values("config_id").reset_index(drop=True)
    if len(frames) != 80 or registry.molecule_id.nunique() != 20 or registry.groupby("regime").size().to_dict() != {"compressed": 20, "cooperative": 20, "equilibrium": 20, "orientation": 20}:
        raise RuntimeError("Panel III geometry count mismatch")
    if registry.config_id.tolist() != [str(atoms.info["config_id"]) for atoms in frames]:
        raise RuntimeError("Frame/registry order mismatch")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    xyz = args.output_dir / "configurations.extxyz"
    csv_path = args.output_dir / "configuration_registry.csv"
    write(xyz, frames, format="extxyz")
    registry.to_csv(csv_path, index=False)
    manifest = {
        "algorithm": "unchanged deterministic ETKDG/MMFF plus frozen structured three-water regime generator",
        "protocol_sha256": protocol_hash,
        "selected_molecules_sha256": sha256(args.selected),
        "n_molecules": 20,
        "n_configurations": 80,
        "counts_by_regime": {key: int(value) for key, value in registry.groupby("regime").size().items()},
        "configuration_sha256": sha256(xyz),
        "registry_sha256": sha256(csv_path),
        "reference_response_labels_exist": False,
        "experimental_hydration_targets_accessed": False,
    }
    (args.output_dir / "GEOMETRY_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
