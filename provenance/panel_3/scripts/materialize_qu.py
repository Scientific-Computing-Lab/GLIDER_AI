#!/usr/bin/env python3
"""Materialize a frozen distributed q/u archive on the fixed probe surface."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read

from glider.inference import dipole_from_sites, esp_from_sites, molecular_surface


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frames = read(args.configurations, index=":")
    frame_map = {str(atoms.info["config_id"]): atoms for atoms in frames}
    values = np.load(args.archive)
    ids = [str(value) for value in values["config_ids"]]
    if set(ids) != set(frame_map):
        raise RuntimeError("Archive and geometry sets differ")
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    for index, config_id in enumerate(ids):
        atoms = frame_map[config_id]
        first, last = int(values["offsets"][index]), int(values["offsets"][index + 1])
        q = np.asarray(values["charges_e"][first:last], dtype=float)
        u = np.asarray(values["dipoles_e_bohr"][first:last], dtype=float)
        points = molecular_surface(atoms.get_chemical_symbols(), atoms.positions)
        filename = f"{hashlib.sha256(config_id.encode()).hexdigest()[:20]}.npz"
        path = args.output / filename
        np.savez_compressed(path, points_angstrom=points, predicted_esp_hartree_per_e=esp_from_sites(atoms.positions, points, q, u), predicted_dipole_debye=dipole_from_sites(atoms.positions, q, u), predicted_charges_e=q, predicted_dipoles_e_bohr=u)
        rows.append({"config_id": config_id, "molecule_id": str(atoms.info["molecule_id"]), "regime": str(atoms.info["regime"]), "prediction_file": filename, "prediction_sha256": sha256(path), "net_charge_e": float(q.sum())})
    registry = args.output / "prediction_registry.csv"
    pd.DataFrame(rows).sort_values("config_id").to_csv(registry, index=False)
    (args.output / "manifest.json").write_text(json.dumps({"method": args.method, "source_archive_sha256": sha256(args.archive), "prediction_registry_sha256": sha256(registry), "n_configurations": len(rows), "reference_labels_accessed": False}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
