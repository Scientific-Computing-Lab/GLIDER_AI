#!/usr/bin/env python3
"""Materialize frozen MACE-POLAR q/u response sites on GLIDER's probe surface."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

from glider.inference import dipole_from_sites, esp_from_sites, molecular_surface


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configurations", type=Path, required=True)
    parser.add_argument("--raw-sites", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--glider-predictions", type=Path, required=True)
    args = parser.parse_args()
    frames = {str(a.info["config_id"]): a for a in read(args.configurations, index=":")}
    with np.load(args.raw_sites, allow_pickle=False) as archive:
        ids = [str(value) for value in archive["config_ids"]]
        offsets = np.asarray(archive["offsets"], dtype=int)
        charges = np.asarray(archive["charges_e"], dtype=float)
        dipoles = np.asarray(archive["dipoles_e_bohr"], dtype=float)
    if len(ids) != 36 or set(ids) != set(frames):
        raise ValueError("MACE-L raw sites do not cover exactly the 36 frozen geometries")
    with (args.glider_predictions / "prediction_registry.csv").open(newline="") as stream:
        glider = {row["config_id"]: row for row in csv.DictReader(stream)}
    if set(glider) != set(ids):
        raise ValueError("GLIDER probe registry differs from the MACE-L cases")
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    for index, config_id in enumerate(ids):
        atoms = frames[config_id]
        first, last = int(offsets[index]), int(offsets[index + 1])
        if last - first != len(atoms):
            raise ValueError(f"Atom offset mismatch for {config_id}")
        q = charges[first:last].copy()
        u = dipoles[first:last].copy()
        q -= q.mean()
        points = molecular_surface(atoms.get_chemical_symbols(), atoms.positions)
        with np.load(
            args.glider_predictions / glider[config_id]["prediction_file"],
            allow_pickle=False,
        ) as reference_probes:
            if not np.array_equal(points, reference_probes["points_angstrom"]):
                raise ValueError(f"Probe grid mismatch for {config_id}")
        response_esp = esp_from_sites(atoms.positions, points, q, u)
        response_dipole = dipole_from_sites(atoms.positions, q, u)
        tag = hashlib.sha256(config_id.encode()).hexdigest()[:20]
        path = args.output / f"{tag}.npz"
        np.savez_compressed(
            path,
            points_angstrom=points,
            predicted_esp_hartree_per_e=response_esp,
            predicted_dipole_debye=response_dipole,
            predicted_charges_e=q,
            predicted_dipoles_e_bohr=u,
        )
        rows.append({
            "config_id": config_id,
            "molecule_id": str(atoms.info["molecule_id"]),
            "perturbant": str(atoms.info["perturbant"]),
            "prediction_file": path.name,
            "prediction_sha256": sha256(path),
        })
    registry = args.output / "prediction_registry.csv"
    with registry.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda item: item["config_id"]))
    (args.output / "manifest.json").write_text(json.dumps({
        "method": "MACE-POLAR-1-L complex-minus-fragments q/u, no response fitting",
        "configurations_sha256": sha256(args.configurations),
        "raw_sites_sha256": sha256(args.raw_sites),
        "prediction_registry_sha256": sha256(registry),
        "n_configurations": len(rows),
        "reference_labels_accessed": False,
    }, indent=2, sort_keys=True) + "\n")
    print(f"MACE-L: {len(rows)} matching probe arrays, registry {sha256(registry)}")


if __name__ == "__main__":
    main()
