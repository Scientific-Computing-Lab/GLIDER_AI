#!/usr/bin/env python3
"""Materialize frozen q/u response archives for the confirmatory panel."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read

ROOT = Path(__file__).resolve().parents[4]
FROZEN = ROOT / "evidence/frozen_source/prospective_1"
sys.path.insert(0, str(FROZEN))
from response_learning import esp_design  # noqa: E402

from glider.inference import (  # noqa: E402
    BOHR_TO_ANGSTROM,
    DEBYE_PER_E_BOHR,
    molecular_surface,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    frames = read(args.configurations, index=":")
    frame_map = {str(atoms.info["config_id"]): atoms for atoms in frames}
    values = np.load(args.predictions, allow_pickle=False)
    ids = [str(value) for value in values["config_ids"]]
    offsets = values["offsets"]
    if set(ids) != set(frame_map):
        raise RuntimeError("configuration and response archive sets differ")

    args.output.mkdir(parents=True, exist_ok=False)
    rows: list[dict[str, object]] = []
    for index, config_id in enumerate(ids):
        atoms = frame_map[config_id]
        first, last = int(offsets[index]), int(offsets[index + 1])
        q = np.asarray(values["charges_e"][first:last], dtype=float)
        u = np.asarray(values["dipoles_e_bohr"][first:last], dtype=float)
        q -= q.mean()
        points = molecular_surface(atoms.get_chemical_symbols(), atoms.positions)
        predicted_esp = esp_design(atoms.positions, points) @ np.r_[q, u.reshape(-1)]
        predicted_dipole = (
            np.sum(q[:, None] * (atoms.positions / BOHR_TO_ANGSTROM), axis=0)
            + np.sum(u, axis=0)
        ) * DEBYE_PER_E_BOHR
        tag = hashlib.sha256(config_id.encode()).hexdigest()[:20]
        output_file = args.output / f"{tag}.npz"
        np.savez_compressed(
            output_file,
            points_angstrom=points,
            predicted_esp_hartree_per_e=predicted_esp,
            predicted_dipole_debye=predicted_dipole,
            predicted_charges_e=q,
            predicted_dipoles_e_bohr=u,
        )
        rows.append(
            {
                "config_id": config_id,
                "molecule_id": str(atoms.info["molecule_id"]),
                "family": str(atoms.info["family"]),
                "regime": str(atoms.info["regime"]),
                "n_waters": int(atoms.info["n_waters"]),
                "source_frame_index": int(atoms.info["source_frame_index"]),
                "prediction_file": output_file.name,
                "prediction_sha256": sha256(output_file),
                "net_charge_e": float(q.sum()),
            }
        )

    registry = args.output / "prediction_registry.csv"
    pd.DataFrame(rows).sort_values("config_id").to_csv(registry, index=False)
    manifest = {
        "method": args.name,
        "source_prediction_sha256": sha256(args.predictions),
        "configuration_sha256": sha256(args.configurations),
        "registry_sha256": sha256(registry),
        "n_configurations": len(rows),
        "probe_definition": "published deterministic two-shell whole-cluster surface",
        "reference_labels_accessed": False,
        "prospective_claim": False,
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
