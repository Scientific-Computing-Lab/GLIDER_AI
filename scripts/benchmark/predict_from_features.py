#!/usr/bin/env python3
"""Run the frozen GLIDER head from cached MACE features and response priors.

The byte-preserved historical implementation remains in ``evidence/frozen_source``.
This release-facing implementation is accepted only because its outputs are
checked against every stored prospective prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch
from ase.io import read

from glider.inference import BOHR_TO_ANGSTROM, DEBYE_PER_E_BOHR, molecular_surface

ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def response_head_class():
    source = ROOT / "evidence/frozen_source/prospective_1/response_learning.py"
    spec = importlib.util.spec_from_file_location("glider_frozen_response_learning", source)
    if spec is None or spec.loader is None:
        raise ImportError(source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.ResponseHead, module.esp_design, module.geometry_basis


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, action="append", required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--base-predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    frames = [atoms for path in args.configurations for atoms in read(path, index=":")]
    frame_map = {str(atoms.info["config_id"]): atoms for atoms in frames}
    features = np.load(args.features)
    feature_ids = [str(value) for value in features["config_ids"]]
    base = np.load(args.base_predictions)
    if not np.array_equal(features["config_ids"], base["config_ids"]):
        raise RuntimeError("feature and prior configuration identifiers differ")

    ResponseHead, esp_design, geometry_basis = response_head_class()
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    models = [
        ResponseHead.from_state_dict(state, checkpoint["head"]).eval()
        for state in checkpoint["states"]
    ]
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    started = time.perf_counter()
    for index, config_id in enumerate(feature_ids):
        atoms = frame_map[config_id]
        first, last = int(features["offsets"][index]), int(features["offsets"][index + 1])
        base_first, base_last = int(base["offsets"][index]), int(base["offsets"][index + 1])
        if last - first != len(atoms) or (base_first, base_last) != (first, last):
            raise RuntimeError(f"atom offsets differ for {config_id}")
        product = np.asarray(features["product_1"][first:last])
        channels = product.shape[1] // 4
        scalar = product[:, :channels]
        vector = product[:, channels:].reshape(len(atoms), channels, 3)[:, :, [2, 0, 1]]
        sample = SimpleNamespace(
            scalars=torch.tensor(scalar, dtype=torch.float32),
            vectors=torch.tensor(vector, dtype=torch.float32),
            geometry_vectors=torch.tensor(
                geometry_basis(atoms.positions, atoms.numbers), dtype=torch.float32
            ),
            positions=torch.tensor(np.asarray(atoms.positions), dtype=torch.float32),
            base_q=torch.tensor(base["charges_e"][first:last], dtype=torch.float32),
            base_u=torch.tensor(base["dipoles_e_bohr"][first:last], dtype=torch.float32),
        )
        with torch.no_grad():
            predictions = [model(sample) for model in models]
        charges = torch.stack([value[0] for value in predictions]).mean(0).numpy()
        site_dipoles = torch.stack([value[1] for value in predictions]).mean(0).numpy()
        charges -= charges.mean()
        points = molecular_surface(atoms.get_chemical_symbols(), atoms.positions)
        response_esp = esp_design(np.asarray(atoms.positions), points) @ np.r_[
            charges, site_dipoles.reshape(-1)
        ]
        response_dipole = (
            np.sum(charges[:, None] * (atoms.positions / BOHR_TO_ANGSTROM), axis=0)
            + np.sum(site_dipoles, axis=0)
        ) * DEBYE_PER_E_BOHR
        filename = f"{hashlib.sha256(config_id.encode()).hexdigest()[:20]}.npz"
        path = args.output / filename
        np.savez_compressed(
            path,
            points_angstrom=points,
            predicted_esp_hartree_per_e=response_esp,
            predicted_dipole_debye=response_dipole,
            predicted_charges_e=charges,
            predicted_dipoles_e_bohr=site_dipoles,
        )
        rows.append(
            {
                "config_id": config_id,
                "molecule_id": str(atoms.info["molecule_id"]),
                "regime": str(atoms.info["regime"]),
                "prediction_file": filename,
                "prediction_sha256": sha256(path),
                "net_charge_e": float(charges.sum()),
            }
        )

    registry = args.output / "prediction_registry.csv"
    pd.DataFrame(rows).sort_values("config_id").to_csv(registry, index=False)
    manifest = {
        "implementation": "clean cached-feature predictor using byte-preserved ResponseHead",
        "checkpoint_sha256": sha256(args.checkpoint),
        "feature_sha256": sha256(args.features),
        "base_predictions_sha256": sha256(args.base_predictions),
        "prediction_registry_sha256": sha256(registry),
        "n_configurations": len(rows),
        "runtime_seconds_total_excluding_encoder": time.perf_counter() - started,
        "reference_labels_accessed": False,
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
