#!/usr/bin/env python3
"""Frozen MACE-POLAR counterpoise-style electronic-response baseline.

The model is evaluated on the complex, solute, and water cluster independently;
its distributed q/u predictions are subtracted on identical coordinates.  No QM
observable is used by this inference program.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import torch
from ase.io import read
from mace.data import AtomicData, KeySpecification, config_from_atoms
from mace.tools import AtomicNumberTable, torch_geometric

BOHR_TO_ANGSTROM = 0.529177210903


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def graph(atoms, model):
    atoms = atoms.copy()
    atoms.info["charge"] = 0
    atoms.info["spin"] = 1
    return AtomicData.from_config(
        config_from_atoms(atoms, key_specification=KeySpecification()),
        z_table=AtomicNumberTable([int(z) for z in model.atomic_numbers]),
        cutoff=float(model.r_max),
    )


def infer(model, structures, device):
    loader = torch_geometric.dataloader.DataLoader(
        [graph(atoms, model) for atoms in structures],
        batch_size=len(structures),
        shuffle=False,
    )
    batch = next(iter(loader)).to(device).to_dict()
    with torch.no_grad():
        output = model(batch, training=False, compute_force=False)
    values = output["density_coefficients"].detach().cpu().numpy()
    ptr = batch["ptr"].detach().cpu().numpy()
    result = []
    for index in range(len(structures)):
        first, last = int(ptr[index]), int(ptr[index + 1])
        # PolarMACE/e3nn Condon-Shortley order: q, uy, uz, ux.
        q = values[first:last, 0]
        u_e_angstrom = values[first:last][:, [3, 1, 2]]
        result.append((q, u_e_angstrom / BOHR_TO_ANGSTROM))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, action="append", required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--environment-batch-size", type=int, default=16)
    args = parser.parse_args()
    frames = []
    for path in args.configurations:
        frames.extend(read(path, index=":"))
    model = (
        torch.load(args.checkpoint, map_location=args.device, weights_only=False)
        .to(args.device)
        .eval()
    )
    started = time.perf_counter()
    ids, q_values, u_values, offsets = [], [], [], [0]
    for begin in range(0, len(frames), args.environment_batch_size):
        chosen = frames[begin : begin + args.environment_batch_size]
        structures = []
        for atoms in chosen:
            nsolute = int(atoms.info["n_solute_atoms"])
            structures.extend([atoms, atoms[:nsolute], atoms[nsolute:]])
        values = infer(model, structures, args.device)
        for local, atoms in enumerate(chosen):
            full, solute, waters = values[3 * local : 3 * local + 3]
            q = full[0] - np.concatenate([solute[0], waters[0]])
            u = full[1] - np.vstack([solute[1], waters[1]])
            q -= q.mean()
            ids.append(str(atoms.info["config_id"]))
            q_values.append(q.astype(np.float32))
            u_values.append(u.astype(np.float32))
            offsets.append(offsets[-1] + len(atoms))
        print(
            f"MACE-POLAR CP response {min(begin + len(chosen), len(frames))}/{len(frames)}",
            flush=True,
        )
    archive = {
        "config_ids": np.asarray(ids),
        "offsets": np.asarray(offsets, dtype=np.int64),
        "charges_e": np.concatenate(q_values),
        "dipoles_e_bohr": np.concatenate(u_values),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **archive)
    manifest = {
        "baseline": "frozen MACE-POLAR complex-minus-fragments q/u response",
        "checkpoint_sha256": sha256(args.checkpoint),
        "prediction_sha256": sha256(args.output),
        "n_configurations": len(frames),
        "runtime_seconds_excluding_checkpoint_load": time.perf_counter() - started,
        "reference_observables_accessed": False,
        "prospective_labels_accessed": False,
        "experimental_hydration_targets_accessed": False,
    }
    args.output.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
