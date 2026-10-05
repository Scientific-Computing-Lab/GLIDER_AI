#!/usr/bin/env python3
"""Cache frozen MACE-POLAR product-basis features for response learning."""

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


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def graph(atoms, model):
    candidate = atoms.copy()
    candidate.info["charge"] = 0
    candidate.info["spin"] = 1
    table = AtomicNumberTable([int(z) for z in model.atomic_numbers])
    return AtomicData.from_config(
        config_from_atoms(candidate, key_specification=KeySpecification()),
        z_table=table,
        cutoff=float(model.r_max),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, action="append", required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--batch-size", type=int, default=8)
    args = parser.parse_args()
    frames = []
    for path in args.configurations:
        frames.extend(read(path, index=":"))
    ids = [str(a.info["config_id"]) for a in frames]
    if len(ids) != len(set(ids)):
        raise RuntimeError("Duplicate configuration identifiers")
    model = (
        torch.load(args.checkpoint, map_location=args.device, weights_only=False)
        .to(args.device)
        .eval()
    )
    captured = [None] * len(model.products)
    handles = []
    for layer, product in enumerate(model.products):

        def hook(_module, _inputs, output, index=layer):
            captured[index] = output

        handles.append(product.register_forward_hook(hook))
    output = {"config_ids": np.asarray(ids), "offsets": [0], "atomic_numbers": []}
    for layer in range(len(model.products)):
        output[f"product_{layer}"] = []
    started = time.perf_counter()
    try:
        for begin in range(0, len(frames), args.batch_size):
            chosen = frames[begin : begin + args.batch_size]
            loader = torch_geometric.dataloader.DataLoader(
                [graph(atoms, model) for atoms in chosen],
                batch_size=len(chosen),
                shuffle=False,
            )
            batch = next(iter(loader)).to(args.device).to_dict()
            captured[:] = [None] * len(captured)
            with torch.no_grad():
                model(batch, training=False, compute_force=False)
            ptr = batch["ptr"].detach().cpu().numpy()
            for local, atoms in enumerate(chosen):
                first, last = int(ptr[local]), int(ptr[local + 1])
                output["atomic_numbers"].append(
                    np.asarray(atoms.numbers, dtype=np.int16)
                )
                for layer, values in enumerate(captured):
                    output[f"product_{layer}"].append(
                        values[first:last].detach().cpu().numpy().astype(np.float32)
                    )
                output["offsets"].append(output["offsets"][-1] + len(atoms))
            print(
                f"features {min(begin + len(chosen), len(frames))}/{len(frames)}",
                flush=True,
            )
    finally:
        for handle in handles:
            handle.remove()
    archive = {
        "config_ids": output["config_ids"],
        "offsets": np.asarray(output["offsets"], dtype=np.int64),
        "atomic_numbers": np.concatenate(output["atomic_numbers"]).astype(np.int16),
    }
    for layer in range(len(model.products)):
        archive[f"product_{layer}"] = np.concatenate(output[f"product_{layer}"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **archive)
    manifest = {
        "encoder": "MACE-POLAR frozen product-basis features",
        "checkpoint_sha256": sha256(args.checkpoint),
        "feature_sha256": sha256(args.output),
        "n_configurations": len(frames),
        "n_atoms": int(archive["offsets"][-1]),
        "runtime_seconds_excluding_checkpoint_load": time.perf_counter() - started,
        "product_shapes": {
            key: list(value.shape)
            for key, value in archive.items()
            if key.startswith("product_")
        },
        "prospective_labels_accessed": False,
        "experimental_hydration_targets_accessed": False,
    }
    path = args.output.with_suffix(".manifest.json")
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
