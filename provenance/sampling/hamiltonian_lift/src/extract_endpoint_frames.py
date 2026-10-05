#!/usr/bin/env python3
"""Extract fully coupled endpoint frames from frozen replica-exchange files."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from openmm import unit
from openmmtools import multistate


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reporter", type=Path, required=True)
    parser.add_argument("--phase", choices=["aqueous", "gas"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint-interval", type=int, default=10)
    parser.add_argument("--thermodynamic-state", type=int, default=0)
    args = parser.parse_args()

    reporter = multistate.MultiStateReporter(str(args.reporter), open_mode="r")
    state_trace = reporter.read_replica_thermodynamic_states()
    final_iteration = state_trace.shape[0] - 1
    iterations = list(range(args.checkpoint_interval, final_iteration + 1, args.checkpoint_interval))
    positions = []
    boxes = []
    replicas = []
    for iteration in iterations:
        state_indices = reporter.read_replica_thermodynamic_states(iteration)
        match = np.flatnonzero(state_indices == args.thermodynamic_state)
        if len(match) != 1:
            raise RuntimeError(f"Iteration {iteration}: expected one replica at coupled state 0, found {match}")
        sampler_states = reporter.read_sampler_states(iteration, analysis_particles_only=False)
        if sampler_states is None:
            raise RuntimeError(f"No full checkpoint at iteration {iteration}")
        replica = int(match[0])
        sampler_state = sampler_states[replica]
        positions.append(np.asarray(sampler_state.positions.value_in_unit(unit.nanometer), dtype=np.float64))
        if sampler_state.box_vectors is None:
            boxes.append(np.full((3, 3), np.nan, dtype=np.float64))
        else:
            boxes.append(np.asarray(sampler_state.box_vectors.value_in_unit(unit.nanometer), dtype=np.float64))
        replicas.append(replica)
    reporter.close()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        positions_nm=np.asarray(positions),
        box_vectors_nm=np.asarray(boxes),
        iteration=np.asarray(iterations, dtype=int),
        replica=np.asarray(replicas, dtype=int),
    )
    manifest = {
        "source_reporter": str(args.reporter),
        "source_reporter_sha256": sha256(args.reporter),
        "phase": args.phase,
        "thermodynamic_state": args.thermodynamic_state,
        "lambda_electrostatics": 1.0,
        "lambda_sterics": 1.0,
        "n_frames": len(iterations),
        "iterations": iterations,
        "output": str(args.output),
        "output_sha256": sha256(args.output),
        "experimental_target_accessed": False,
    }
    args.output.with_suffix(".json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
