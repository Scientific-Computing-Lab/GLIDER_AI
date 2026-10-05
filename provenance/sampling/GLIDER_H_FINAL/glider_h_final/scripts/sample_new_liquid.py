#!/usr/bin/env python3
"""Generate four independent-block H0 liquid snapshots for new development solutes."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import openmm
import pandas as pd
from openmm import app, unit

ROOT = Path(__file__).resolve().parents[1]
PRIOR = Path("/home/galoren/freesolv_agent_benchmark/hamiltonian_lift/src/h0_systems.py")


def load_h0_module():
    ambertools = "/home/galoren/freesolv_agent_benchmark/.conda_hlift/bin"
    os.environ["PATH"] = ambertools + ":" + os.environ.get("PATH", "")
    spec = importlib.util.spec_from_file_location("glider_h_frozen_h0", PRIOR)
    if spec is None or spec.loader is None:
        raise ImportError(PRIOR)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--selection", type=Path, default=ROOT / "liquid/DEVELOPMENT_SELECTION.csv"
    )
    parser.add_argument("--output", type=Path, default=ROOT / "liquid/raw")
    parser.add_argument("--equilibration-steps", type=int, default=100_000)
    parser.add_argument("--block-steps", type=int, default=50_000)
    args = parser.parse_args()
    h0 = load_h0_module()
    selected = pd.read_csv(args.selection)
    selected = selected[selected.trajectory_source == "new_H0_target_independent_trajectory"]
    args.output.mkdir(parents=True, exist_ok=True)
    for ordinal, row in enumerate(selected.itertuples(), start=1):
        destination = args.output / row.compound_id
        if (destination / "trajectory_manifest.json").exists():
            print(f"already complete: {row.compound_id}", flush=True)
            continue
        destination.mkdir(parents=True, exist_ok=True)
        seed = 2026081500 + ordinal
        started = time.perf_counter()
        bundle = h0.build_h0(row.canonical_isomeric_smiles, seed)
        h0.serialize_h0(
            bundle,
            destination,
            {
                "compound_id": row.compound_id,
                "canonical_isomeric_smiles": row.canonical_isomeric_smiles,
                "seed": seed,
                "experimental_hfe_accessed": False,
            },
        )
        integrator = openmm.LangevinMiddleIntegrator(
            298.15 * unit.kelvin,
            1.0 / unit.picosecond,
            1.0 * unit.femtosecond,
        )
        integrator.setRandomNumberSeed(seed)
        for force in bundle.aqueous_system.getForces():
            if isinstance(force, openmm.MonteCarloBarostat):
                force.setRandomNumberSeed(seed)
        platform, properties = h0.cuda_platform()
        simulation = app.Simulation(
            bundle.topology, bundle.aqueous_system, integrator, platform, properties
        )
        simulation.context.setPositions(bundle.positions)
        simulation.minimizeEnergy(maxIterations=1000)
        simulation.context.setVelocitiesToTemperature(298.15 * unit.kelvin, seed)
        simulation.step(args.equilibration_steps)
        positions, boxes, energies = [], [], []
        for block in range(4):
            simulation.step(args.block_steps)
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            positions.append(
                np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            )
            boxes.append(
                np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
            )
            energies.append(
                float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            )
            print(f"{row.compound_id}: block {block + 1}/4", flush=True)
        trajectory = destination / "trajectory.npz"
        np.savez_compressed(
            trajectory,
            positions_nm=np.asarray(positions),
            box_vectors_nm=np.asarray(boxes),
            potential_energy_kj_mol=np.asarray(energies),
            block_index=np.arange(4, dtype=int),
        )
        manifest = {
            "compound_id": row.compound_id,
            "canonical_isomeric_smiles": row.canonical_isomeric_smiles,
            "h0": "OpenFF Sage 2.2.1/AM1-BCC/TIP3P",
            "seed": seed,
            "equilibration_steps": args.equilibration_steps,
            "block_steps": args.block_steps,
            "block_time_ps": args.block_steps * 0.001,
            "n_blocks": 4,
            "trajectory_sha256": sha256(trajectory),
            "topology_sha256": sha256(destination / "initial_solvated.pdb"),
            "wall_seconds": time.perf_counter() - started,
            "experimental_hfe_accessed": False,
        }
        (destination / "trajectory_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )


if __name__ == "__main__":
    main()
