#!/usr/bin/env python3
"""Run and analyze a reproducible H0 absolute hydration free-energy cycle."""

from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
import openmm
from openmm import app, unit
from openmmtools import alchemy, cache, mcmc, multistate, states

from h0_systems import PRESSURE, TEMPERATURE, TIMESTEP, cuda_platform


KB_KCAL_PER_MOL_K = 0.00198720425864083


def schedule(phase: str) -> list[tuple[float, float]]:
    electrostatics = [1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0.0]
    if phase == "gas":
        return [(x, 1.0) for x in electrostatics]
    sterics = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    return [(x, 1.0) for x in electrostatics] + [(0.0, x) for x in sterics]


def create_alchemical_system(reference: openmm.System, n_solute_atoms: int) -> openmm.System:
    region = alchemy.AlchemicalRegion(
        alchemical_atoms=list(range(n_solute_atoms)),
        annihilate_electrostatics=True,
        annihilate_sterics=False,
    )
    factory = alchemy.AbsoluteAlchemicalFactory(
        alchemical_pme_treatment="exact",
        disable_alchemical_dispersion_correction=False,
        split_alchemical_forces=True,
    )
    return factory.create_alchemical_system(reference, region)


def create_states(system: openmm.System, phase: str) -> list[states.CompoundThermodynamicState]:
    pressure = PRESSURE if phase == "aqueous" else None
    thermo = states.ThermodynamicState(system=system, temperature=TEMPERATURE, pressure=pressure)
    result = []
    for lambda_electrostatics, lambda_sterics in schedule(phase):
        alchemical_state = alchemy.AlchemicalState.from_system(system)
        alchemical_state.lambda_electrostatics = lambda_electrostatics
        alchemical_state.lambda_sterics = lambda_sterics
        result.append(states.CompoundThermodynamicState(copy.deepcopy(thermo), [alchemical_state]))
    return result


def analyze(reporter: multistate.MultiStateReporter, phase: str) -> dict:
    analyzer = multistate.MultiStateSamplerAnalyzer(reporter, n_equilibration_iterations=0)
    delta_f, d_delta_f = analyzer.get_free_energy()
    kT = KB_KCAL_PER_MOL_K * TEMPERATURE.value_in_unit(unit.kelvin)
    overlap = analyzer.mbar.compute_overlap()
    matrix = np.asarray(overlap["matrix"], dtype=float)
    adjacent = [float(matrix[i, i + 1] + matrix[i + 1, i]) for i in range(len(matrix) - 1)]
    return {
        "phase": phase,
        "n_states": int(len(matrix)),
        "decoupling_free_energy_kcal_mol": float(delta_f[0, -1] * kT),
        "decoupling_uncertainty_kcal_mol": float(d_delta_f[0, -1] * kT),
        "minimum_bidirectional_adjacent_overlap": float(min(adjacent)),
        "adjacent_bidirectional_overlap": adjacent,
        "overlap_scalar": float(overlap["scalar"]),
        "overlap_matrix": matrix.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--system-dir", type=Path, required=True)
    parser.add_argument("--phase", choices=["aqueous", "gas"], required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--equil-iterations", type=int, default=100)
    parser.add_argument("--iterations", type=int, default=500)
    parser.add_argument("--steps-per-iteration", type=int, default=500)
    args = parser.parse_args()

    manifest = json.loads((args.system_dir / "system_manifest.json").read_text())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference_name = "aqueous_system.xml" if args.phase == "aqueous" else "gas_system.xml"
    reference = openmm.XmlSerializer.deserialize((args.system_dir / reference_name).read_text())
    alchemical_system = create_alchemical_system(reference, int(manifest["n_solute_atoms"]))
    (args.output_dir / f"{args.phase}_alchemical_system.xml").write_text(openmm.XmlSerializer.serialize(alchemical_system))

    equilibrated = openmm.XmlSerializer.deserialize((args.system_dir / "equilibrated_state.xml").read_text())
    if args.phase == "aqueous":
        positions = equilibrated.getPositions()
        box_vectors = equilibrated.getPeriodicBoxVectors()
    else:
        positions = equilibrated.getPositions()[: int(manifest["n_solute_atoms"])]
        box_vectors = None
    sampler_state = states.SamplerState(positions=positions, box_vectors=box_vectors)
    thermodynamic_states = create_states(alchemical_system, args.phase)

    platform, properties = cuda_platform()
    cache.global_context_cache = cache.ContextCache(platform=platform, platform_properties=properties, capacity=None)
    move = mcmc.LangevinDynamicsMove(
        timestep=TIMESTEP,
        collision_rate=1.0 / unit.picosecond,
        n_steps=args.steps_per_iteration,
        reassign_velocities=False,
        n_restart_attempts=6,
        constraint_tolerance=1.0e-6,
    )
    sampler = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.iterations)
    storage_path = args.output_dir / f"{args.phase}.nc"
    reporter = multistate.MultiStateReporter(
        str(storage_path),
        checkpoint_interval=10,
        analysis_particle_indices=list(range(int(manifest["n_solute_atoms"]))),
        position_interval=10,
        velocity_interval=0,
    )
    start = time.time()
    sampler.create(thermodynamic_states=thermodynamic_states, sampler_states=sampler_state, storage=reporter)
    sampler.equilibrate(args.equil_iterations)
    sampler.run()
    reporter.close()

    reporter = multistate.MultiStateReporter(str(storage_path), open_mode="r")
    report = analyze(reporter, args.phase)
    reporter.close()
    report.update({
        "compound_id": manifest["compound_id"],
        "equilibration_iterations": args.equil_iterations,
        "production_iterations": args.iterations,
        "steps_per_iteration": args.steps_per_iteration,
        "aggregate_production_ns": args.iterations * args.steps_per_iteration * float(TIMESTEP.value_in_unit(unit.picosecond)) * len(thermodynamic_states) / 1000.0,
        "wall_seconds": time.time() - start,
        "experimental_target_accessed": False,
    })
    report["overlap_pass"] = bool(report["minimum_bidirectional_adjacent_overlap"] >= 0.03)
    (args.output_dir / f"{args.phase}_result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
