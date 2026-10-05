#!/usr/bin/env python3
"""Acquire rigid probe variants for confirmatory torque and orientation endpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
COMMON = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(1, str(COMMON))
from common import sha256  # noqa: E402
from run_qm_densities import METHOD, RESULT_ROOT, gpu_scf, molecule  # noqa: E402

AXES = np.eye(3)
PRIMARY_STEP_DEGREES = 0.1
CONVERGENCE_STEPS_DEGREES = (0.05, 0.2)


def canonical_quaternion(quaternion: np.ndarray) -> np.ndarray:
    value = np.asarray(quaternion, dtype=float).copy()
    for component in value:
        if abs(component) > 1.0e-12:
            if component < 0:
                value *= -1.0
            break
    return value


def octahedral_rotations() -> list[tuple[str, np.ndarray, np.ndarray]]:
    values = []
    for rotation in Rotation.create_group("O"):
        quaternion = canonical_quaternion(rotation.as_quat())
        values.append((quaternion, rotation.as_matrix()))
    values.sort(key=lambda item: tuple(np.round(item[0], 14)))
    return [
        (f"O_{index:02d}", quaternion, matrix)
        for index, (quaternion, matrix) in enumerate(values)
    ]


def rotate_about_oxygen(positions: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    oxygen = np.asarray(positions[0])
    return oxygen + (np.asarray(positions) - oxygen) @ np.asarray(matrix).T


def geometry_invariants(positions: np.ndarray) -> tuple[np.ndarray, float]:
    vectors = np.asarray(positions[1:]) - np.asarray(positions[0])
    lengths = np.linalg.norm(vectors, axis=1)
    cosine = np.dot(vectors[0], vectors[1]) / np.prod(lengths)
    angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
    return lengths, angle


def torque_variants(probe_id: str, positions: np.ndarray, step: float, kind: str):
    for axis_index, axis in enumerate(AXES):
        for sign, label in ((-1.0, "minus"), (1.0, "plus")):
            matrix = Rotation.from_rotvec(axis * np.radians(sign * step)).as_matrix()
            yield {
                "probe_id": probe_id,
                "variant_id": f"torque_{step:g}deg_axis{axis_index}_{label}",
                "variant_type": kind,
                "axis": axis_index,
                "sign": int(sign),
                "step_degrees": step,
                "rotation_matrix": matrix,
                "positions": rotate_about_oxygen(positions, matrix),
            }


def variants_for_probe(probe, positions: np.ndarray):
    probe_id = str(probe.info["probe_id"])
    yield from torque_variants(
        probe_id, positions, PRIMARY_STEP_DEGREES, "torque_primary"
    )
    if int(probe.info["water_rank"]) == 4:
        for step in CONVERGENCE_STEPS_DEGREES:
            yield from torque_variants(
                probe_id, positions, step, "torque_convergence"
            )
        for rotation_id, quaternion, matrix in octahedral_rotations():
            yield {
                "probe_id": probe_id,
                "variant_id": f"orientation_{rotation_id}",
                "variant_type": "orientation",
                "rotation_id": rotation_id,
                "quaternion_xyzw": quaternion,
                "rotation_matrix": matrix,
                "positions": rotate_about_oxygen(positions, matrix),
            }


def load_initial_density_map(directory: Path) -> dict[str, Path]:
    table = pd.read_csv(directory / "density_registry.csv")
    return {
        str(row.probe_id): directory / str(row.density_file)
        for row in table.itertuples()
    }


def variant_tag(probe_id: str, variant_id: str) -> str:
    return hashlib.sha256(f"{probe_id}::{variant_id}".encode()).hexdigest()[:24]


def acquire_variant(
    variant: dict[str, object],
    initial_positions: np.ndarray,
    initial_density_file: Path,
    output: Path,
) -> dict[str, object]:
    probe_id = str(variant["probe_id"])
    variant_id = str(variant["variant_id"])
    tag = variant_tag(probe_id, variant_id)
    array_path = output / f"{tag}.npz"
    record_path = output / f"{tag}.json"
    if array_path.exists() and record_path.exists():
        return json.loads(record_path.read_text())

    positions = np.asarray(variant["positions"], dtype=float)
    original_lengths, original_angle = geometry_invariants(initial_positions)
    lengths, angle = geometry_invariants(positions)
    bond_discrepancy = float(np.max(np.abs(lengths - original_lengths)))
    angle_discrepancy = abs(angle - original_angle)
    if bond_discrepancy > 1.0e-10 or angle_discrepancy > 1.0e-10:
        raise RuntimeError(f"rigid geometry check failed: {probe_id}::{variant_id}")

    identity = np.max(np.abs(np.asarray(variant["rotation_matrix"]) - np.eye(3))) < 1.0e-13
    if identity:
        with np.load(initial_density_file, allow_pickle=False) as initial:
            density = np.asarray(initial["probe_density_matrix"])
            energy = float(initial["probe_energy_hartree"])
            dipole = np.asarray(initial["probe_dipole_debye"])
        runtime = 0.0
        cycles = 0
        density_source = "reused_initial_identity_orientation"
    else:
        result = gpu_scf(molecule(["O", "H", "H"], positions))
        density = np.asarray(result["density"])
        energy = float(result["energy_hartree"])
        dipole = np.asarray(result["dipole_debye"])
        runtime = float(result["runtime_seconds"])
        cycles = int(result["scf_cycles"])
        density_source = "independent_isolated_probe_SCF"

    payload = {
        key: value
        for key, value in variant.items()
        if key not in {"positions", "rotation_matrix", "quaternion_xyzw"}
    }
    temporary = array_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        positions_angstrom=positions,
        density_matrix=density,
        energy_hartree=energy,
        dipole_debye=dipole,
        rotation_matrix=np.asarray(variant["rotation_matrix"]),
        quaternion_xyzw=np.asarray(variant.get("quaternion_xyzw", np.full(4, np.nan))),
    )
    temporary.replace(array_path)
    record = {
        **payload,
        "density_source": density_source,
        "runtime_seconds": runtime,
        "scf_cycles": cycles,
        "bond_length_max_discrepancy_angstrom": bond_discrepancy,
        "angle_discrepancy_degrees": angle_discrepancy,
        "density_file": array_path.name,
        "density_file_sha256": sha256(array_path),
        "converged": True,
    }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probes", type=Path, default=RESULT_ROOT / "configurations/outer_w4_w12.extxyz")
    parser.add_argument("--initial-densities", type=Path, default=RESULT_ROOT / "qm/probe_densities")
    parser.add_argument("--output", type=Path, default=RESULT_ROOT / "qm/probe_variants")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    probes = sorted(read(args.probes, index=":"), key=lambda item: (str(item.info["case_id"]), int(item.info["water_rank"])))
    if len(probes) != 90:
        raise RuntimeError("frozen confirmatory protocol requires all 90 probes")
    initial = load_initial_density_map(args.initial_densities)
    variants = []
    initial_positions = {}
    for probe in probes:
        probe_id = str(probe.info["probe_id"])
        positions = np.asarray(probe.positions)
        initial_positions[probe_id] = positions
        for variant in variants_for_probe(probe, positions):
            variant["case_id"] = str(probe.info["case_id"])
            variant["molecule_id"] = str(probe.info["molecule_id"])
            variant["water_rank"] = int(probe.info["water_rank"])
            variants.append(variant)
    variants.sort(key=lambda item: (str(item["probe_id"]), str(item["variant_id"])))

    stop = len(variants) if args.limit is None else min(len(variants), args.start + args.limit)
    for index in range(args.start, stop):
        variant = variants[index]
        probe_id = str(variant["probe_id"])
        record = acquire_variant(
            variant, initial_positions[probe_id], initial[probe_id], args.output
        )
        print(f"variant [{index + 1}/{len(variants)}] {probe_id}::{record['variant_id']}", flush=True)

    records = []
    for variant in variants:
        path = args.output / f"{variant_tag(str(variant['probe_id']), str(variant['variant_id']))}.json"
        if path.exists():
            records.append(json.loads(path.read_text()))
    if len(records) == len(variants):
        registry = args.output / "variant_registry.csv"
        pd.DataFrame(records).sort_values(["probe_id", "variant_id"]).to_csv(registry, index=False)
        manifest = {
            "experiment": "confirmatory outer-water rigid torque and W4 orientation variants",
            "protocol": "results/posthoc_downstream_response_confirmatory/protocol_freeze.json",
            "method": METHOD,
            "n_variants": len(records),
            "n_primary_torque_variants": sum(record["variant_type"] == "torque_primary" for record in records),
            "n_convergence_torque_variants": sum(record["variant_type"] == "torque_convergence" for record in records),
            "n_orientation_variants": sum(record["variant_type"] == "orientation" for record in records),
            "variant_registry_sha256": sha256(registry),
            "maximum_bond_discrepancy_angstrom": max(float(record["bond_length_max_discrepancy_angstrom"]) for record in records),
            "maximum_angle_discrepancy_degrees": max(float(record["angle_discrepancy_degrees"]) for record in records),
            "all_converged": all(bool(record["converged"]) for record in records),
            "prospective_claim": False,
        }
        (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
