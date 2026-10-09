#!/usr/bin/env python3
"""Compute CP-consistent QM response for the archived cyclic-carbamate separation scan.

This is a post hoc diagnostic. It does not alter the frozen predictions or the
prospective benchmarks. Run from the GLIDER_AI environment and save outputs to
an untracked build directory until the scientific checks are complete.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pyscf
from ase.io import iread
from pyscf.scf import hf

ROOT = Path(__file__).resolve().parents[2]
SOLUTE = "dev_cyclic_carbamate"
DISTANCES = (3, 4, 5, 6, 8, 10, 15, 20, 50, 100)
SYSTEMS = ("single_water", "whole_environment")
METHOD = "DF-RKS omegaB97X-D3(BJ)/def2-TZVPD grid4 SCF 1e-10"
COMPONENTS = ("full_complex", "solute_with_neighbour_ghosts",
              "neighbour_with_solute_ghosts")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def backend_functions(backend: str, gpu_helper: Path | None):
    if backend == "cpu":
        from acquire_nonwater_contact_qm_cpu import esp, solve
        return solve, esp, {"backend": "cpu_pyscf"}
    if gpu_helper is None or not gpu_helper.is_file():
        raise ValueError("--gpu-helper must identify the original GPU reference script")
    spec = importlib.util.spec_from_file_location("gpu_reference_helper", gpu_helper)
    if spec is None or spec.loader is None:
        raise ValueError(f"Cannot load GPU helper: {gpu_helper}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def solve_gpu(atoms, real, name, initial_density=None):
        molecule = module.molecule_with_ghosts(
            list(atoms.get_chemical_symbols()), np.asarray(atoms.positions, dtype=float), real
        )
        return module.solve(molecule, name, initial_density)

    return solve_gpu, module.esp, {
        "backend": "gpu4pyscf",
        "gpu_reference_helper_sha256": sha256(gpu_helper),
        "gpu4pyscf_version": module.gpu4pyscf.__version__,
    }


def source_for(distance: int) -> Path:
    group = "dissociation_extended" if distance >= 50 else "dissociation"
    return ROOT / "experiments" / group


def frame_for(source: Path, config_id: str):
    matches = (atoms for atoms in iread(source / "configurations.extxyz")
               if str(atoms.info.get("config_id")) == config_id)
    atoms = next(matches, None)
    if atoms is None or next(matches, None) is not None:
        raise ValueError(f"Expected exactly one geometry for {config_id}")
    if str(atoms.info["molecule_id"]) != SOLUTE:
        raise ValueError(f"Unexpected molecule for {config_id}")
    return atoms


def calculate(source: Path, config_id: str, output: Path,
              solver: Callable, esp_fn: Callable, backend_meta: dict) -> dict[str, object]:
    atoms = frame_for(source, config_id)
    with np.load(source / "solute_probe_points.npz", allow_pickle=False) as data:
        points = np.asarray(data[SOLUTE], dtype=float)
    n_solute = int(atoms.info["n_solute_atoms"])
    real_solute = np.arange(len(atoms)) < n_solute
    real_environment = ~real_solute
    masks = (real_solute, real_environment, np.ones(len(atoms), dtype=bool))
    names = (COMPONENTS[1], COMPONENTS[2], COMPONENTS[0])
    results = {}
    for name, mask in zip(names, masks):
        initial = (results[names[0]]["density"] + results[names[1]]["density"]
                   if name == COMPONENTS[0] else None)
        mf, density, attempts = solver(atoms, mask, name, initial_density=initial)
        results[name] = {
            "mean_field": mf,
            "density": density,
            "attempts": attempts,
            "energy_hartree": float(mf.e_tot),
            "dipole_debye": np.asarray(
                hf.dip_moment(mf.mol, density, unit="Debye", verbose=0), dtype=float
            ),
            "esp_hartree_per_e": esp_fn(mf.mol, density, points),
        }
    full, solute, environment = (results[name] for name in COMPONENTS)
    delta_esp = (full["esp_hartree_per_e"] - solute["esp_hartree_per_e"]
                 - environment["esp_hartree_per_e"])
    delta_dipole = (full["dipole_debye"] - solute["dipole_debye"]
                    - environment["dipole_debye"])
    with np.load(source / "predicted_probe_potentials.npz", allow_pickle=False) as data:
        glider = np.asarray(data[f"{config_id}__glider"], dtype=float)
        prior = np.asarray(data[f"{config_id}__averaged_prior"], dtype=float)
    if not (len(points) == len(delta_esp) == len(glider) == len(prior)):
        raise ValueError("Probe-point mismatch")
    arrays = output / f"{config_id}.npz"
    tmp = output / f"{config_id}.tmp.npz"
    np.savez_compressed(
        tmp, points_angstrom=points, qm_response_esp_hartree_per_e=delta_esp,
        qm_response_dipole_debye=delta_dipole,
        glider_esp_hartree_per_e=glider,
        averaged_prior_esp_hartree_per_e=prior,
        component_energy_hartree=np.asarray([x["energy_hartree"] for x in
                                             (full, solute, environment)]),
        component_esp_hartree_per_e=np.asarray([x["esp_hartree_per_e"] for x in
                                                (full, solute, environment)]),
    )
    tmp.replace(arrays)
    rms = lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)) * 1000)
    record = {
        "config_id": config_id,
        "method": METHOD,
        "n_atoms": len(atoms),
        "n_solute_atoms": n_solute,
        "n_probe_points": len(points),
        "geometry_source": str((source / "configurations.extxyz").relative_to(ROOT)),
        "geometry_source_sha256": sha256(source / "configurations.extxyz"),
        "probe_source_sha256": sha256(source / "solute_probe_points.npz"),
        "prediction_source_sha256": sha256(source / "predicted_probe_potentials.npz"),
        "qm_response_rms_mEh_per_e": rms(delta_esp),
        "glider_rms_mEh_per_e": rms(glider),
        "prior_rms_mEh_per_e": rms(prior),
        "glider_error_rms_mEh_per_e": rms(glider - delta_esp),
        "prior_error_rms_mEh_per_e": rms(prior - delta_esp),
        "qm_response_dipole_norm_debye": float(np.linalg.norm(delta_dipole)),
        "cp_interaction_energy_kcal_mol": float(
            (full["energy_hartree"] - solute["energy_hartree"]
             - environment["energy_hartree"]) * 627.5094740631
        ),
        "component_energy_hartree": [x["energy_hartree"] for x in
                                     (full, solute, environment)],
        "attempts": {name: results[name]["attempts"] for name in COMPONENTS},
        "arrays_file": arrays.name,
        "arrays_sha256": sha256(arrays),
        "pyscf_version": pyscf.__version__,
        **backend_meta,
    }
    (output / f"{config_id}.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n"
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "build/separation_qm_dev_cyclic_carbamate")
    parser.add_argument("--systems", nargs="+", choices=SYSTEMS, default=SYSTEMS)
    parser.add_argument("--distances", nargs="+", type=int, choices=DISTANCES,
                        default=DISTANCES)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--backend", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--gpu-helper", type=Path,
                        help="Path to the frozen run_reference_qm.py GPU implementation")
    args = parser.parse_args()
    pyscf.lib.num_threads(args.threads)
    solver, esp_fn, backend_meta = backend_functions(args.backend, args.gpu_helper)
    args.output.mkdir(parents=True, exist_ok=True)
    for system in args.systems:
        for distance in args.distances:
            source = source_for(distance)
            config_id = f"dissociation__{SOLUTE}__{system}__{distance}A"
            record_file = args.output / f"{config_id}.json"
            if record_file.exists():
                record = json.loads(record_file.read_text())
                if record["config_id"] != config_id or sha256(
                    args.output / record["arrays_file"]
                ) != record["arrays_sha256"]:
                    raise ValueError(f"Existing output is incomplete: {config_id}")
                print(f"Verified existing {config_id}", flush=True)
            else:
                print(f"Computing {config_id}", flush=True)
                started = time.perf_counter()
                record = calculate(source, config_id, args.output,
                                   solver, esp_fn, backend_meta)
                print(f"Finished {config_id} in {time.perf_counter()-started:.1f}s; "
                      f"QM={record['qm_response_rms_mEh_per_e']:.4f}, "
                      f"GLIDER error={record['glider_error_rms_mEh_per_e']:.4f}, "
                      f"prior error={record['prior_error_rms_mEh_per_e']:.4f} mEh/e",
                      flush=True)
    rows = []
    for system in SYSTEMS:
        for distance in DISTANCES:
            config_id = f"dissociation__{SOLUTE}__{system}__{distance}A"
            path = args.output / f"{config_id}.json"
            if path.exists():
                row = json.loads(path.read_text())
                rows.append({"system": system, "distance_A": distance,
                             **{key: row[key] for key in (
                                 "qm_response_rms_mEh_per_e", "glider_rms_mEh_per_e",
                                 "prior_rms_mEh_per_e", "glider_error_rms_mEh_per_e",
                                 "prior_error_rms_mEh_per_e",
                                 "qm_response_dipole_norm_debye",
                                 "cp_interaction_energy_kcal_mol")}})
    if rows:
        with (args.output / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(f"Complete: {len(rows)}/20 conditions", flush=True)


if __name__ == "__main__":
    main()
