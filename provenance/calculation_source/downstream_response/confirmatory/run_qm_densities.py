#!/usr/bin/env python3
"""Acquire confirmatory base-response and isolated W4--W12 QM densities."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "8")

import cupy as cp  # noqa: E402
import gpu4pyscf  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from ase.io import read  # noqa: E402
from pyscf import dft, scf  # noqa: E402

COMMON = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(COMMON))
from common import DEBYE_PER_E_BOHR, molecule, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[4]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response_confirmatory"
METHOD = "DF-RKS omegaB97X-D3(BJ)/def2-TZVPD grid4 SCF 1e-10 GPU4PySCF-1.8.1"


def gpu_scf(molecule_value):
    mean_field = dft.RKS(molecule_value, xc="wb97x-d3bj").density_fit().to_gpu()
    mean_field.grids.level = 4
    mean_field.conv_tol = 1.0e-10
    mean_field.max_cycle = 160
    mean_field.init_guess = "minao"
    started = time.perf_counter()
    energy = float(mean_field.kernel())
    cp.cuda.runtime.deviceSynchronize()
    runtime = time.perf_counter() - started
    if not mean_field.converged:
        raise RuntimeError("SCF did not converge under the frozen protocol")
    density = mean_field.make_rdm1()
    if hasattr(density, "get"):
        density = density.get()
    density = np.asarray(density)
    result = {
        "energy_hartree": energy,
        "density": density,
        "dipole_debye": np.asarray(
            scf.hf.dip_moment(molecule_value, density, unit="Debye", verbose=0)
        ),
        "runtime_seconds": runtime,
        "scf_cycles": int(getattr(mean_field, "cycles", -1)),
    }
    del mean_field
    gc.collect()
    cp.get_default_memory_pool().free_all_blocks()
    return result


def tag(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()[:24]


def acquire_base(base, output: Path) -> dict[str, object]:
    case_id = str(base.info["case_id"])
    array_path = output / f"{tag(case_id)}.npz"
    record_path = output / f"{tag(case_id)}.json"
    if array_path.exists() and record_path.exists():
        return json.loads(record_path.read_text())

    symbols = list(base.get_chemical_symbols())
    positions = np.asarray(base.positions)
    n_solute = int(base.info["n_solute_atoms"])
    solute = np.arange(len(base)) < n_solute
    masks = (np.ones(len(base), dtype=bool), solute, ~solute)
    component_names = (
        "full_base_complex",
        "solute_with_water_ghosts",
        "waters_with_solute_ghosts",
    )
    component_results = []
    for name, mask in zip(component_names, masks, strict=True):
        print(f"  {case_id}: {name}", flush=True)
        component_results.append(gpu_scf(molecule(symbols, positions, mask)))

    component_densities = np.stack([result["density"] for result in component_results])
    response_density = component_densities[0] - component_densities[1] - component_densities[2]
    component_dipoles = np.stack([result["dipole_debye"] for result in component_results])
    response_dipole = component_dipoles[0] - component_dipoles[1] - component_dipoles[2]

    full_molecule = molecule(symbols, positions)
    overlap = full_molecule.intor("int1e_ovlp")
    response_electron_number = float(np.einsum("ij,ji->", response_density, overlap))
    position_integrals = full_molecule.intor("int1e_r", comp=3)
    density_dipole = -np.einsum("xij,ji->x", position_integrals, response_density)
    density_dipole_debye = density_dipole * DEBYE_PER_E_BOHR
    dipole_discrepancy = float(np.max(np.abs(density_dipole_debye - response_dipole)))
    if abs(response_electron_number) > 1.0e-6:
        raise RuntimeError(f"response charge check failed for {case_id}: {response_electron_number}")
    if dipole_discrepancy > 1.0e-5:
        raise RuntimeError(
            f"response-density dipole check failed for {case_id}: {dipole_discrepancy:.3e} D"
        )

    temporary = array_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        base_symbols=np.asarray(symbols),
        base_positions_angstrom=positions,
        component_density_matrices=component_densities,
        response_density_matrix=response_density,
        component_energy_hartree=np.asarray([result["energy_hartree"] for result in component_results]),
        component_dipole_debye=component_dipoles,
        response_dipole_debye=response_dipole,
        response_density_dipole_debye=density_dipole_debye,
        response_electron_number=response_electron_number,
        component_runtime_seconds=np.asarray([result["runtime_seconds"] for result in component_results]),
        component_scf_cycles=np.asarray([result["scf_cycles"] for result in component_results]),
    )
    temporary.replace(array_path)
    record = {
        "case_id": case_id,
        "molecule_id": str(base.info["molecule_id"]),
        "family": str(base.info["family"]),
        "source_frame_index": int(base.info["source_frame_index"]),
        "n_base_atoms": len(base),
        "n_solute_atoms": n_solute,
        "base_nao": int(response_density.shape[0]),
        "method": METHOD,
        "response_electron_number": response_electron_number,
        "response_density_dipole_max_abs_discrepancy_debye": dipole_discrepancy,
        "component_runtime_seconds": [result["runtime_seconds"] for result in component_results],
        "component_scf_cycles": [result["scf_cycles"] for result in component_results],
        "density_file": array_path.name,
        "density_file_sha256": sha256(array_path),
        "all_outer_probes_in_base_response": False,
        "converged": True,
        "prospective_claim": False,
    }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def acquire_probe(probe, output: Path) -> dict[str, object]:
    probe_id = str(probe.info["probe_id"])
    array_path = output / f"{tag(probe_id)}.npz"
    record_path = output / f"{tag(probe_id)}.json"
    if array_path.exists() and record_path.exists():
        return json.loads(record_path.read_text())
    symbols = list(probe.get_chemical_symbols())
    positions = np.asarray(probe.positions)
    print(f"  {probe_id}: isolated probe", flush=True)
    result = gpu_scf(molecule(symbols, positions))
    temporary = array_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        probe_symbols=np.asarray(symbols),
        probe_positions_angstrom=positions,
        probe_density_matrix=result["density"],
        probe_energy_hartree=result["energy_hartree"],
        probe_dipole_debye=result["dipole_debye"],
        runtime_seconds=result["runtime_seconds"],
        scf_cycles=result["scf_cycles"],
    )
    temporary.replace(array_path)
    record = {
        "probe_id": probe_id,
        "case_id": str(probe.info["case_id"]),
        "molecule_id": str(probe.info["molecule_id"]),
        "water_rank": int(probe.info["water_rank"]),
        "source_water_index": int(probe.info["source_water_index"]),
        "oxygen_distance_A": float(probe.info["oxygen_distance_A"]),
        "probe_nao": int(result["density"].shape[0]),
        "method": METHOD,
        "runtime_seconds": result["runtime_seconds"],
        "scf_cycles": result["scf_cycles"],
        "density_file": array_path.name,
        "density_file_sha256": sha256(array_path),
        "converged": True,
        "prospective_claim": False,
    }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def runtime_manifest() -> dict[str, object]:
    gpu = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=name,driver_version,memory.total",
            "--format=csv,noheader",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return {
        "method": METHOD,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pyscf_version": importlib.metadata.version("pyscf"),
        "gpu4pyscf_version": gpu4pyscf.__version__,
        "cupy_version": cp.__version__,
        "gpu": gpu,
        "protocol": "results/posthoc_downstream_response_confirmatory/protocol_freeze.json",
        "prospective_claim": False,
    }


def write_base_manifest(output: Path, records: list[dict[str, object]], geometry: Path) -> None:
    registry = output / "density_registry.csv"
    pd.DataFrame(records).sort_values("case_id").to_csv(registry, index=False)
    manifest = {
        **runtime_manifest(),
        "stage": "base_response_densities",
        "n_cases": len(records),
        "geometry_sha256": sha256(geometry),
        "density_registry_sha256": sha256(registry),
        "all_converged": all(bool(record["converged"]) for record in records),
        "maximum_abs_response_electron_number": max(abs(float(record["response_electron_number"])) for record in records),
        "maximum_dipole_reconstruction_discrepancy_debye": max(float(record["response_density_dipole_max_abs_discrepancy_debye"]) for record in records),
        "all_outer_probes_in_base_response": False,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_probe_manifest(output: Path, records: list[dict[str, object]], geometry: Path) -> None:
    registry = output / "density_registry.csv"
    pd.DataFrame(records).sort_values(["case_id", "water_rank"]).to_csv(registry, index=False)
    manifest = {
        **runtime_manifest(),
        "stage": "isolated_outer_probe_densities",
        "n_probes": len(records),
        "geometry_sha256": sha256(geometry),
        "density_registry_sha256": sha256(registry),
        "all_converged": all(bool(record["converged"]) for record in records),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bases", type=Path, default=RESULT_ROOT / "configurations/base_3water.extxyz")
    parser.add_argument("--probes", type=Path, default=RESULT_ROOT / "configurations/outer_w4_w12.extxyz")
    parser.add_argument("--output", type=Path, default=RESULT_ROOT / "qm")
    parser.add_argument("--stage", choices=("base", "probe", "all"), default="all")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    bases = sorted(read(args.bases, index=":"), key=lambda item: str(item.info["case_id"]))
    probes = sorted(read(args.probes, index=":"), key=lambda item: (str(item.info["case_id"]), int(item.info["water_rank"])))
    if len(bases) != 10 or len(probes) != 90:
        raise RuntimeError("frozen confirmatory protocol requires 10 bases and 90 probes")

    if args.stage in ("base", "all"):
        directory = args.output / "base_densities"
        directory.mkdir(parents=True, exist_ok=True)
        stop = len(bases) if args.limit is None else min(len(bases), args.start + args.limit)
        for index in range(args.start, stop):
            record = acquire_base(bases[index], directory)
            print(f"base [{index + 1}/{len(bases)}] {record['case_id']}", flush=True)
        records = []
        for base in bases:
            path = directory / f"{tag(str(base.info['case_id']))}.json"
            if path.exists():
                records.append(json.loads(path.read_text()))
        if len(records) == len(bases):
            write_base_manifest(directory, records, args.bases)

    if args.stage in ("probe", "all"):
        directory = args.output / "probe_densities"
        directory.mkdir(parents=True, exist_ok=True)
        stop = len(probes) if args.limit is None else min(len(probes), args.start + args.limit)
        for index in range(args.start, stop):
            record = acquire_probe(probes[index], directory)
            print(f"probe [{index + 1}/{len(probes)}] {record['probe_id']}", flush=True)
        records = []
        for probe in probes:
            path = directory / f"{tag(str(probe.info['probe_id']))}.json"
            if path.exists():
                records.append(json.loads(path.read_text()))
        if len(records) == len(probes):
            write_probe_manifest(directory, records, args.probes)


if __name__ == "__main__":
    main()
