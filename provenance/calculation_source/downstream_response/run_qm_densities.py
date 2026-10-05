#!/usr/bin/env python3
"""Reconstruct the frozen base response and isolated W4 density matrices."""

from __future__ import annotations

import argparse
import csv
import gc
import importlib.metadata
import json
import os
import platform
import subprocess
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "8")

import cupy as cp  # noqa: E402
import gpu4pyscf  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from ase.io import read  # noqa: E402
from common import molecule, sha256  # noqa: E402
from pyscf import dft, scf  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
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
    result = {
        "energy_hartree": energy,
        "density": np.asarray(density),
        "dipole_debye": np.asarray(
            scf.hf.dip_moment(
                molecule_value, np.asarray(density), unit="Debye", verbose=0
            )
        ),
        "runtime_seconds": runtime,
        "scf_cycles": int(getattr(mean_field, "cycles", -1)),
    }
    del mean_field
    gc.collect()
    cp.get_default_memory_pool().free_all_blocks()
    return result


def reference_map() -> dict[str, Path]:
    root = ROOT / "evidence/liquid_bridge/references"
    with (root / "REFERENCE_REGISTRY.csv").open() as handle:
        return {
            row["config_id"]: root / row["observable_file"]
            for row in csv.DictReader(handle)
        }


def acquire_case(base, probe, reference_path: Path, output: Path) -> dict[str, object]:
    case_id = str(base.info["case_id"])
    base_config_id = str(base.info["base_config_id"])
    tag = __import__("hashlib").sha256(case_id.encode()).hexdigest()[:20]
    array_path = output / f"{tag}.npz"
    record_path = output / f"{tag}.json"
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

    probe_symbols = list(probe.get_chemical_symbols())
    probe_positions = np.asarray(probe.positions)
    print(f"  {case_id}: isolated_W4", flush=True)
    probe_result = gpu_scf(molecule(probe_symbols, probe_positions))

    component_densities = np.stack(
        [result["density"] for result in component_results]
    )
    response_density = (
        component_densities[0] - component_densities[1] - component_densities[2]
    )
    component_dipoles = np.stack(
        [result["dipole_debye"] for result in component_results]
    )
    response_dipole = component_dipoles[0] - component_dipoles[1] - component_dipoles[2]
    with np.load(reference_path, allow_pickle=False) as stored:
        stored_response_dipole = np.asarray(stored["delta_dipole_debye"])
    dipole_discrepancy = float(
        np.max(np.abs(response_dipole - stored_response_dipole))
    )
    if dipole_discrepancy > 1.0e-5:
        raise RuntimeError(
            f"recomputed response dipole disagrees with stored reference by "
            f"{dipole_discrepancy:.3e} D: {case_id}"
        )

    temporary = array_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        base_symbols=np.asarray(symbols),
        base_positions_angstrom=positions,
        component_density_matrices=component_densities,
        response_density_matrix=response_density,
        component_energy_hartree=np.asarray(
            [result["energy_hartree"] for result in component_results]
        ),
        component_dipole_debye=component_dipoles,
        response_dipole_debye=response_dipole,
        stored_response_dipole_debye=stored_response_dipole,
        probe_symbols=np.asarray(probe_symbols),
        probe_positions_angstrom=probe_positions,
        probe_density_matrix=probe_result["density"],
        probe_energy_hartree=probe_result["energy_hartree"],
        probe_dipole_debye=probe_result["dipole_debye"],
        component_runtime_seconds=np.asarray(
            [result["runtime_seconds"] for result in component_results]
        ),
        component_scf_cycles=np.asarray(
            [result["scf_cycles"] for result in component_results]
        ),
        probe_runtime_seconds=probe_result["runtime_seconds"],
        probe_scf_cycles=probe_result["scf_cycles"],
    )
    temporary.replace(array_path)
    record = {
        "case_id": case_id,
        "base_config_id": base_config_id,
        "molecule_id": str(base.info["molecule_id"]),
        "source_frame_index": int(base.info["source_frame_index"]),
        "n_base_atoms": len(base),
        "n_solute_atoms": n_solute,
        "base_nao": int(response_density.shape[0]),
        "probe_nao": int(probe_result["density"].shape[0]),
        "method": METHOD,
        "response_dipole_max_abs_discrepancy_debye": dipole_discrepancy,
        "component_runtime_seconds": [
            result["runtime_seconds"] for result in component_results
        ],
        "component_scf_cycles": [
            result["scf_cycles"] for result in component_results
        ],
        "probe_runtime_seconds": probe_result["runtime_seconds"],
        "probe_scf_cycles": probe_result["scf_cycles"],
        "stored_reference": str(reference_path.relative_to(ROOT)),
        "stored_reference_sha256": sha256(reference_path),
        "density_file": array_path.name,
        "density_file_sha256": sha256(array_path),
        "converged": True,
        "W4_in_base_response": False,
        "prospective_claim": False,
    }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def write_manifest(
    output: Path, records: list[dict[str, object]], base_path: Path, probe_path: Path
) -> None:
    table = output / "density_registry.csv"
    pd.DataFrame(records).sort_values("case_id").to_csv(table, index=False)
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
    manifest = {
        "experiment": "post hoc held-out-W4 downstream response",
        "protocol": "results/posthoc_downstream_response/protocol_freeze.json",
        "method": METHOD,
        "n_cases": len(records),
        "base_geometry_sha256": sha256(base_path),
        "probe_geometry_sha256": sha256(probe_path),
        "density_registry_sha256": sha256(table),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pyscf_version": importlib.metadata.version("pyscf"),
        "gpu4pyscf_version": gpu4pyscf.__version__,
        "cupy_version": cp.__version__,
        "gpu": gpu,
        "all_converged": all(bool(record["converged"]) for record in records),
        "maximum_stored_dipole_discrepancy_debye": max(
            float(record["response_dipole_max_abs_discrepancy_debye"])
            for record in records
        ),
        "W4_in_base_response": False,
        "prospective_claim": False,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bases",
        type=Path,
        default=ROOT
        / "results/posthoc_downstream_response/configurations/base_3water.extxyz",
    )
    parser.add_argument(
        "--probes",
        type=Path,
        default=ROOT
        / "results/posthoc_downstream_response/configurations/heldout_w4.extxyz",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results/posthoc_downstream_response/qm_densities",
    )
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    bases = read(args.bases, index=":")
    probes = {str(probe.info["case_id"]): probe for probe in read(args.probes, index=":")}
    references = reference_map()
    if len(bases) != 24 or len(probes) != 24:
        raise RuntimeError("the frozen protocol requires all 24 cases")
    stop = len(bases) if args.limit is None else min(len(bases), args.start + args.limit)
    for index in range(args.start, stop):
        base = bases[index]
        case_id = str(base.info["case_id"])
        base_config_id = str(base.info["base_config_id"])
        record = acquire_case(
            base, probes[case_id], references[base_config_id], args.output
        )
        print(
            f"[{index + 1}/{len(bases)}] {case_id} "
            f"dipole_check={record['response_dipole_max_abs_discrepancy_debye']:.3e} D",
            flush=True,
        )
    records = []
    for base in bases:
        case_id = str(base.info["case_id"])
        tag = __import__("hashlib").sha256(case_id.encode()).hexdigest()[:20]
        record_path = args.output / f"{tag}.json"
        if record_path.exists():
            records.append(json.loads(record_path.read_text()))
    if len(records) == len(bases):
        write_manifest(args.output, records, args.bases, args.probes)


if __name__ == "__main__":
    main()
