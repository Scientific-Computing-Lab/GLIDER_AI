#!/usr/bin/env python3
"""Compute new CP-consistent response references for a frozen contact panel.

This CPU implementation uses the manuscript's density-functional, basis,
grid, convergence tolerance and non-self-consistent D3(BJ) correction.  The
official frozen GLIDER predictions provide only the geometric probe points.
New response labels are never used to select configurations or a checkpoint.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import ase
import dftd3
import numpy as np
import pyscf
from ase.io import read
from dftd3.interface import DispersionModel, RationalDampingParam
from pyscf import df, dft, gto
from pyscf.scf import hf
from pyscf.scf.addons import project_dm_nr2nr


ROOT = Path(__file__).resolve().parents[2]
BOHR_TO_ANGSTROM = 0.529177210903
COMPONENTS = ("full_complex", "solute_with_neighbour_ghosts", "neighbour_with_solute_ghosts")
DAMPING = RationalDampingParam(method="wb97x")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def make_molecule(atoms: ase.Atoms, real: np.ndarray) -> gto.Mole:
    entries = [
        (symbol if bool(real[i]) else f"ghost-{symbol}", tuple(map(float, xyz)))
        for i, (symbol, xyz) in enumerate(zip(atoms.get_chemical_symbols(), atoms.positions))
    ]
    return gto.M(
        atom=entries, basis="def2-tzvpd", unit="Angstrom", charge=0, spin=0,
        verbose=0, max_memory=4000,
    )


def d3bj_correction(atoms: ase.Atoms, real: np.ndarray) -> float:
    model = DispersionModel(
        np.asarray(atoms.numbers[real], dtype=int),
        np.asarray(atoms.positions[real], dtype=float) / BOHR_TO_ANGSTROM,
    )
    return float(model.get_dispersion(DAMPING, grad=False)["energy"])


def solve(atoms: ase.Atoms, real: np.ndarray, name: str, initial_density=None):
    molecule = make_molecule(atoms, real)
    d3 = d3bj_correction(atoms, real)
    previous_density = initial_density
    attempts = []
    if name == COMPONENTS[1] and np.any(real & (atoms.numbers == 53)):
        # The ordinary MINAO start can occupy the wrong iodine state.  A
        # fragment-only atomic guess supplies an initial density, which is
        # projected into the full ghost basis.  The reported energy and
        # observables still come solely from the unchanged ghost-basis SCF.
        isolated = atoms[real]
        seed_molecule = make_molecule(isolated, np.ones(len(isolated), dtype=bool))
        seed = dft.RKS(seed_molecule, xc="wb97x-d3bj").density_fit()
        seed.grids.level = 4
        seed.conv_tol = 1e-10
        seed.max_cycle = 160
        seed.init_guess = "atom"
        seed.get_dispersion = lambda *args, **kwargs: d3
        print(f"  {name}: isolated atomic-density seed", flush=True)
        started = time.perf_counter()
        seed_energy = float(seed.kernel())
        seed_elapsed = time.perf_counter() - started
        attempts.append({
            "component": name,
            "solver": "isolated_atomic_projection_seed",
            "role": "initial_density_only; not a reported reference component",
            "energy_hartree": seed_energy,
            "converged": bool(seed.converged),
            "cycles": int(getattr(seed, "cycles", -1)),
            "runtime_seconds": seed_elapsed,
        })
        if seed.converged:
            previous_density = project_dm_nr2nr(
                seed_molecule, np.asarray(seed.make_rdm1()), molecule
            )
        print(
            f"  {name}: seed converged={bool(seed.converged)} "
            f"elapsed={seed_elapsed:.1f}s", flush=True,
        )
    for maximum_cycles in (160, 320):
        mean_field = dft.RKS(molecule, xc="wb97x-d3bj").density_fit()
        mean_field.grids.level = 4
        mean_field.conv_tol = 1e-10
        mean_field.level_shift = 0.0
        mean_field.init_guess = "minao"
        mean_field.max_cycle = maximum_cycles
        # PySCF's D3(BJ) interface is not distributed for macOS ARM.  The
        # official s-dftd3 library supplies the identical two-body correction.
        # D3(BJ) changes total energy but is not part of the SCF potential.
        mean_field.get_dispersion = lambda *args, **kwargs: d3
        started = time.perf_counter()
        print(f"  {name}: CPU DIIS, at most {maximum_cycles} cycles", flush=True)
        energy = float(mean_field.kernel(dm0=previous_density))
        elapsed = time.perf_counter() - started
        previous_density = np.asarray(mean_field.make_rdm1())
        attempts.append({
            "component": name,
            "solver": f"cpu_diis_{maximum_cycles}",
            "energy_hartree": energy,
            "d3bj_correction_hartree": d3,
            "converged": bool(mean_field.converged),
            "cycles": int(getattr(mean_field, "cycles", -1)),
            "runtime_seconds": elapsed,
        })
        print(
            f"  {name}: converged={bool(mean_field.converged)} "
            f"energy={energy:.9f} Eh elapsed={elapsed:.1f}s",
            flush=True,
        )
        if mean_field.converged:
            return mean_field, previous_density, attempts
    # Match the original panel's final same-Hamiltonian recovery stage.  Only
    # the numerical optimizer changes; the XC functional, basis and grid do not.
    base = dft.RKS(molecule, xc="wb97x-d3bj").density_fit()
    base.grids.level = 4
    base.conv_tol = 1e-10
    base.level_shift = 0.0
    base.get_dispersion = lambda *args, **kwargs: d3
    mean_field = base.newton()
    mean_field.max_cycle = 200
    mean_field.conv_tol = 1e-10
    mean_field.get_dispersion = lambda *args, **kwargs: d3
    print(f"  {name}: CPU second-order, at most 200 cycles", flush=True)
    started = time.perf_counter()
    energy = float(mean_field.kernel(dm0=previous_density))
    elapsed = time.perf_counter() - started
    previous_density = np.asarray(mean_field.make_rdm1())
    attempts.append({
        "component": name,
        "solver": "cpu_second_order_200",
        "energy_hartree": energy,
        "d3bj_correction_hartree": d3,
        "converged": bool(mean_field.converged),
        "cycles": int(getattr(mean_field, "cycles", -1)),
        "runtime_seconds": elapsed,
    })
    print(
        f"  {name}: converged={bool(mean_field.converged)} "
        f"energy={energy:.9f} Eh elapsed={elapsed:.1f}s", flush=True,
    )
    if mean_field.converged:
        return mean_field, previous_density, attempts
    raise RuntimeError(f"{name} did not converge under the unchanged Hamiltonian")


def esp(molecule: gto.Mole, density: np.ndarray, points: np.ndarray) -> np.ndarray:
    coords = np.asarray(points, dtype=float) / BOHR_TO_ANGSTROM
    nuclear = np.zeros(len(coords))
    for atom in range(molecule.natm):
        displacement = molecule.atom_coord(atom) - coords
        nuclear += molecule.atom_charge(atom) / np.linalg.norm(displacement, axis=1)
    electronic = np.empty(len(coords))
    for start in range(0, len(coords), 64):
        stop = min(start + 64, len(coords))
        fake = gto.fakemol_for_charges(coords[start:stop])
        integrals = df.incore.aux_e2(molecule, fake)
        electronic[start:stop] = np.einsum("ijp,ij->p", integrals, density)
    return nuclear - electronic


def compute_one(atoms: ase.Atoms, prediction_file: Path, output: Path) -> dict[str, object]:
    config_id = str(atoms.info["config_id"])
    tag = hashlib.sha256(config_id.encode()).hexdigest()[:20]
    array_path = output / f"{tag}.npz"
    record_path = output / f"{tag}.json"
    if array_path.exists() and record_path.exists():
        record = json.loads(record_path.read_text())
        if record["config_id"] == config_id and record["observable_sha256"] == sha256(array_path):
            return record
        raise ValueError(f"Incomplete or changed existing reference for {config_id}")
    if array_path.exists() or record_path.exists():
        raise ValueError(f"Half-written reference for {config_id}")
    with np.load(prediction_file, allow_pickle=False) as prediction:
        points = np.asarray(prediction["points_angstrom"], dtype=float)
    n_solute = int(atoms.info["n_solute_atoms"])
    solute = np.arange(len(atoms)) < n_solute
    neighbour = ~solute
    full = np.ones(len(atoms), dtype=bool)
    solute_mf, solute_dm, solute_attempts = solve(atoms, solute, COMPONENTS[1])
    neighbour_mf, neighbour_dm, neighbour_attempts = solve(atoms, neighbour, COMPONENTS[2])
    full_mf, full_dm, full_attempts = solve(
        atoms, full, COMPONENTS[0], initial_density=solute_dm + neighbour_dm
    )
    mean_fields = [full_mf, solute_mf, neighbour_mf]
    densities = [full_dm, solute_dm, neighbour_dm]
    dipoles = np.asarray([
        hf.dip_moment(mf.mol, density, unit="Debye", verbose=0)
        for mf, density in zip(mean_fields, densities)
    ])
    potentials = np.asarray([
        esp(mf.mol, density, points)
        for mf, density in zip(mean_fields, densities)
    ])
    energies = np.asarray([float(mf.e_tot) for mf in mean_fields])
    response_dipole = dipoles[0] - dipoles[1] - dipoles[2]
    response_esp = potentials[0] - potentials[1] - potentials[2]
    temporary = array_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        points_angstrom=points,
        component_energy_hartree=energies,
        component_dipole_debye=dipoles,
        component_esp_hartree_per_e=potentials,
        delta_dipole_debye=response_dipole,
        delta_esp_hartree_per_e=response_esp,
    )
    temporary.replace(array_path)
    attempts = full_attempts + solute_attempts + neighbour_attempts
    record: dict[str, object] = {
        "config_id": config_id,
        "molecule_id": str(atoms.info["molecule_id"]),
        "perturbant": str(atoms.info["perturbant"]),
        "component_order": list(COMPONENTS),
        "component_energy_hartree": energies.tolist(),
        "cp_interaction_energy_kcal_mol": float(
            (energies[0] - energies[1] - energies[2]) * 627.509474
        ),
        "all_components_converged": True,
        "attempts": attempts,
        "prediction_probe_sha256": sha256(prediction_file),
        "observable_file": array_path.name,
        "observable_sha256": sha256(array_path),
        "response_esp_rms_mEh_per_e": float(np.sqrt(np.mean(response_esp**2)) * 1000),
        "response_dipole_norm_debye": float(np.linalg.norm(response_dipole)),
    }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", type=Path, default=ROOT / "build/nonwater_contact_panel")
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    manifest_path = args.experiment / "prediction_freeze.json"
    if not manifest_path.exists():
        raise RuntimeError("Refusing QM before a prediction-freeze manifest exists")
    manifest = json.loads(manifest_path.read_text())
    geometry_file = args.experiment / "configurations.extxyz"
    if manifest["geometry_sha256"] != sha256(geometry_file):
        raise RuntimeError("Frozen geometry digest changed")
    for method, key in (("glider", "glider_predictions"),
                        ("mace_polar_l", "mace_polar_l_predictions")):
        method_dir = args.experiment / "predictions" / method
        registry_file = method_dir / "prediction_registry.csv"
        frozen = manifest[key]
        if sha256(registry_file) != frozen["registry_sha256"]:
            raise RuntimeError(f"Frozen {method} registry changed")
        current = {row["config_id"]: row for row in read_csv(registry_file)}
        if set(current) != set(frozen["case_file_sha256"]):
            raise RuntimeError(f"Frozen {method} case set changed")
        for config_id, row in current.items():
            if (row["prediction_sha256"] != frozen["case_file_sha256"][config_id]
                    or sha256(method_dir / row["prediction_file"])
                    != frozen["case_file_sha256"][config_id]):
                raise RuntimeError(f"Frozen {method} prediction changed: {config_id}")
    pred_registry = {
        row["config_id"]: row for row in read_csv(args.predictions / "prediction_registry.csv")
    }
    frames = read(geometry_file, index=":")
    if len(frames) != 36 or set(pred_registry) != {str(atoms.info["config_id"]) for atoms in frames}:
        raise RuntimeError("Predictions do not cover all frozen geometries")
    output = args.experiment / "references"
    output.mkdir(exist_ok=True)
    stop = len(frames) if args.limit is None else min(len(frames), args.start + args.limit)
    for index in range(args.start, stop):
        atoms = frames[index]
        config_id = str(atoms.info["config_id"])
        entry = pred_registry[config_id]
        prediction_file = args.predictions / entry["prediction_file"]
        if sha256(prediction_file) != entry["prediction_sha256"]:
            raise RuntimeError(f"Prediction digest changed for {config_id}")
        started = time.perf_counter()
        record = compute_one(atoms, prediction_file, output)
        print(
            f"[{index+1}/{len(frames)}] {config_id} "
            f"CP={record['cp_interaction_energy_kcal_mol']:.3f} kcal/mol "
            f"elapsed={time.perf_counter()-started:.1f}s",
            flush=True,
        )
    finished = [
        json.loads(path.read_text())
        for path in output.glob("*.json")
        if path.name != "REFERENCE_MANIFEST.json"
    ]
    if len(finished) == len(frames):
        rows = [{
            "config_id": row["config_id"],
            "molecule_id": row["molecule_id"],
            "perturbant": row["perturbant"],
            "cp_interaction_energy_kcal_mol": row["cp_interaction_energy_kcal_mol"],
            "response_esp_rms_mEh_per_e": row["response_esp_rms_mEh_per_e"],
            "observable_file": row["observable_file"],
            "observable_sha256": row["observable_sha256"],
        } for row in sorted(finished, key=lambda item: item["config_id"])]
        with (output / "REFERENCE_REGISTRY.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (output / "REFERENCE_MANIFEST.json").write_text(json.dumps({
            "status": "all 36 new response references computed after prediction freeze",
            "n_configurations": len(rows),
            "geometry_sha256": sha256(geometry_file),
            "prediction_freeze_sha256": sha256(manifest_path),
            "reference_registry_sha256": sha256(output / "REFERENCE_REGISTRY.csv"),
            "method": "DF-RKS omegaB97X-D3(BJ)/def2-TZVPD grid4 SCF 1e-10",
            "pyscf_version": pyscf.__version__,
            "dftd3_version": dftd3.__version__,
            "ase_version": ase.__version__,
        }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
