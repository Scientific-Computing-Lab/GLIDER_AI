#!/usr/bin/env python3
"""Compute counterpoise interaction energies for archived water-panel geometries.

The original response NPZ files do not store total component energies. This
energy-only calculation uses the same DFT functional, basis, grid, convergence
criterion and D3(BJ) correction as the response reference calculation. It does
not change any geometry, response reference, prediction or primary score.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import ase
import dftd3
import numpy as np
import pyscf
from ase.io import read

from acquire_nonwater_contact_qm_cpu import COMPONENTS, solve


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "experiments" / "water_contact_audit" / "energies"
HARTREE_TO_KCAL_MOL = 627.5094740631
PANELS = ("panel_1", "panel_2", "panel_3")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def calculate(atoms: ase.Atoms, source: Path) -> dict[str, object]:
    n_solute = int(atoms.info["n_solute_atoms"])
    n_atoms = len(atoms)
    masks = (
        np.ones(n_atoms, dtype=bool),
        np.arange(n_atoms) < n_solute,
        np.arange(n_atoms) >= n_solute,
    )
    names = COMPONENTS
    energies = {}
    attempts = {}
    for name, mask in zip(names, masks):
        mean_field, _density, record = solve(atoms, mask, name)
        energies[name] = float(mean_field.e_tot)
        attempts[name] = record
    cp_hartree = energies[names[0]] - energies[names[1]] - energies[names[2]]
    return {
        "config_id": str(atoms.info["config_id"]),
        "panel_source": str(source.relative_to(ROOT)),
        "geometry_sha256": sha256(source),
        "component_energies_hartree": energies,
        "cp_interaction_energy_hartree": cp_hartree,
        "cp_interaction_energy_kcal_mol": cp_hartree * HARTREE_TO_KCAL_MOL,
        "method": "DF-RKS omegaB97X-D3(BJ)/def2-TZVPD grid4 SCF 1e-10",
        "pyscf_version": pyscf.__version__,
        "dftd3_version": dftd3.__version__,
        "ase_version": ase.__version__,
        "attempts": attempts,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="Run all 224 archived configurations")
    parser.add_argument("--limit", type=int, help="Maximum new configurations this invocation")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    pyscf.lib.num_threads(args.threads)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    newly_run = 0
    expected: list[tuple[str, str, str, float]] = []
    for panel in PANELS:
        folder = ROOT / "experiments" / panel
        registry = folder / "geometries" / "configuration_registry.csv"
        with registry.open(newline="") as stream:
            geometry = {row["config_id"]: row for row in csv.DictReader(stream)}
        source = folder / "geometries" / "configurations.extxyz"
        for atoms in read(source, index=":"):
            config_id = str(atoms.info["config_id"])
            if not args.all and float(
                geometry[config_id]["minimum_solute_water_distance_A"]
            ) >= 1.5:
                continue
            expected.append((panel, config_id, str(atoms.info["regime"]),
                             float(geometry[config_id]["minimum_solute_water_distance_A"])))
            output = OUTPUT / f"{hashlib.sha256(config_id.encode()).hexdigest()[:20]}.json"
            if output.exists():
                existing = json.loads(output.read_text())
                if existing["config_id"] != config_id or existing["geometry_sha256"] != sha256(source):
                    raise ValueError(f"Stale output: {output}")
                continue
            if args.limit is not None and newly_run >= args.limit:
                return
            print(f"Calculating {panel}: {config_id}", flush=True)
            result = calculate(atoms, source)
            output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
            newly_run += 1
            print(
                f"CP interaction energy: {result['cp_interaction_energy_kcal_mol']:.3f} kcal/mol",
                flush=True,
            )
    summary_rows = []
    for panel, config_id, regime, distance in expected:
        output = OUTPUT / f"{hashlib.sha256(config_id.encode()).hexdigest()[:20]}.json"
        if not output.exists():
            raise RuntimeError(f"Missing calculated component energies for {config_id}")
        record = json.loads(output.read_text())
        components = record["component_energies_hartree"]
        summary_rows.append({
            "panel": panel,
            "config_id": config_id,
            "regime": regime,
            "minimum_distance_A": distance,
            "cp_interaction_energy_kcal_mol": record["cp_interaction_energy_kcal_mol"],
            "full_complex_energy_hartree": components[COMPONENTS[0]],
            "solute_ghost_energy_hartree": components[COMPONENTS[1]],
            "waters_ghost_energy_hartree": components[COMPONENTS[2]],
        })
    summary_path = OUTPUT.parent / (
        "cp_interaction_energies_all.csv" if args.all else "cp_interaction_energies_short.csv"
    )
    with summary_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Verified {len(summary_rows)} complete CP energy records: {summary_path}")


if __name__ == "__main__":
    main()
