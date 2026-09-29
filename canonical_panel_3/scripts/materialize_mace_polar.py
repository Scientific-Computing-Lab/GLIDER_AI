#!/usr/bin/env python3
"""Apply the frozen faithful Gaussian MACE-POLAR decoder to label-free grids."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def frozen_decoder():
    source = ROOT / "evidence/frozen_source/comparator_verification/evaluate_mace_faithful.py"
    spec = importlib.util.spec_from_file_location("frozen_mace_faithful", source)
    if spec is None or spec.loader is None:
        raise ImportError(source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, required=True)
    parser.add_argument("--grid-predictions", type=Path, required=True)
    parser.add_argument("--response-archive", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    args = parser.parse_args()
    decoder = frozen_decoder()
    frames = read(args.configurations, index=":")
    frame_map = {str(atoms.info["config_id"]): atoms for atoms in frames}
    grid_table = pd.read_csv(args.grid_predictions / "prediction_registry.csv")
    grids = {row.config_id: args.grid_predictions / row.prediction_file for row in grid_table.itertuples()}
    response = np.load(args.response_archive)
    ids = [str(value) for value in response["config_ids"]]
    if set(ids) != set(frame_map) or set(ids) != set(grids):
        raise RuntimeError("Geometry, grid and MACE response sets differ")
    sigma, max_l, kcut, model_type = decoder.load_sigma_and_metadata(args.checkpoint)
    args.output.mkdir(parents=True, exist_ok=False)
    rows, first_values = [], None
    started = time.perf_counter()
    for index, config_id in enumerate(ids):
        atoms = frame_map[config_id]
        first, last = int(response["offsets"][index]), int(response["offsets"][index + 1])
        q = np.asarray(response["charges_e"][first:last], dtype=np.float64)
        u = np.asarray(response["dipoles_e_bohr"][first:last], dtype=np.float64)
        points = np.asarray(np.load(grids[config_id])["points_angstrom"], dtype=np.float64)
        potential = decoder.gaussian_multipole_potential(atoms.positions, points, q, u, sigma)
        dipole = (np.sum(q[:, None] * (atoms.positions / decoder.BOHR_TO_ANGSTROM), axis=0) + np.sum(u, axis=0)) * decoder.DEBYE_PER_E_BOHR
        filename = f"{hashlib.sha256(config_id.encode()).hexdigest()[:20]}.npz"
        path = args.output / filename
        np.savez_compressed(path, points_angstrom=points, predicted_esp_hartree_per_e=potential, predicted_dipole_debye=dipole, predicted_charges_e=q, predicted_dipoles_e_bohr=u)
        rows.append({"config_id": config_id, "molecule_id": str(atoms.info["molecule_id"]), "regime": str(atoms.info["regime"]), "prediction_file": filename, "prediction_sha256": sha256(path), "net_charge_e": float(q.sum())})
        if first_values is None:
            first_values = (atoms.positions.copy(), points.copy(), q.copy(), u.copy())
    registry = args.output / "prediction_registry.csv"
    pd.DataFrame(rows).sort_values("config_id").to_csv(registry, index=False)
    audit = {
        "method": args.method,
        "decoder": "unchanged analytic open-boundary Coulomb integral of official normalized Gaussian l<=1 density",
        "frozen_decoder_sha256": sha256(ROOT / "evidence/frozen_source/comparator_verification/evaluate_mace_faithful.py"),
        "checkpoint_sha256": sha256(args.checkpoint),
        "response_archive_sha256": sha256(args.response_archive),
        "configuration_sha256": sha256(args.configurations),
        "grid_registry_sha256": sha256(args.grid_predictions / "prediction_registry.csv"),
        "density_smearing_width_angstrom": sigma,
        "density_max_l": max_l,
        "kspace_cutoff_inverse_angstrom": kcut,
        "checkpoint_model_type": model_type,
        "n_configurations": len(rows),
        "max_abs_net_response_charge_e": float(max(abs(row["net_charge_e"]) for row in rows)),
        "covariance_check": decoder.covariance_check(*first_values, sigma),
        "runtime_seconds": time.perf_counter() - started,
        "reference_values_or_labels_accessed": False,
        "probe_points_source": "frozen GLIDER label-free prediction files",
    }
    args.audit_output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    (args.output / "manifest.json").write_text(json.dumps({"method": args.method, "prediction_registry_sha256": sha256(registry), "audit_sha256": sha256(args.audit_output), "n_configurations": len(rows)}, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
