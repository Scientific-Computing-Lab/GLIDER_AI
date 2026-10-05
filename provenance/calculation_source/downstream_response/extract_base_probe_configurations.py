#!/usr/bin/env python3
"""Materialize the fixed liquid-bridge base clusters and held-out fourth waters.

This script performs geometry selection only.  It never reads a prediction or a
QM response value.  The parent snapshots are the complete, pre-existing
6-solute x 4-snapshot liquid bridge, and the ordering is the original
minimum-image oxygen-to-nearest-solute-atom rule with source water index as the
deterministic tie-break.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[3]
WORKSPACE = ROOT.parent
SOURCE_REGISTRY = ROOT / "evidence/liquid_bridge/source/LIQUID_CLUSTER_REGISTRY.csv"
SOURCE_GEOMETRIES = ROOT / "evidence/liquid_bridge/source/liquid_3water.extxyz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    fractional = np.linalg.solve(box.T, np.asarray(delta).T).T
    fractional -= np.rint(fractional)
    return fractional @ box


def topology_paths(trajectory: Path, molecule_id: str) -> tuple[Path, Path]:
    if "hamiltonian_lift/data/endpoints" in str(trajectory):
        directory = WORKSPACE / "hamiltonian_lift/data/h0" / molecule_id
    else:
        directory = trajectory.parent
    return directory / "initial_solvated.pdb", directory / "system_manifest.json"


def ordered_waters(
    positions_angstrom: np.ndarray,
    box_angstrom: np.ndarray,
    n_solute: int,
) -> tuple[np.ndarray, list[tuple[float, int, np.ndarray]]]:
    solute = positions_angstrom[:n_solute].copy()
    solvent = positions_angstrom[n_solute:]
    if len(solvent) % 3:
        raise RuntimeError("expected solute followed by O-H-H waters")
    waters = solvent.reshape(-1, 3, 3)
    ranked: list[tuple[float, int, np.ndarray]] = []
    for source_index, water in enumerate(waters):
        displacement = minimum_image(water[0] - solute, box_angstrom)
        nearest_solute = int(np.argmin(np.linalg.norm(displacement, axis=1)))
        oxygen = solute[nearest_solute] + displacement[nearest_solute]
        hydrogen_1 = oxygen + minimum_image(water[1] - water[0], box_angstrom)
        hydrogen_2 = oxygen + minimum_image(water[2] - water[0], box_angstrom)
        distance = float(np.linalg.norm(oxygen - solute[nearest_solute]))
        ranked.append(
            (distance, source_index, np.vstack((oxygen, hydrogen_1, hydrogen_2)))
        )
    return solute, sorted(ranked, key=lambda item: (item[0], item[1]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results/posthoc_downstream_response/configurations",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    registry = pd.read_csv(SOURCE_REGISTRY)
    registry = registry.loc[registry["shell_waters"] == 3].copy()
    if len(registry) != 24 or registry["molecule_id"].nunique() != 6:
        raise RuntimeError("the authoritative liquid bridge is not 6 solutes x 4 snapshots")
    source_frames = {
        str(atoms.info["config_id"]): atoms
        for atoms in read(SOURCE_GEOMETRIES, index=":")
    }
    if set(source_frames) != set(registry["config_id"]):
        raise RuntimeError("registry and frozen 3-water geometry identities differ")

    bases: list[Atoms] = []
    probes: list[Atoms] = []
    combined: list[Atoms] = []
    rows: list[dict[str, object]] = []
    maximum_base_discrepancy = 0.0

    for row in registry.sort_values(["molecule_id", "source_frame_index"]).itertuples():
        source_trajectory = Path(row.source_trajectory)
        topology, system_manifest = topology_paths(source_trajectory, row.molecule_id)
        if sha256(source_trajectory) != row.source_trajectory_sha256:
            raise RuntimeError(f"trajectory hash mismatch: {source_trajectory}")
        manifest = json.loads(system_manifest.read_text())
        n_solute = int(manifest["n_solute_atoms"])
        if n_solute != int(row.n_solute_atoms):
            raise RuntimeError(f"solute atom count mismatch: {row.config_id}")
        values = np.load(source_trajectory, allow_pickle=False)
        positions = np.asarray(values["positions_nm"][row.source_frame_index]) * 10.0
        box = np.asarray(values["box_vectors_nm"][row.source_frame_index]) * 10.0
        solute, ranking = ordered_waters(positions, box, n_solute)
        selected = ranking[:4]

        reconstructed_base_positions = np.vstack(
            (solute, selected[0][2], selected[1][2], selected[2][2])
        )
        source_base = source_frames[row.config_id].copy()
        discrepancy = float(
            np.max(np.abs(reconstructed_base_positions - source_base.positions))
        )
        maximum_base_discrepancy = max(maximum_base_discrepancy, discrepancy)
        if discrepancy > 1.0e-7:
            raise RuntimeError(
                f"reconstructed base differs from frozen geometry by {discrepancy}: "
                f"{row.config_id}"
            )

        case_id = str(row.config_id).replace("__3w", "__base3_probe4")
        water_indices = [int(item[1]) for item in selected]
        water_distances = [float(item[0]) for item in selected]
        base = source_base.copy()
        base.info.update(
            {
                "case_id": case_id,
                "base_config_id": str(row.config_id),
                "heldout_water_rank": 4,
                "heldout_source_water_index": water_indices[3],
            }
        )
        probe = Atoms(("O", "H", "H"), positions=selected[3][2])
        probe.info.update(
            {
                "case_id": case_id,
                "base_config_id": str(row.config_id),
                "molecule_id": str(row.molecule_id),
                "source_frame_index": int(row.source_frame_index),
                "heldout_water_rank": 4,
                "heldout_source_water_index": water_indices[3],
                "heldout_oxygen_distance_A": water_distances[3],
            }
        )
        joined = base + probe
        joined.info = dict(base.info)
        joined.info.update(
            {
                "n_base_atoms": len(base),
                "n_probe_atoms": 3,
                "n_solute_atoms": n_solute,
                "source_water_indices_rank1_to_4": ";".join(map(str, water_indices)),
                "water_distances_A_rank1_to_4": ";".join(
                    f"{value:.10f}" for value in water_distances
                ),
            }
        )
        bases.append(base)
        probes.append(probe)
        combined.append(joined)
        rows.append(
            {
                "case_id": case_id,
                "base_config_id": str(row.config_id),
                "molecule_id": str(row.molecule_id),
                "source_frame_index": int(row.source_frame_index),
                "n_solute_atoms": n_solute,
                "n_base_atoms": len(base),
                "w1_source_index": water_indices[0],
                "w2_source_index": water_indices[1],
                "w3_source_index": water_indices[2],
                "w4_source_index": water_indices[3],
                "w1_distance_A": water_distances[0],
                "w2_distance_A": water_distances[1],
                "w3_distance_A": water_distances[2],
                "w4_distance_A": water_distances[3],
                "source_trajectory": str(source_trajectory),
                "source_trajectory_sha256": sha256(source_trajectory),
                "source_topology": str(topology),
                "source_topology_sha256": sha256(topology),
                "base_geometry_source": str(SOURCE_GEOMETRIES.relative_to(ROOT)),
                "base_geometry_reconstruction_max_abs_A": discrepancy,
                "selection_used_model_or_qm_performance": False,
            }
        )

    base_path = args.output / "base_3water.extxyz"
    probe_path = args.output / "heldout_w4.extxyz"
    combined_path = args.output / "base_plus_probe.extxyz"
    registry_path = args.output / "configurations.csv"
    write(base_path, bases, format="extxyz")
    write(probe_path, probes, format="extxyz")
    write(combined_path, combined, format="extxyz")
    pd.DataFrame(rows).to_csv(registry_path, index=False)

    output_manifest = {
        "experiment": "post hoc downstream response coupling to held-out fourth water",
        "n_cases": len(rows),
        "n_solutes": len({row["molecule_id"] for row in rows}),
        "n_snapshots_per_solute": 4,
        "base_definition": "solute plus nearest waters ranked 1, 2, and 3",
        "probe_definition": "held-out water ranked 4 from the same parent snapshot",
        "ordering_rule": (
            "minimum-image oxygen distance to nearest solute atom; "
            "source water index tie-break"
        ),
        "maximum_base_reconstruction_discrepancy_A": maximum_base_discrepancy,
        "source_registry": str(SOURCE_REGISTRY.relative_to(ROOT)),
        "source_registry_sha256": sha256(SOURCE_REGISTRY),
        "source_geometries": str(SOURCE_GEOMETRIES.relative_to(ROOT)),
        "source_geometries_sha256": sha256(SOURCE_GEOMETRIES),
        "outputs": {
            path.name: sha256(path)
            for path in (base_path, probe_path, combined_path, registry_path)
        },
        "performance_based_selection": False,
        "prospective_claim": False,
    }
    (args.output / "extraction_manifest.json").write_text(
        json.dumps(output_manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(output_manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
