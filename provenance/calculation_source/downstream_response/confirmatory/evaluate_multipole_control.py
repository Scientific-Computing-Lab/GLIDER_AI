#!/usr/bin/env python3
"""Evaluate the frozen exact-QM global-multipole hierarchy control."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from multipole_control import (
    CONTROL_ROOT,
    HIERARCHIES,
    ORIGINS,
    MultipoleData,
    case_origins,
    hierarchy_couplings,
    source_extent_for_case,
)


class CachedEvaluator:
    def __init__(self, data: MultipoleData, output: Path) -> None:
        self.data = data
        self.cache = output / "raw_evaluations"
        self.cache.mkdir(parents=True, exist_ok=True)

    def evaluate(
        self,
        probe_id: str,
        evaluation_id: str,
        positions: np.ndarray,
        density: np.ndarray,
    ) -> list[dict[str, object]]:
        specification = f"{probe_id}::{evaluation_id}::h0.01::seven-point"
        tag = hashlib.sha256(specification.encode()).hexdigest()[:24]
        path = self.cache / f"{tag}.json"
        if path.exists():
            return json.loads(path.read_text())
        case_id = str(self.data.probes[probe_id].info["case_id"])
        values = hierarchy_couplings(
            base=self.data.bases[case_id],
            density_values=self.data.base_values(case_id),
            probe_positions_angstrom=positions,
            probe_density=density,
            step_bohr=0.01,
        )
        rows = [
            {
                **self.data.metadata(probe_id),
                "evaluation_id": evaluation_id,
                "hierarchy": hierarchy,
                "origin": origin,
                **result,
            }
            for (hierarchy, origin), result in values.items()
        ]
        path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
        return rows


def run_initial(data: MultipoleData, evaluator: CachedEvaluator, output: Path) -> None:
    rows = []
    for index, probe_id in enumerate(sorted(data.probes), 1):
        values = data.probe_values(probe_id)
        rows.extend(
            evaluator.evaluate(
                probe_id,
                "initial",
                values["probe_positions_angstrom"],
                values["probe_density_matrix"],
            )
        )
        print(f"multipole initial [{index}/90] {probe_id}", flush=True)
    table = pd.DataFrame(rows).sort_values(
        ["molecule_id", "water_rank", "hierarchy", "origin"]
    )
    table.to_csv(output / "W4_W12_energy_all_origins.csv", index=False)
    table[table.water_rank == 4].to_csv(
        output / "W4_energy_all_origins.csv", index=False
    )


def run_torque(data: MultipoleData, evaluator: CachedEvaluator, output: Path) -> None:
    existing = pd.read_csv(
        data.base_density_files[next(iter(data.base_density_files))].parents[2]
        / "evaluated/breadth_torque_per_solute.csv"
    )
    rows = []
    for index, probe_id in enumerate(data.w4_ids(), 1):
        reference = existing[
            (existing.probe_id == probe_id) & (existing.method == "qm")
        ][["x", "y", "z"]].to_numpy()[0]
        energy = {
            (hierarchy, origin): np.empty((3, 2))
            for hierarchy in HIERARCHIES
            for origin in ORIGINS
        }
        for axis in range(3):
            for sign_index, label in enumerate(("minus", "plus")):
                variant_id = f"torque_0.1deg_axis{axis}_{label}"
                variant = data.variant_values(probe_id, variant_id)
                evaluated = evaluator.evaluate(
                    probe_id,
                    variant_id,
                    variant["positions_angstrom"],
                    variant["density_matrix"],
                )
                for row in evaluated:
                    energy[(row["hierarchy"], row["origin"])][
                        axis, sign_index
                    ] = row["energy_kcal_mol"]
        for (hierarchy, origin), values in energy.items():
            vector = -(values[:, 1] - values[:, 0]) / (2.0 * np.radians(0.1))
            error = vector - reference
            reference_norm = float(np.linalg.norm(reference))
            vector_norm = float(np.linalg.norm(vector))
            if reference_norm > 0 and vector_norm > 0:
                cosine = float(
                    np.clip(
                        np.dot(vector, reference) / (vector_norm * reference_norm),
                        -1.0,
                        1.0,
                    )
                )
                angle = float(np.degrees(np.arccos(cosine)))
            else:
                cosine = float("nan")
                angle = float("nan")
            rows.append(
                {
                    **data.metadata(probe_id),
                    "hierarchy": hierarchy,
                    "origin": origin,
                    "x": vector[0],
                    "y": vector[1],
                    "z": vector[2],
                    "magnitude": vector_norm,
                    "reference_x": reference[0],
                    "reference_y": reference[1],
                    "reference_z": reference[2],
                    "reference_magnitude": reference_norm,
                    "vector_error_norm": float(np.linalg.norm(error)),
                    "magnitude_error": vector_norm - reference_norm,
                    "cosine_similarity": cosine,
                    "angular_error_degrees": angle,
                }
            )
        print(f"multipole torque [{index}/10] {probe_id}", flush=True)
    pd.DataFrame(rows).sort_values(
        ["molecule_id", "hierarchy", "origin"]
    ).to_csv(output / "W4_torque_all_origins.csv", index=False)


def run_orientation(
    data: MultipoleData, evaluator: CachedEvaluator, output: Path
) -> None:
    existing = pd.read_csv(
        data.base_density_files[next(iter(data.base_density_files))].parents[2]
        / "evaluated/orientation_profiles.csv"
    )
    reference = existing[existing.method == "qm"][
        [
            "probe_id",
            "rotation_id",
            "rotation_index",
            "quaternion_x",
            "quaternion_y",
            "quaternion_z",
            "quaternion_w",
            "energy_kcal_mol",
        ]
    ].rename(columns={"energy_kcal_mol": "reference_energy_kcal_mol"})
    rows = []
    for case_index, probe_id in enumerate(data.w4_ids(), 1):
        for rotation_index in range(24):
            rotation_id = f"O_{rotation_index:02d}"
            variant_id = f"orientation_{rotation_id}"
            variant = data.variant_values(probe_id, variant_id)
            evaluated = evaluator.evaluate(
                probe_id,
                variant_id,
                variant["positions_angstrom"],
                variant["density_matrix"],
            )
            reference_row = reference[
                (reference.probe_id == probe_id)
                & (reference.rotation_id == rotation_id)
            ].iloc[0]
            for row in evaluated:
                rows.append(
                    {
                        **row,
                        "rotation_id": rotation_id,
                        "rotation_index": rotation_index,
                        "quaternion_x": reference_row.quaternion_x,
                        "quaternion_y": reference_row.quaternion_y,
                        "quaternion_z": reference_row.quaternion_z,
                        "quaternion_w": reference_row.quaternion_w,
                        "reference_energy_kcal_mol": reference_row.reference_energy_kcal_mol,
                    }
                )
        print(f"multipole orientation [{case_index}/10] {probe_id}", flush=True)
    pd.DataFrame(rows).sort_values(
        ["molecule_id", "hierarchy", "origin", "rotation_index"]
    ).to_csv(output / "W4_orientation_profiles_all_origins.csv", index=False)


def run_source_extent(data: MultipoleData, output: Path) -> None:
    rows = []
    diagnostics = []
    probes_by_case: dict[str, list[str]] = {}
    for probe_id, probe in data.probes.items():
        probes_by_case.setdefault(str(probe.info["case_id"]), []).append(probe_id)
    for index, (case_id, base) in enumerate(sorted(data.bases.items()), 1):
        extents, case_diagnostic = source_extent_for_case(
            base, data.base_values(case_id)
        )
        diagnostics.append(
            {
                "case_id": case_id,
                "molecule_id": str(base.info["molecule_id"]),
                "family": str(base.info["family"]),
                **case_diagnostic,
            }
        )
        origins = case_origins(base)
        for probe_id in sorted(probes_by_case[case_id]):
            probe_values = data.probe_values(probe_id)
            positions = probe_values["probe_positions_angstrom"]
            for origin_name, origin in origins.items():
                oxygen_distance = float(np.linalg.norm(positions[0] - origin))
                minimum_nuclear = float(
                    np.min(np.linalg.norm(positions - origin, axis=1))
                )
                extent = extents[origin_name]
                rows.append(
                    {
                        **data.metadata(probe_id),
                        "origin": origin_name,
                        **extent,
                        "probe_oxygen_to_origin_angstrom": oxygen_distance,
                        "probe_min_nuclear_to_origin_angstrom": minimum_nuclear,
                        "oxygen_over_r90": oxygen_distance
                        / extent["r90_angstrom"],
                        "oxygen_over_r95": oxygen_distance
                        / extent["r95_angstrom"],
                        "oxygen_over_r99": oxygen_distance
                        / extent["r99_angstrom"],
                        "oxygen_over_max_base_nuclear_radius": oxygen_distance
                        / extent["max_base_nuclear_radius_angstrom"],
                        "min_nuclear_over_r99": minimum_nuclear
                        / extent["r99_angstrom"],
                    }
                )
        print(f"source extent [{index}/10] {case_id}", flush=True)
    pd.DataFrame(rows).sort_values(
        ["molecule_id", "water_rank", "origin"]
    ).to_csv(output / "source_extent_geometry.csv", index=False)
    pd.DataFrame(diagnostics).sort_values("molecule_id").to_csv(
        output / "source_extent_grid_diagnostics.csv", index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=CONTROL_ROOT)
    parser.add_argument(
        "--stage",
        choices=("initial", "torque", "orientation", "extent", "all"),
        default="all",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = MultipoleData()
    evaluator = CachedEvaluator(data, args.output)
    if args.stage in ("initial", "all"):
        run_initial(data, evaluator, args.output)
    if args.stage in ("torque", "all"):
        run_torque(data, evaluator, args.output)
    if args.stage in ("orientation", "all"):
        run_orientation(data, evaluator, args.output)
    if args.stage in ("extent", "all"):
        run_source_extent(data, args.output)


if __name__ == "__main__":
    main()
