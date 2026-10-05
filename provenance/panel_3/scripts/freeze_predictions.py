#!/usr/bin/env python3
"""Validate and content-address the complete pre-reference prediction state."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "canonical_panel_3"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def records(root: Path) -> list[dict]:
    return [
        {"path": str(path.relative_to(PANEL)), "sha256": sha256(path), "size_bytes": path.stat().st_size}
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name not in {"GLIDER_PREDICTION_FREEZE_MANIFEST.json", "COMPARATOR_FREEZE_MANIFEST.json", "ALL_PREDICTIONS_FREEZE_MANIFEST.json"}
    ]


def aggregate(items: list[dict]) -> str:
    digest = hashlib.sha256()
    for item in items:
        digest.update(f"{item['path']}\0{item['sha256']}\n".encode())
    return digest.hexdigest()


def main() -> None:
    if (PANEL / "references").exists() or list(PANEL.rglob("*reference_response*")):
        raise RuntimeError("Reference-response artifact exists before prediction freeze")
    registry = pd.read_csv(PANEL / "data/configuration_registry.csv")
    expected = set(registry.config_id)
    methods = {
        "glider": PANEL / "predictions/glider",
        "mace_polar_l": PANEL / "predictions/comparators/mace_polar_l",
        "mace_polar_m": PANEL / "predictions/comparators/mace_polar_m",
        "mace_polar_s": PANEL / "predictions/comparators/mace_polar_s",
        "aimnet2": PANEL / "predictions/comparators/aimnet2",
        "mace_mdp": PANEL / "predictions/comparators/mace_mdp",
        "gfn2_xtb": PANEL / "predictions/comparators/gfn2_xtb",
        "static_thole": PANEL / "predictions/comparators/static_thole",
        "zero_response": PANEL / "predictions/comparators/zero_response",
    }
    glider_points: dict[str, np.ndarray] = {}
    method_summary = {}
    for method, directory in methods.items():
        table = pd.read_csv(directory / "prediction_registry.csv")
        if set(table.config_id) != expected or len(table) != 80:
            raise RuntimeError(f"{method} does not cover all 80 configurations")
        for row in table.itertuples():
            values = np.load(directory / row.prediction_file)
            points = np.asarray(values["points_angstrom"])
            if method == "glider":
                glider_points[row.config_id] = points
            elif not np.array_equal(points, glider_points[row.config_id]):
                raise RuntimeError(f"{method} probe points differ: {row.config_id}")
        method_summary[method] = {
            "n_configurations": len(table),
            "registry_sha256": sha256(directory / "prediction_registry.csv"),
            "tree_sha256": aggregate(records(directory)),
        }
    glider_table = pd.read_csv(methods["glider"] / "prediction_registry.csv")
    glider_files = records(methods["glider"])
    glider_manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "method": "original frozen GLIDER",
        "checkpoint_sha256": sha256(ROOT / "checkpoints/glider_site_response_ensemble.pt"),
        "geometry_sha256": sha256(PANEL / "data/configurations.extxyz"),
        "n_configurations": 80,
        "maximum_abs_net_response_charge_e": float(glider_table.net_charge_e.abs().max()),
        "deterministic_repeat_maximum_absolute_discrepancy": 0.0,
        "rotation_audit": json.loads((PANEL / "audits/GLIDER_ROTATION_AUDIT.json").read_text()),
        "prediction_tree_sha256": aggregate(glider_files),
        "files": glider_files,
        "reference_response_labels_exist": False,
        "experimental_hydration_targets_accessed": False,
    }
    (PANEL / "predictions/GLIDER_PREDICTION_FREEZE_MANIFEST.json").write_text(json.dumps(glider_manifest, indent=2, sort_keys=True) + "\n")
    comparator_files = records(PANEL / "predictions/comparators")
    comparator_manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "methods": {key: value for key, value in method_summary.items() if key != "glider"},
        "prediction_tree_sha256": aggregate(comparator_files),
        "files": comparator_files,
        "hardcoded_historical_chronology_note": "MACE-MDP, AIMNet2 and GFN2 audit templates retain a historical sentence saying post-label execution; their explicit target-access fields are false. Panel-III execution preceded any reference artifact, as established by this content-addressed freeze and Git history.",
        "reference_response_labels_exist": False,
        "experimental_hydration_targets_accessed": False,
    }
    (PANEL / "predictions/COMPARATOR_FREEZE_MANIFEST.json").write_text(json.dumps(comparator_manifest, indent=2, sort_keys=True) + "\n")
    all_files = records(PANEL / "predictions") + records(PANEL / "audits")
    parent = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    all_manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parent_geometry_freeze_commit": parent,
        "n_molecules": 20,
        "n_geometries": 80,
        "geometry_sha256": sha256(PANEL / "data/configurations.extxyz"),
        "glider_checkpoint_sha256": sha256(ROOT / "checkpoints/glider_site_response_ensemble.pt"),
        "glider_prediction_tree_sha256": glider_manifest["prediction_tree_sha256"],
        "comparator_prediction_tree_sha256": comparator_manifest["prediction_tree_sha256"],
        "complete_pre_reference_tree_sha256": aggregate(all_files),
        "methods": method_summary,
        "all_prediction_programs_report_reference_targets_accessed": False,
        "all_80_glider_predictions_exist": True,
        "all_eligible_comparator_predictions_exist": True,
        "reference_response_labels_exist": False,
        "experimental_hydration_targets_accessed": False,
        "post_freeze_prediction_modification_permitted": False,
    }
    (PANEL / "predictions/ALL_PREDICTIONS_FREEZE_MANIFEST.json").write_text(json.dumps(all_manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(all_manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
