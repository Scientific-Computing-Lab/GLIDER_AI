#!/usr/bin/env python3
"""Freeze the label-blind non-water contact panel before new QM acquisition."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from ase.io import read


ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT = ROOT / "build/nonwater_contact_panel"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def prediction_digest(path: Path, case_ids: set[str]) -> dict[str, object]:
    registry = path / "prediction_registry.csv"
    with registry.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 36 or {row["config_id"] for row in rows} != case_ids:
        raise ValueError(f"Prediction coverage mismatch in {path}")
    for row in rows:
        prediction = path / row["prediction_file"]
        if sha256(prediction) != row["prediction_sha256"]:
            raise ValueError(f"Prediction digest mismatch: {prediction}")
    return {
        "registry_sha256": sha256(registry),
        "case_count": len(rows),
        "case_file_sha256": {
            row["config_id"]: row["prediction_sha256"] for row in rows
        },
    }


def main() -> None:
    freeze_path = EXPERIMENT / "prediction_freeze.json"
    if freeze_path.exists():
        raise RuntimeError("A prediction freeze already exists; do not overwrite it")
    if (EXPERIMENT / "references").exists():
        raise RuntimeError("Reference directory already exists; freeze must precede QM")
    geometry = EXPERIMENT / "configurations.extxyz"
    frames = read(geometry, index=":")
    case_ids = {str(atoms.info["config_id"]) for atoms in frames}
    if len(frames) != 36 or len(case_ids) != 36:
        raise ValueError("Expected exactly 36 unique selected contact geometries")
    design = json.loads((EXPERIMENT / "design_manifest.json").read_text())
    if design["geometry_sha256"] != sha256(geometry):
        raise ValueError("Geometry changed since label-blind design")
    glider = EXPERIMENT / "predictions/glider"
    baseline = EXPERIMENT / "predictions/mace_polar_l"
    with (glider / "manifest.json").open() as stream:
        glider_manifest = json.load(stream)
    with (baseline / "manifest.json").open() as stream:
        baseline_manifest = json.load(stream)
    if baseline_manifest.get("reference_labels_accessed") is not False:
        raise ValueError("MACE-L provenance does not attest label blindness")
    record = {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "36 geometry selections and both prediction sets frozen before new QM labels",
        "geometry_sha256": sha256(geometry),
        "design_manifest_sha256": sha256(EXPERIMENT / "design_manifest.json"),
        "glider_checkpoint_sha256": glider_manifest["checkpoint_sha256"],
        "mace_polar_l_checkpoint_sha256": json.loads(
            (EXPERIMENT / "mace_l_raw.manifest.json").read_text()
        )["checkpoint_sha256"],
        "glider_predictions": prediction_digest(glider, case_ids),
        "mace_polar_l_predictions": prediction_digest(baseline, case_ids),
        "new_qm_reference_labels_accessed": False,
    }
    freeze_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(f"Frozen {len(frames)} cases: {sha256(freeze_path)}")


if __name__ == "__main__":
    main()
