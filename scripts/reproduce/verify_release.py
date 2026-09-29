#!/usr/bin/env python3
"""Lightweight release verification: hashes, registries, physics and metrics."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from glider.data import validate_registry
from glider.inference import dipole_from_sites

ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_predictions(panel: str, method: str) -> tuple[float, float]:
    base = ROOT / (f"data/{panel}" if panel != "canonical_panel_3" else "canonical_panel_3")
    registry = pd.read_csv(base / "data/configuration_registry.csv" if panel == "canonical_panel_3" else base / "configuration_registry.csv")
    refs = pd.read_csv(base / "references/observable_registry.csv")
    max_charge = 0.0
    max_dipole = 0.0
    for row in registry.itertuples(index=False):
        config_id = str(row.config_id)
        # Observable filenames are content identifiers; resolve through the reference registry.
        filename = refs.loc[refs.config_id == config_id, "observable_file"].iloc[0]
        prediction = np.load(base / f"predictions/{method}" / filename)
        max_charge = max(max_charge, abs(float(prediction["predicted_charges_e"].sum())))
        rebuilt = dipole_from_sites(
            np.asarray(read_positions(panel, config_id)),
            prediction["predicted_charges_e"],
            prediction["predicted_dipoles_e_bohr"],
        )
        max_dipole = max(
            max_dipole, float(np.max(np.abs(rebuilt - prediction["predicted_dipole_debye"])))
        )
    return max_charge, max_dipole


def read_positions(panel: str, config_id: str) -> np.ndarray:
    """Minimal extxyz reader sufficient for immutable release geometries."""
    path = ROOT / (f"data/{panel}/configurations.extxyz" if panel != "canonical_panel_3" else "canonical_panel_3/data/configurations.extxyz")
    with path.open() as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            n_atoms = int(line.strip())
            comment = handle.readline()
            rows = [handle.readline().split() for _ in range(n_atoms)]
            if f"config_id={config_id}" in comment or f'config_id="{config_id}"' in comment:
                return np.asarray([[float(x[1]), float(x[2]), float(x[3])] for x in rows])
    raise KeyError(config_id)


def main() -> None:
    checkpoint_manifest = json.loads(
        (ROOT / "checkpoints/glider_site_response_ensemble.manifest.json").read_text()
    )
    observed = digest(ROOT / "checkpoints/glider_site_response_ensemble.pt")
    assert observed == checkpoint_manifest["checkpoint_sha256"], (
        observed,
        checkpoint_manifest["checkpoint_sha256"],
    )

    freeze_manifest = json.loads(
        (ROOT / "evidence/manifests/prospective_freeze_manifest.json").read_text()
    )
    development_root = ROOT / "results/development/original"
    for relative, expected in freeze_manifest["development_artifact_hashes"].items():
        copied = development_root / relative
        if copied.exists():
            assert digest(copied) == expected, relative
    validate_registry(
        ROOT / "data/prospective_1/configuration_registry.csv", molecules=12, configurations=48
    )
    validate_registry(
        ROOT / "data/prospective_2/configuration_registry.csv", molecules=24, configurations=96
    )
    validate_registry(
        ROOT / "canonical_panel_3/data/configuration_registry.csv", molecules=20, configurations=80
    )
    p3_freeze = json.loads((ROOT / "canonical_panel_3/predictions/ALL_PREDICTIONS_FREEZE_MANIFEST.json").read_text())
    assert p3_freeze["glider_checkpoint_sha256"] == observed
    assert p3_freeze["reference_response_labels_exist"] is False
    p3_ref = json.loads((ROOT / "canonical_panel_3/references/REFERENCE_MANIFEST.json").read_text())
    assert p3_ref["reference_audit_pass"] is True
    assert p3_ref["n_complete_references"] == 80 and p3_ref["n_excluded"] == 0
    source = json.loads((ROOT / "canonical_panel_3/source/SOURCE_MANIFEST.json").read_text())
    assert source["hydration_targets_accessed"] is False
    assert source["fields_after_index_1_parsed_or_stored"] is False
    for panel, method in [
        ("prospective_1", "glider"),
        ("prospective_2", "previous_candidate"),
        ("canonical_panel_3", "glider"),
    ]:
        charge, dipole = validate_predictions(panel, method)
        assert charge < 1e-5, charge
        assert dipole < 1e-6, dipole
        print(f"{panel}: max |sum dq|={charge:.3e} e; max dipole reconstruction={dipole:.3e} D")
    print("checkpoint and registry verification: PASS")


if __name__ == "__main__":
    main()
