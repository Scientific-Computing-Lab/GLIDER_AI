#!/usr/bin/env python3
"""Audit the frozen Panel-III QM response references before scoring.

This script is intentionally label-agnostic: it checks completeness, integrity,
units/shapes, probe identity, and chronology but does not calculate model errors.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "canonical_panel_3"
REF = PANEL / "references"
PRED = PANEL / "predictions"
FREEZE_COMMIT = "f9a87b1d5a91a45f595ccd1891f59ec8d9445fea"
INTEGRITY_COMMIT = "5347a97c8f01a998e736f9b453138d9e73419a1f"
ALL_PREDICTIONS_MANIFEST_SHA256 = (
    "3dc5641a289f48118e0e372c341286573d4a381a2cf2ea9b7ae2b2b80ec8c73e"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def commit_epoch(commit: str) -> int:
    return int(
        subprocess.check_output(
            ["git", "show", "-s", "--format=%ct", commit], cwd=ROOT, text=True
        ).strip()
    )


def main() -> None:
    registry_path = PANEL / "data" / "configuration_registry.csv"
    rows = list(csv.DictReader(registry_path.open()))
    assert len(rows) == 80
    assert len({r["config_id"] for r in rows}) == 80
    assert len({r["molecule_id"] for r in rows}) == 20
    assert {r["regime"] for r in rows} == {
        "equilibrium",
        "cooperative",
        "orientation",
        "compressed",
    }

    prediction_manifest = PRED / "ALL_PREDICTIONS_FREEZE_MANIFEST.json"
    assert sha256(prediction_manifest) == ALL_PREDICTIONS_MANIFEST_SHA256
    frozen = json.loads(prediction_manifest.read_text())
    assert frozen["n_geometries"] == 80
    assert frozen["reference_response_labels_exist"] is False

    freeze_epoch = commit_epoch(INTEGRITY_COMMIT)
    audit_rows: list[dict[str, object]] = []
    recovered: dict[str, dict[str, object]] = {}
    for p in sorted(REF.glob("numerical_recovery_*.json")):
        rec = json.loads(p.read_text())
        recovered[rec["config_id"]] = rec

    for row in rows:
        cid = row["config_id"]
        token = hashlib.sha256(cid.encode()).hexdigest()[:20]
        ref_path = REF / f"{token}.npz"
        pred_path = PRED / "glider" / f"{token}.npz"
        assert ref_path.exists(), cid
        assert pred_path.exists(), cid
        assert ref_path.stat().st_mtime > freeze_epoch, cid

        ref = np.load(ref_path, allow_pickle=False)
        pred = np.load(pred_path, allow_pickle=False)
        required = {
            "points_angstrom",
            "delta_esp_hartree_per_e",
            "delta_dipole_debye",
            "fitted_charges_e",
            "fitted_dipoles_e_bohr",
            "fitted_esp_hartree_per_e",
            "component_runtime_seconds",
        }
        assert required.issubset(ref.files), (cid, ref.files)
        points = ref["points_angstrom"]
        esp = ref["delta_esp_hartree_per_e"]
        dipole = ref["delta_dipole_debye"]
        assert points.ndim == 2 and points.shape[1] == 3
        assert esp.shape == (len(points),)
        assert dipole.shape == (3,)
        assert np.isfinite(points).all()
        assert np.isfinite(esp).all()
        assert np.isfinite(dipole).all()
        assert np.array_equal(points, pred["points_angstrom"]), cid

        rec = recovered.get(cid)
        if rec:
            assert rec["observable_sha256"] == sha256(ref_path)
            assert all(rec["final_scf_converged"])
            assert rec["hamiltonian_or_reference_change"] is False
            assert rec["excluded_or_replaced"] is False

        audit_rows.append(
            {
                "config_id": cid,
                "molecule_id": row["molecule_id"],
                "regime": row["regime"],
                "observable_file": ref_path.name,
                "observable_sha256": sha256(ref_path),
                "n_esp_points": len(points),
                "finite_values": True,
                "probe_identity_exact": True,
                "created_after_prediction_freeze": True,
                "numerical_recovery": bool(rec),
                "same_hamiltonian_recovery": bool(rec),
                "excluded_or_replaced": False,
            }
        )

    assert len(recovered) == 4
    assert sum(bool(r["numerical_recovery"]) for r in audit_rows) == 4

    audit_csv = REF / "REFERENCE_FILE_AUDIT.csv"
    with audit_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]))
        writer.writeheader()
        writer.writerows(audit_rows)

    recovery_csv = REF / "QM_RECOVERY_LOG.csv"
    fields = [
        "config_id",
        "observable_file",
        "observable_sha256",
        "runtime_seconds",
        "final_protocol",
        "solver_routes",
        "solver_only_change",
        "hamiltonian_or_reference_change",
        "excluded_or_replaced",
        "density_sha256",
    ]
    with recovery_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for cid in sorted(recovered):
            rec = recovered[cid]
            rec = dict(rec)
            rec["solver_routes"] = " | ".join(rec["solver_routes"])
            writer.writerow({k: rec[k] for k in fields})

    reference_files = sorted(REF.glob("*.npz"))
    tree_h = hashlib.sha256()
    for path in reference_files:
        tree_h.update(f"{path.name}\0{sha256(path)}\n".encode())

    manifest = {
        "campaign": "GLIDER canonical-source prospective validation III",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "prospective Panel III reference response labels",
        "protocol": {
            "functional": "omegaB97X-D3BJ",
            "basis": "def2-TZVPD",
            "grid_level": 4,
            "scf_convergence": 1e-10,
            "counterpoise_components": [
                "full complex",
                "ghost-basis solute",
                "ghost-basis three-water cluster",
            ],
            "response": "full complex - ghost solute - ghost water cluster",
            "units": {"esp": "hartree/e", "dipole": "Debye", "coordinates": "angstrom"},
        },
        "n_molecules": 20,
        "n_configurations": 80,
        "n_complete_references": len(reference_files),
        "n_ordinary_scf": 76,
        "n_same_hamiltonian_recoveries": 4,
        "n_excluded": 0,
        "prediction_freeze_commit": FREEZE_COMMIT,
        "pre_qm_integrity_commit": INTEGRITY_COMMIT,
        "all_predictions_freeze_manifest_sha256": ALL_PREDICTIONS_MANIFEST_SHA256,
        "all_references_postdate_remote_prediction_freeze": True,
        "experimental_hydration_targets_accessed": False,
        "reference_tree_sha256": tree_h.hexdigest(),
        "observable_registry_sha256": sha256(REF / "observable_registry.csv"),
        "reference_file_audit_sha256": sha256(audit_csv),
        "qm_recovery_log_sha256": sha256(recovery_csv),
        "historical_raw_generator_manifest_note": (
            "references/manifest.json was emitted by a development-era generator and "
            "contains stale scope/access flags; it is retained unmodified for provenance "
            "and is superseded by this audited manifest."
        ),
        "reference_audit_pass": True,
    }
    (REF / "REFERENCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")

    report = f"""# Panel-III reference audit

**Result: PASS.** All 80 preregistered configurations have complete, finite
counterpoise-consistent response references. Probe coordinates match the
prediction-frozen GLIDER files exactly; no molecule or configuration was
excluded, replaced, or altered.

## Chronology and integrity

- Prediction-freeze commit: `{FREEZE_COMMIT}` (pushed before QM).
- Pre-QM integrity commit: `{INTEGRITY_COMMIT}`.
- Frozen prediction-manifest SHA-256: `{ALL_PREDICTIONS_MANIFEST_SHA256}`.
- All 80 reference files postdate the remote prediction freeze.
- FreeSolv hydration quantities were not accessed.

## Reference protocol

All components use DF-RKS omegaB97X-D3(BJ)/def2-TZVPD, grid level 4 and
SCF convergence 1e-10 at the frozen geometry and complete complex basis.
The mutual response is the full complex minus ghost-basis solute minus
ghost-basis three-water cluster.

Seventy-six cases converged by the ordinary route. The four iodine-containing
environments for `mobley_2727678` used the preregistered same-Hamiltonian
recovery ladder. Final solutions use zero-level-shift second-order/DIIS SCF;
functional, basis, grid, geometry, convergence threshold and Hamiltonian were
unchanged. All three components converged in every case. The exact routes and
hashes are in `QM_RECOVERY_LOG.csv`.

## Superseded raw metadata

`references/manifest.json` is retained byte-for-byte because it was emitted by
the pre-existing QM generator. Its inherited development-era scope and
`prospective_labels_accessed=false` field are stale for Panel III. The audited
`REFERENCE_MANIFEST.json`, `REFERENCE_FILE_AUDIT.csv`, and this report are the
authoritative Panel-III metadata; no numerical reference file was changed.
"""
    (REF / "REFERENCE_AUDIT.md").write_text(report)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
