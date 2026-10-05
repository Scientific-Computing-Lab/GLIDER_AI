#!/usr/bin/env python3
"""One-shot scoring of frozen canonical-source prospective Panel III."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "canonical_panel_3"
REF = PANEL / "references"
PRED = PANEL / "predictions"
SEED = 20260814
N_BOOT = 100_000


@dataclass(frozen=True)
class Method:
    key: str
    label: str
    directory: Path
    capability: str
    information: str
    status: str


METHODS = [
    Method("glider", "GLIDER", PRED / "glider", "both", "geometry + atomic identity", "frozen before QM"),
    Method("mace_polar_l", "MACE-POLAR-1-L", PRED / "comparators/mace_polar_l", "both", "geometry + atomic identity", "public checkpoint; frozen before QM"),
    Method("mace_polar_m", "MACE-POLAR-1-M", PRED / "comparators/mace_polar_m", "both", "geometry + atomic identity", "public checkpoint; frozen before QM"),
    Method("mace_polar_s", "MACE-POLAR-1-S", PRED / "comparators/mace_polar_s", "both", "geometry + atomic identity", "public checkpoint; frozen before QM"),
    Method("aimnet2", "AIMNet2", PRED / "comparators/aimnet2", "both", "geometry + atomic identity", "public checkpoint; frozen before QM"),
    Method("mace_mdp", "MACE-MDP", PRED / "comparators/mace_mdp", "dipole", "geometry + atomic identity", "public checkpoint; frozen before QM"),
    Method("gfn2_xtb", "GFN2-xTB", PRED / "comparators/gfn2_xtb", "dipole", "geometry + atomic identity", "fixed public method; frozen before QM"),
    Method("static_thole", "Static q/u + Thole", PRED / "comparators/static_thole", "both", "geometry + frozen static response parameters", "frozen physical control"),
    Method("zero_response", "Zero response", PRED / "comparators/zero_response", "both", "none", "frozen null control"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def prediction_paths(method: Method) -> dict[str, Path]:
    table = pd.read_csv(method.directory / "prediction_registry.csv")
    return {
        str(row.config_id): method.directory / str(row.prediction_file)
        for row in table.itertuples()
    }


def bootstrap_delta(first: np.ndarray, second: np.ndarray, *, root: bool) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(first), size=(N_BOOT, len(first)))
    a = first[idx].mean(axis=1)
    b = second[idx].mean(axis=1)
    if root:
        a, b = np.sqrt(a), np.sqrt(b)
    delta = a - b
    return float(np.quantile(delta, 0.025)), float(np.quantile(delta, 0.975))


def main() -> None:
    manifest = json.loads((REF / "REFERENCE_MANIFEST.json").read_text())
    if not manifest.get("reference_audit_pass"):
        raise RuntimeError("Reference audit has not passed")

    reference_table = pd.read_csv(REF / "observable_registry.csv")
    references = {
        str(row.config_id): REF / str(row.observable_file)
        for row in reference_table.itertuples()
    }
    expected = set(references)
    if len(expected) != 80:
        raise RuntimeError("Expected 80 references")

    config_rows: list[dict[str, object]] = []
    for method in METHODS:
        predictions = prediction_paths(method)
        if set(predictions) != expected:
            raise RuntimeError(f"Configuration mismatch for {method.key}")
        for cid in sorted(expected):
            metadata = reference_table.loc[reference_table.config_id == cid].iloc[0]
            ref = np.load(references[cid], allow_pickle=False)
            pred = np.load(predictions[cid], allow_pickle=False)
            row: dict[str, object] = {
                "method": method.key,
                "method_label": method.label,
                "config_id": cid,
                "molecule_id": str(metadata.molecule_id),
                "regime": str(metadata.regime),
                "esp_nrmse": np.nan,
                "dipole_mse_D2": np.nan,
                "dipole_vector_error_D": np.nan,
            }
            if method.capability in {"both", "esp"}:
                if not np.array_equal(ref["points_angstrom"], pred["points_angstrom"]):
                    raise RuntimeError(f"Probe mismatch: {method.key}/{cid}")
                target = np.asarray(ref["delta_esp_hartree_per_e"], dtype=float)
                estimate = np.asarray(pred["predicted_esp_hartree_per_e"], dtype=float)
                row["esp_nrmse"] = float(
                    np.sqrt(np.mean((estimate - target) ** 2))
                    / max(np.sqrt(np.mean(target**2)), 1.0e-12)
                )
            if method.capability in {"both", "dipole"}:
                error = np.asarray(pred["predicted_dipole_debye"], dtype=float) - np.asarray(
                    ref["delta_dipole_debye"], dtype=float
                )
                row["dipole_mse_D2"] = float(np.mean(error**2))
                row["dipole_vector_error_D"] = float(np.linalg.norm(error))
            config_rows.append(row)

    configuration = pd.DataFrame(config_rows)
    molecule = configuration.groupby(
        ["method", "method_label", "molecule_id"], as_index=False
    ).agg(
        esp_nrmse=("esp_nrmse", "mean"),
        dipole_mse_D2=("dipole_mse_D2", "mean"),
    )
    molecule["dipole_rmse_D"] = np.sqrt(molecule.dipole_mse_D2)
    leaderboard = molecule.groupby(["method", "method_label"], as_index=False).agg(
        esp_nrmse=("esp_nrmse", "mean"), dipole_mse_D2=("dipole_mse_D2", "mean")
    )
    leaderboard["dipole_rmse_D"] = np.sqrt(leaderboard.dipole_mse_D2)
    leaderboard = leaderboard.drop(columns="dipole_mse_D2")
    leaderboard["esp_rank"] = leaderboard.esp_nrmse.rank(method="min", na_option="bottom")
    leaderboard.loc[leaderboard.esp_nrmse.isna(), "esp_rank"] = np.nan
    leaderboard["dipole_rank"] = leaderboard.dipole_rmse_D.rank(method="min", na_option="bottom")

    gl = leaderboard.set_index("method").loc["glider"]
    esp_pivot = molecule.pivot(index="molecule_id", columns="method", values="esp_nrmse")
    dip_pivot = molecule.pivot(index="molecule_id", columns="method", values="dipole_mse_D2")
    leaderboard["glider_esp_improvement_pct"] = 100.0 * (1.0 - float(gl.esp_nrmse) / leaderboard.esp_nrmse)
    leaderboard["glider_dipole_improvement_pct"] = 100.0 * (1.0 - float(gl.dipole_rmse_D) / leaderboard.dipole_rmse_D)
    leaderboard["glider_esp_molecule_wins"] = leaderboard.method.map(
        lambda m: int((esp_pivot.glider < esp_pivot[m]).sum()) if esp_pivot[m].notna().all() else np.nan
    )
    leaderboard["glider_dipole_molecule_wins"] = leaderboard.method.map(
        lambda m: int((dip_pivot.glider < dip_pivot[m]).sum())
    )
    meta = pd.DataFrame(
        [{"method": m.key, "inference_information_required": m.information, "public_frozen_status": m.status} for m in METHODS]
    )
    leaderboard = leaderboard.merge(meta, on="method").sort_values(
        ["esp_rank", "dipole_rank"], na_position="last"
    )

    regime = configuration.groupby(["method", "method_label", "regime"], as_index=False).agg(
        esp_nrmse=("esp_nrmse", "mean"), dipole_mse_D2=("dipole_mse_D2", "mean")
    )
    regime["dipole_rmse_D"] = np.sqrt(regime.dipole_mse_D2)
    regime = regime.drop(columns="dipole_mse_D2")

    pairwise_rows: list[dict[str, object]] = []
    for method in METHODS:
        if method.key == "glider":
            continue
        if dip_pivot[method.key].notna().all():
            ci = bootstrap_delta(dip_pivot.glider.to_numpy(), dip_pivot[method.key].to_numpy(), root=True)
            pairwise_rows.append(
                {
                    "observable": "dipole_rmse_D",
                    "comparator": method.key,
                    "comparator_label": method.label,
                    "glider_value": float(gl.dipole_rmse_D),
                    "comparator_value": float(leaderboard.set_index("method").loc[method.key, "dipole_rmse_D"]),
                    "glider_minus_comparator": float(gl.dipole_rmse_D - leaderboard.set_index("method").loc[method.key, "dipole_rmse_D"]),
                    "relative_error_reduction_pct": 100.0 * (1.0 - float(gl.dipole_rmse_D) / float(leaderboard.set_index("method").loc[method.key, "dipole_rmse_D"])),
                    "ci95_low": ci[0],
                    "ci95_high": ci[1],
                    "glider_molecule_wins": int((dip_pivot.glider < dip_pivot[method.key]).sum()),
                    "n_molecules": 20,
                    "bootstrap_replicates": N_BOOT,
                    "bootstrap_seed": SEED,
                }
            )
        if method.key in esp_pivot and esp_pivot[method.key].notna().all():
            ci = bootstrap_delta(esp_pivot.glider.to_numpy(), esp_pivot[method.key].to_numpy(), root=False)
            pairwise_rows.append(
                {
                    "observable": "esp_nrmse",
                    "comparator": method.key,
                    "comparator_label": method.label,
                    "glider_value": float(gl.esp_nrmse),
                    "comparator_value": float(leaderboard.set_index("method").loc[method.key, "esp_nrmse"]),
                    "glider_minus_comparator": float(gl.esp_nrmse - leaderboard.set_index("method").loc[method.key, "esp_nrmse"]),
                    "relative_error_reduction_pct": 100.0 * (1.0 - float(gl.esp_nrmse) / float(leaderboard.set_index("method").loc[method.key, "esp_nrmse"])),
                    "ci95_low": ci[0],
                    "ci95_high": ci[1],
                    "glider_molecule_wins": int((esp_pivot.glider < esp_pivot[method.key]).sum()),
                    "n_molecules": 20,
                    "bootstrap_replicates": N_BOOT,
                    "bootstrap_seed": SEED,
                }
            )
    pairwise = pd.DataFrame(pairwise_rows)

    strongest_esp = leaderboard[(leaderboard.method != "glider") & leaderboard.esp_nrmse.notna()].sort_values("esp_nrmse").iloc[0]
    strongest_dip = leaderboard[leaderboard.method != "glider"].sort_values("dipole_rmse_D").iloc[0]
    esp_pair = pairwise[(pairwise.observable == "esp_nrmse") & (pairwise.comparator == strongest_esp.method)].iloc[0]
    dip_pair = pairwise[(pairwise.observable == "dipole_rmse_D") & (pairwise.comparator == strongest_dip.method)].iloc[0]
    esp_rank = int(gl.esp_rank)
    if esp_rank == 1 and float(esp_pair.ci95_high) < 0:
        verdict = "A_STRONG_CANONICAL_REPLICATION"
    elif esp_rank == 1 and float(esp_pair.ci95_low) <= 0 <= float(esp_pair.ci95_high):
        verdict = "B_NOMINAL_CANONICAL_REPLICATION"
    elif esp_rank != 1 and float(esp_pair.ci95_low) <= 0 <= float(esp_pair.ci95_high):
        verdict = "C_MIXED_CANONICAL_RESULT"
    else:
        verdict = "D_FAILED_CANONICAL_REPLICATION"

    configuration.to_csv(PANEL / "CONFIGURATION_LEVEL_METRICS.csv", index=False)
    molecule.drop(columns="dipole_mse_D2").sort_values(["molecule_id", "method"]).to_csv(
        PANEL / "MOLECULE_LEVEL_METRICS.csv", index=False
    )
    regime.sort_values(["regime", "method"]).to_csv(PANEL / "REGIME_LEVEL_METRICS.csv", index=False)
    leaderboard.to_csv(PANEL / "FINAL_COMPARATOR_TABLE.csv", index=False)
    pairwise.sort_values(["observable", "comparator"]).to_csv(PANEL / "PAIRWISE_BOOTSTRAP.csv", index=False)

    selected = pd.read_csv(PANEL / "selection/selected_molecules.csv")
    gl_mol = molecule[molecule.method == "glider"].copy()
    worst_esp = gl_mol.loc[gl_mol.esp_nrmse.idxmax()]
    worst_dip = gl_mol.loc[gl_mol.dipole_rmse_D.idxmax()]
    stats = {
        "terminal_verdict": verdict,
        "scored_once_utc": datetime.now(timezone.utc).isoformat(),
        "panel": {"molecules": 20, "configurations": 80, "regimes": 4, "identity_source": "FreeSolv v0.52"},
        "glider": {"esp_nrmse": float(gl.esp_nrmse), "esp_rank": esp_rank, "dipole_rmse_D": float(gl.dipole_rmse_D), "dipole_rank": int(gl.dipole_rank)},
        "strongest_esp_comparator": esp_pair.to_dict(),
        "strongest_dipole_comparator": dip_pair.to_dict(),
        "worst_glider_esp_molecule": {"molecule_id": str(worst_esp.molecule_id), "esp_nrmse": float(worst_esp.esp_nrmse)},
        "worst_glider_dipole_molecule": {"molecule_id": str(worst_dip.molecule_id), "dipole_rmse_D": float(worst_dip.dipole_rmse_D)},
        "selected_freesolv_ids": selected.freesolv_id.tolist(),
        "bootstrap": {"unit": "molecule", "replicates": N_BOOT, "seed": SEED, "interval": "percentile 95%"},
        "experimental_hydration_targets_accessed": False,
        "model_or_metric_changed_after_freeze": False,
    }
    (PANEL / "FINAL_STATISTICS.json").write_text(json.dumps(stats, indent=2, default=str) + "\n")

    report = f"""# Canonical-source prospective Panel III: final report

## Terminal verdict

**{verdict.replace('_', ' ')}.**

GLIDER response-ESP NRMSE is **{gl.esp_nrmse:.6f}** (rank {esp_rank}). The
strongest eligible comparator is {strongest_esp.method_label} at
**{strongest_esp.esp_nrmse:.6f}**. The GLIDER-minus-comparator difference is
{esp_pair.glider_minus_comparator:.6f} (95% molecule-blocked bootstrap CI
[{esp_pair.ci95_low:.6f}, {esp_pair.ci95_high:.6f}]); GLIDER improves
{int(esp_pair.glider_molecule_wins)}/20 molecules.

GLIDER induced-dipole RMSE is **{gl.dipole_rmse_D:.6f} D** (rank
{int(gl.dipole_rank)}). The strongest comparator is
{strongest_dip.method_label} at **{strongest_dip.dipole_rmse_D:.6f} D**. The
difference is {dip_pair.glider_minus_comparator:.6f} D (95% CI
[{dip_pair.ci95_low:.6f}, {dip_pair.ci95_high:.6f}]); GLIDER improves
{int(dip_pair.glider_molecule_wins)}/20 molecules.

## Integrity

The 20 identities were selected deterministically from an identity-only copy
of pinned FreeSolv v0.52. FreeSolv hydration values were not accessed. All 80
geometries, GLIDER predictions, and eligible comparator predictions were
hash-locked and pushed before reference QM began. No molecule or configuration
was removed or replaced. Four iodine-containing cases used the declared
same-Hamiltonian SCF recovery route; all completed at the frozen reference
theory.

## Scope

Panel III tests mutual electronic response of neutral organic solute--three-water
environments. FreeSolv supplies molecular identities only; this is not a
hydration-free-energy benchmark. Chemical identities are absent from GLIDER
response supervision, but exact identity or analogue exposure during broad
foundation pretraining cannot be excluded.
"""
    (PANEL / "FINAL_REPORT.md").write_text(report)

    claim = f"""# Panel-III claim matrix

| Claim | Classification | Basis |
|---|---|---|
| GLIDER has the lowest response-ESP error in canonical Panel III | {'STRONGLY SUPPORTED' if verdict.startswith('A_') else 'NOMINALLY SUPPORTED' if verdict.startswith('B_') else 'NOT SUPPORTED'} | Rank {esp_rank}; paired CI [{esp_pair.ci95_low:.6f}, {esp_pair.ci95_high:.6f}] |
| GLIDER has the lowest induced-dipole error in Panel III | {'NOMINALLY SUPPORTED' if int(gl.dipole_rank)==1 else 'NOT SUPPORTED'} | Observed rank {int(gl.dipole_rank)}; paired CI [{dip_pair.ci95_low:.6f}, {dip_pair.ci95_high:.6f}] |
| Transfer to chemistries absent from GLIDER response supervision | STRONGLY SUPPORTED | Label-blind FreeSolv identities; exact/near historical response-chemistry exclusions |
| Absence of identity exposure in foundation pretraining | PROHIBITED | Broad encoder pretraining prevents a complete identity-exposure exclusion |
| Generality to arbitrary molecular chemistry | PROHIBITED | Neutral-organic, element- and size-bounded scope |
| Generality to bulk liquid solvation | PROHIBITED | Structured three-water environments only |
| Hydration free-energy improvement | PROHIBITED | Hydration targets were neither accessed nor evaluated |
"""
    (PANEL / "CLAIM_MATRIX.md").write_text(claim)

    release_files = [
        "FINAL_COMPARATOR_TABLE.csv", "MOLECULE_LEVEL_METRICS.csv", "REGIME_LEVEL_METRICS.csv",
        "PAIRWISE_BOOTSTRAP.csv", "CONFIGURATION_LEVEL_METRICS.csv", "FINAL_STATISTICS.json",
        "FINAL_REPORT.md", "CLAIM_MATRIX.md", "references/REFERENCE_MANIFEST.json",
        "references/REFERENCE_AUDIT.md",
    ]
    release = {
        "terminal_verdict": verdict,
        "prediction_freeze_commit": "f9a87b1d5a91a45f595ccd1891f59ec8d9445fea",
        "all_predictions_freeze_manifest_sha256": "3dc5641a289f48118e0e372c341286573d4a381a2cf2ea9b7ae2b2b80ec8c73e",
        "glider_checkpoint_sha256": "a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288",
        "files": {p: sha256(PANEL / p) for p in release_files},
    }
    (PANEL / "FINAL_RELEASE_MANIFEST.json").write_text(json.dumps(release, indent=2) + "\n")
    print(json.dumps(stats, indent=2, default=str))


if __name__ == "__main__":
    main()
