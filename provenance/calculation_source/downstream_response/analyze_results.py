#!/usr/bin/env python3
"""Aggregate the frozen downstream-response experiment and make figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from common import sha256  # noqa: E402
from evaluate_couplings import METHODS  # noqa: E402
from scipy.spatial.transform import Rotation  # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response"
BOOTSTRAP_SEED = 2026081901
N_BOOTSTRAP = 100_000
LABELS = {
    "qm": "QM",
    "glider": "GLIDER",
    "mace_polar_l": "MACE-POLAR-1-L",
    "exact_dipole_base_com": "Exact dipole (base COM)",
    "exact_dipole_solute_com": "Exact dipole (solute COM)",
    "exact_dipole_nuclear_center": "Exact dipole (nuclear centre)",
    "zero": "Zero response",
}
COLORS = {
    "qm": "#202124",
    "glider": "#1565C0",
    "mace_polar_l": "#D97706",
    "exact_dipole_base_com": "#218739",
    "zero": "#777777",
}
PRIMARY_COMPARATORS = (
    "glider",
    "mace_polar_l",
    "exact_dipole_base_com",
    "zero",
)


def safe_correlations(reference: np.ndarray, predicted: np.ndarray) -> tuple[float, float]:
    reference = np.asarray(reference, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if np.std(reference) == 0 or np.std(predicted) == 0:
        return float("nan"), float("nan")
    return float(pearsonr(reference, predicted).statistic), float(
        spearmanr(reference, predicted).statistic
    )


def interval(values: np.ndarray) -> tuple[float, float]:
    return tuple(float(value) for value in np.quantile(values, (0.025, 0.975)))


def energy_statistics(table: pd.DataFrame, output: Path) -> dict[str, object]:
    reference = table["qm_energy_kcal_mol"].to_numpy()
    aggregate = []
    per_solute = []
    for method in METHODS:
        predicted = table[f"{method}_energy_kcal_mol"].to_numpy()
        error = predicted - reference
        pearson, spearman = safe_correlations(reference, predicted)
        aggregate.append(
            {
                "method": method,
                "label": LABELS[method],
                "n_cases": len(table),
                "MAE_kcal_mol": float(np.mean(np.abs(error))),
                "RMSE_kcal_mol": float(np.sqrt(np.mean(error**2))),
                "signed_bias_kcal_mol": float(np.mean(error)),
                "Pearson_r": pearson,
                "Spearman_rho": spearman,
            }
        )
        table[f"{method}_absolute_error_kcal_mol"] = np.abs(error)
        for molecule_id, group in table.groupby("molecule_id", sort=True):
            group_reference = group["qm_energy_kcal_mol"].to_numpy()
            group_prediction = group[f"{method}_energy_kcal_mol"].to_numpy()
            group_error = group_prediction - group_reference
            per_solute.append(
                {
                    "molecule_id": molecule_id,
                    "method": method,
                    "n_cases": len(group),
                    "MAE_kcal_mol": float(np.mean(np.abs(group_error))),
                    "RMSE_kcal_mol": float(np.sqrt(np.mean(group_error**2))),
                    "signed_bias_kcal_mol": float(np.mean(group_error)),
                }
            )
    table.to_csv(output / "experiment1_per_case.csv", index=False)
    aggregate_table = pd.DataFrame(aggregate)
    per_solute_table = pd.DataFrame(per_solute)
    aggregate_table.to_csv(output / "energy_aggregate_statistics.csv", index=False)
    per_solute_table.to_csv(output / "energy_per_solute_statistics.csv", index=False)

    solutes = sorted(table["molecule_id"].unique())
    rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    sampled = rng.integers(0, len(solutes), size=(N_BOOTSTRAP, len(solutes)))
    bootstrap_rows = []
    per_method_bootstrap: dict[str, dict[str, np.ndarray]] = {}
    for method in METHODS:
        solute_mae = []
        solute_mse = []
        solute_bias = []
        for molecule_id in solutes:
            group = table[table["molecule_id"] == molecule_id]
            error = (
                group[f"{method}_energy_kcal_mol"]
                - group["qm_energy_kcal_mol"]
            ).to_numpy()
            solute_mae.append(np.mean(np.abs(error)))
            solute_mse.append(np.mean(error**2))
            solute_bias.append(np.mean(error))
        mae_values = np.mean(np.asarray(solute_mae)[sampled], axis=1)
        rmse_values = np.sqrt(np.mean(np.asarray(solute_mse)[sampled], axis=1))
        bias_values = np.mean(np.asarray(solute_bias)[sampled], axis=1)
        per_method_bootstrap[method] = {
            "MAE_kcal_mol": mae_values,
            "RMSE_kcal_mol": rmse_values,
            "signed_bias_kcal_mol": bias_values,
        }
        for metric, values in per_method_bootstrap[method].items():
            low, high = interval(values)
            bootstrap_rows.append(
                {
                    "quantity": "method_metric",
                    "method": method,
                    "comparator": "",
                    "metric": metric,
                    "point_estimate": float(np.mean(values)),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    for comparator in METHODS:
        if comparator == "glider":
            continue
        for metric in ("MAE_kcal_mol", "RMSE_kcal_mol"):
            delta = (
                per_method_bootstrap["glider"][metric]
                - per_method_bootstrap[comparator][metric]
            )
            low, high = interval(delta)
            bootstrap_rows.append(
                {
                    "quantity": "paired_glider_minus_comparator",
                    "method": "glider",
                    "comparator": comparator,
                    "metric": metric,
                    "point_estimate": float(np.mean(delta)),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    pd.DataFrame(bootstrap_rows).to_csv(
        output / "blocked_bootstrap_intervals.csv", index=False
    )

    leave_one_out = []
    for held_out in solutes:
        subset = table[table["molecule_id"] != held_out]
        subset_reference = subset["qm_energy_kcal_mol"].to_numpy()
        for method in METHODS:
            error = (
                subset[f"{method}_energy_kcal_mol"].to_numpy() - subset_reference
            )
            leave_one_out.append(
                {
                    "held_out_solute": held_out,
                    "method": method,
                    "n_solutes": len(solutes) - 1,
                    "n_cases": len(subset),
                    "MAE_kcal_mol": float(np.mean(np.abs(error))),
                    "RMSE_kcal_mol": float(np.sqrt(np.mean(error**2))),
                    "signed_bias_kcal_mol": float(np.mean(error)),
                }
            )
    pd.DataFrame(leave_one_out).to_csv(
        output / "leave_one_solute_out.csv", index=False
    )

    full = table["full_frozen_electrostatic_kcal_mol"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.abs(reference) / np.abs(full)
    physical = {
        "n_cases": len(table),
        "reference_rms_kcal_mol": float(np.sqrt(np.mean(reference**2))),
        "reference_median_absolute_kcal_mol": float(np.median(np.abs(reference))),
        "reference_min_kcal_mol": float(np.min(reference)),
        "reference_max_kcal_mol": float(np.max(reference)),
        "reference_max_absolute_kcal_mol": float(np.max(np.abs(reference))),
        "kBT_298.15K_kcal_mol": 0.592484949,
        "reference_rms_over_kBT": float(
            np.sqrt(np.mean(reference**2)) / 0.592484949
        ),
        "median_abs_response_to_abs_full_frozen_electrostatic": float(
            np.nanmedian(ratios)
        ),
        "response_to_full_ratio_interpretation": "secondary scale only; the denominator is not a total interaction energy",
    }
    (output / "physical_scale.json").write_text(
        json.dumps(physical, indent=2, sort_keys=True) + "\n"
    )
    return {
        "aggregate": aggregate_table,
        "per_solute": per_solute_table,
        "physical": physical,
    }


def vector_statistics(
    table: pd.DataFrame,
    output: Path,
    observable: str,
    angular_threshold: float,
) -> pd.DataFrame:
    rows = []
    per_solute_rows = []
    for method, group in table.groupby("method", sort=False):
        meaningful = group[
            (group["reference_magnitude"] >= angular_threshold)
            & np.isfinite(group["angular_error_degrees"])
        ]
        rows.append(
            {
                "observable": observable,
                "method": method,
                "n_cases": len(group),
                "vector_RMSE": float(
                    np.sqrt(np.mean(group["vector_error_norm"] ** 2))
                ),
                "magnitude_MAE": float(np.mean(np.abs(group["magnitude_error"]))),
                "magnitude_RMSE": float(
                    np.sqrt(np.mean(group["magnitude_error"] ** 2))
                ),
                "mean_cosine_similarity": float(
                    meaningful["cosine_similarity"].mean()
                ),
                "median_angular_error_degrees": float(
                    meaningful["angular_error_degrees"].median()
                ),
                "n_directionally_meaningful": len(meaningful),
                "reference_RMS_magnitude": float(
                    np.sqrt(np.mean(group["reference_magnitude"] ** 2))
                ),
                "reference_median_magnitude": float(
                    np.median(group["reference_magnitude"])
                ),
            }
        )
        for molecule_id, solute_group in group.groupby("molecule_id", sort=True):
            per_solute_rows.append(
                {
                    "observable": observable,
                    "molecule_id": molecule_id,
                    "method": method,
                    "n_cases": len(solute_group),
                    "vector_RMSE": float(
                        np.sqrt(np.mean(solute_group["vector_error_norm"] ** 2))
                    ),
                    "vector_MSE": float(
                        np.mean(solute_group["vector_error_norm"] ** 2)
                    ),
                    "magnitude_MAE": float(
                        np.mean(np.abs(solute_group["magnitude_error"]))
                    ),
                    "magnitude_MSE": float(
                        np.mean(solute_group["magnitude_error"] ** 2)
                    ),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output / f"{observable}_aggregate_statistics.csv", index=False)
    per_solute = pd.DataFrame(per_solute_rows)
    per_solute.to_csv(
        output / f"{observable}_per_solute_statistics.csv", index=False
    )

    solutes = sorted(per_solute["molecule_id"].unique())
    rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    sampled = rng.integers(0, len(solutes), size=(N_BOOTSTRAP, len(solutes)))
    method_bootstrap: dict[str, dict[str, np.ndarray]] = {}
    bootstrap_rows = []
    for method in METHODS:
        values = per_solute[per_solute["method"] == method].set_index("molecule_id")
        values = values.loc[solutes]
        vector_rmse = np.sqrt(
            np.mean(values["vector_MSE"].to_numpy()[sampled], axis=1)
        )
        magnitude_mae = np.mean(
            values["magnitude_MAE"].to_numpy()[sampled], axis=1
        )
        method_bootstrap[method] = {
            "vector_RMSE": vector_rmse,
            "magnitude_MAE": magnitude_mae,
        }
        direct = result[result["method"] == method].iloc[0]
        for metric, bootstrap_values in method_bootstrap[method].items():
            low, high = interval(bootstrap_values)
            bootstrap_rows.append(
                {
                    "observable": observable,
                    "quantity": "method_metric",
                    "method": method,
                    "comparator": "",
                    "metric": metric,
                    "point_estimate": float(direct[metric]),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    for comparator in METHODS:
        if comparator == "glider":
            continue
        for metric in ("vector_RMSE", "magnitude_MAE"):
            delta = (
                method_bootstrap["glider"][metric]
                - method_bootstrap[comparator][metric]
            )
            low, high = interval(delta)
            point = float(
                result.loc[result["method"] == "glider", metric].iloc[0]
                - result.loc[result["method"] == comparator, metric].iloc[0]
            )
            bootstrap_rows.append(
                {
                    "observable": observable,
                    "quantity": "paired_glider_minus_comparator",
                    "method": "glider",
                    "comparator": comparator,
                    "metric": metric,
                    "point_estimate": point,
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    pd.DataFrame(bootstrap_rows).to_csv(
        output / f"{observable}_blocked_bootstrap_intervals.csv", index=False
    )
    return result


def orientation_statistics(table: pd.DataFrame, output: Path) -> pd.DataFrame:
    rows = []
    for case_id, case_group in table.groupby("case_id", sort=True):
        reference_table = case_group[case_group["method"] == "qm"].sort_values(
            "rotation_index"
        )
        reference = reference_table["energy_kcal_mol"].to_numpy()
        reference_amplitude = float(np.ptp(reference))
        rotations = {
            row.rotation_id: Rotation.from_quat(
                [row.quaternion_x, row.quaternion_y, row.quaternion_z, row.quaternion_w]
            )
            for row in reference_table.itertuples()
        }
        reference_minimum = reference_table.iloc[int(np.argmin(reference))][
            "rotation_id"
        ]
        reference_maximum = reference_table.iloc[int(np.argmax(reference))][
            "rotation_id"
        ]
        for method, method_group in case_group.groupby("method", sort=False):
            method_group = method_group.sort_values("rotation_index")
            predicted = method_group["energy_kcal_mol"].to_numpy()
            error = predicted - reference
            pearson, spearman = safe_correlations(reference, predicted)
            predicted_minimum = method_group.iloc[int(np.argmin(predicted))][
                "rotation_id"
            ]
            predicted_maximum = method_group.iloc[int(np.argmax(predicted))][
                "rotation_id"
            ]
            if reference_amplitude >= 0.05:
                minimum_angle = float(
                    (rotations[reference_minimum].inv() * rotations[predicted_minimum]).magnitude()
                    * 180.0
                    / np.pi
                )
                maximum_angle = float(
                    (rotations[reference_maximum].inv() * rotations[predicted_maximum]).magnitude()
                    * 180.0
                    / np.pi
                )
            else:
                minimum_angle = float("nan")
                maximum_angle = float("nan")
            rows.append(
                {
                    "case_id": case_id,
                    "molecule_id": method_group.iloc[0]["molecule_id"],
                    "method": method,
                    "profile_RMSE_kcal_mol": float(np.sqrt(np.mean(error**2))),
                    "profile_MAE_kcal_mol": float(np.mean(np.abs(error))),
                    "profile_bias_kcal_mol": float(np.mean(error)),
                    "Pearson_r": pearson,
                    "Spearman_rho": spearman,
                    "reference_amplitude_kcal_mol": reference_amplitude,
                    "predicted_amplitude_kcal_mol": float(np.ptp(predicted)),
                    "amplitude_error_kcal_mol": float(
                        np.ptp(predicted) - reference_amplitude
                    ),
                    "minimum_geodesic_error_degrees": minimum_angle,
                    "maximum_geodesic_error_degrees": maximum_angle,
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output / "orientation_statistics.csv", index=False)
    aggregate = (
        result.groupby("method", sort=False)
        .agg(
            n_cases=("case_id", "size"),
            mean_profile_RMSE_kcal_mol=("profile_RMSE_kcal_mol", "mean"),
            rms_profile_error_kcal_mol=(
                "profile_RMSE_kcal_mol",
                lambda value: float(np.sqrt(np.mean(np.asarray(value) ** 2))),
            ),
            mean_absolute_amplitude_error_kcal_mol=(
                "amplitude_error_kcal_mol",
                lambda value: float(np.mean(np.abs(value))),
            ),
            median_Pearson_r=("Pearson_r", "median"),
            median_Spearman_rho=("Spearman_rho", "median"),
        )
        .reset_index()
    )
    aggregate.to_csv(output / "orientation_aggregate_statistics.csv", index=False)
    solutes = sorted(result["molecule_id"].unique())
    rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    sampled = rng.integers(0, len(solutes), size=(N_BOOTSTRAP, len(solutes)))
    bootstrap_values = {}
    bootstrap_rows = []
    for method in METHODS:
        method_values = (
            result[result["method"] == method]
            .set_index("molecule_id")
            .loc[solutes, "profile_RMSE_kcal_mol"]
            .to_numpy()
        )
        values = np.mean(method_values[sampled], axis=1)
        bootstrap_values[method] = values
        low, high = interval(values)
        point = float(
            aggregate.loc[
                aggregate["method"] == method, "mean_profile_RMSE_kcal_mol"
            ].iloc[0]
        )
        bootstrap_rows.append(
            {
                "quantity": "method_metric",
                "method": method,
                "comparator": "",
                "metric": "mean_profile_RMSE_kcal_mol",
                "point_estimate": point,
                "ci95_low": low,
                "ci95_high": high,
            }
        )
    for comparator in METHODS:
        if comparator == "glider":
            continue
        delta = bootstrap_values["glider"] - bootstrap_values[comparator]
        low, high = interval(delta)
        glider_point = aggregate.loc[
            aggregate["method"] == "glider", "mean_profile_RMSE_kcal_mol"
        ].iloc[0]
        comparator_point = aggregate.loc[
            aggregate["method"] == comparator, "mean_profile_RMSE_kcal_mol"
        ].iloc[0]
        bootstrap_rows.append(
            {
                "quantity": "paired_glider_minus_comparator",
                "method": "glider",
                "comparator": comparator,
                "metric": "mean_profile_RMSE_kcal_mol",
                "point_estimate": float(glider_point - comparator_point),
                "ci95_low": low,
                "ci95_high": high,
            }
        )
    pd.DataFrame(bootstrap_rows).to_csv(
        output / "orientation_blocked_bootstrap_intervals.csv", index=False
    )
    return result


def save_figure(figure: plt.Figure, stem: Path) -> None:
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def make_energy_figures(table: pd.DataFrame, figures: Path) -> None:
    methods = ("glider", "mace_polar_l", "exact_dipole_base_com")
    reference = table["qm_energy_kcal_mol"].to_numpy()
    all_values = [reference]
    for method in methods:
        all_values.append(table[f"{method}_energy_kcal_mol"].to_numpy())
    lower = min(np.min(values) for values in all_values)
    upper = max(np.max(values) for values in all_values)
    padding = 0.08 * (upper - lower or 1.0)
    lower -= padding
    upper += padding
    figure, axis = plt.subplots(figsize=(5.2, 4.7))
    for method, marker in zip(methods, ("o", "s", "^"), strict=True):
        axis.scatter(
            reference,
            table[f"{method}_energy_kcal_mol"],
            label=LABELS[method],
            color=COLORS[method],
            marker=marker,
            s=36,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.4,
        )
    axis.plot([lower, upper], [lower, upper], color="#444444", lw=1, ls="--")
    axis.axhline(0, color="#BBBBBB", lw=0.7)
    axis.axvline(0, color="#BBBBBB", lw=0.7)
    axis.set(xlim=(lower, upper), ylim=(lower, upper))
    axis.set_xlabel("QM response coupling (kcal/mol)")
    axis.set_ylabel("Predicted response coupling (kcal/mol)")
    axis.legend(frameon=False, fontsize=8)
    axis.set_aspect("equal", adjustable="box")
    figure.tight_layout()
    save_figure(figure, figures / "energy_scatter")

    figure, axis = plt.subplots(figsize=(5.8, 4.6))
    x = np.arange(len(methods))
    errors = np.column_stack(
        [
            np.abs(table[f"{method}_energy_kcal_mol"].to_numpy() - reference)
            for method in methods
        ]
    )
    for values in errors:
        axis.plot(x, values, color="#B8BDC5", lw=0.7, alpha=0.55)
    for index, method in enumerate(methods):
        axis.scatter(
            np.full(len(table), index),
            errors[:, index],
            color=COLORS[method],
            s=20,
            zorder=3,
        )
        axis.scatter(
            index,
            np.mean(errors[:, index]),
            color="white",
            edgecolor="black",
            marker="D",
            s=44,
            zorder=4,
        )
    axis.set_xticks(x, [LABELS[method] for method in methods], rotation=12)
    axis.set_ylabel("Absolute error (kcal/mol)")
    axis.set_title("Paired held-out-water response-coupling errors")
    figure.tight_layout()
    save_figure(figure, figures / "paired_energy_errors")


def make_vector_figure(
    force_stats: pd.DataFrame, torque_stats: pd.DataFrame, figures: Path
) -> None:
    methods = ("glider", "mace_polar_l", "exact_dipole_base_com")
    figure, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
    for axis, table, title, unit in (
        (axes[0], force_stats, "Response force", "kcal/mol/Angstrom"),
        (axes[1], torque_stats, "Response torque", "kcal/mol"),
    ):
        subset = table.set_index("method").loc[list(methods)]
        positions = np.arange(len(methods))
        axis.bar(
            positions - 0.18,
            subset["vector_RMSE"],
            width=0.36,
            color=[COLORS[method] for method in methods],
            label="Vector RMSE",
        )
        axis.bar(
            positions + 0.18,
            subset["magnitude_MAE"],
            width=0.36,
            color=[COLORS[method] for method in methods],
            alpha=0.45,
            hatch="//",
            label="Magnitude MAE",
        )
        axis.set_xticks(positions, ["GLIDER", "MACE-L", "Exact dipole"])
        axis.set_ylabel(unit)
        axis.set_title(title)
    axes[0].legend(frameon=False, fontsize=8)
    figure.tight_layout()
    save_figure(figure, figures / "force_torque_diagnostics")


def make_orientation_figure(
    profiles: pd.DataFrame, statistics: pd.DataFrame, figures: Path
) -> None:
    case_ids = sorted(profiles["case_id"].unique())
    figure = plt.figure(figsize=(10.0, 8.0))
    grid = figure.add_gridspec(3, 3, height_ratios=(1, 1, 0.85))
    methods = ("qm", "glider", "mace_polar_l", "exact_dipole_base_com")
    for index, case_id in enumerate(case_ids):
        axis = figure.add_subplot(grid[index // 3, index % 3])
        subset = profiles[profiles["case_id"] == case_id]
        for method in methods:
            values = subset[subset["method"] == method].sort_values("rotation_index")
            axis.plot(
                values["rotation_index"],
                values["energy_kcal_mol"],
                color=COLORS[method],
                lw=1.25,
                label=LABELS[method],
            )
        axis.axhline(0, color="#CCCCCC", lw=0.6)
        axis.set_title(str(subset.iloc[0]["molecule_id"]).replace("MNSOL-CONN-", ""), fontsize=8)
        axis.set_xlabel("Octahedral rotation index", fontsize=7)
        axis.set_ylabel("kcal/mol", fontsize=7)
        axis.tick_params(labelsize=7)
        if index == 0:
            axis.legend(frameon=False, fontsize=6)
    summary_axis = figure.add_subplot(grid[2, :])
    summary = (
        statistics[statistics["method"].isin(methods[1:])]
        .groupby("method", sort=False)["profile_RMSE_kcal_mol"]
        .mean()
        .reindex(methods[1:])
    )
    summary_axis.bar(
        np.arange(len(summary)),
        summary.to_numpy(),
        color=[COLORS[method] for method in summary.index],
    )
    summary_axis.set_xticks(
        np.arange(len(summary)), [LABELS[method] for method in summary.index]
    )
    summary_axis.set_ylabel("Mean profile RMSE (kcal/mol)")
    summary_axis.set_title("Orientation-profile error across six fixed parent-00 cases")
    figure.tight_layout()
    save_figure(figure, figures / "orientation_profiles")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluated", type=Path, default=RESULT_ROOT / "evaluated")
    parser.add_argument("--output", type=Path, default=RESULT_ROOT / "analysis")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    figures = args.output / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    energy = pd.read_csv(args.evaluated / "experiment1_per_case.csv")
    energy_result = energy_statistics(energy, args.output)
    force = pd.read_csv(args.evaluated / "force_per_case.csv")
    torque = pd.read_csv(args.evaluated / "torque_per_case.csv")
    force_stats = vector_statistics(force, args.output, "force", 0.05)
    torque_stats = vector_statistics(torque, args.output, "torque", 0.01)
    profiles = pd.read_csv(args.evaluated / "orientation_profiles.csv")
    orientation_stats = orientation_statistics(profiles, args.output)

    make_energy_figures(energy, figures)
    make_vector_figure(force_stats, torque_stats, figures)
    make_orientation_figure(profiles, orientation_stats, figures)
    outputs = [path for path in args.output.rglob("*") if path.is_file()]
    manifest = {
        "experiment": "post hoc downstream response analysis",
        "protocol": "results/posthoc_downstream_response/protocol_freeze.json",
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_resamples": N_BOOTSTRAP,
        "n_cases": len(energy),
        "n_solutes": int(energy["molecule_id"].nunique()),
        "physical_scale": energy_result["physical"],
        "outputs": {
            str(path.relative_to(args.output)): sha256(path)
            for path in sorted(outputs)
            if path.name != "analysis_manifest.json"
        },
        "prospective_claim": False,
    }
    (args.output / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
