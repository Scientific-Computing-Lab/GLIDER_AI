#!/usr/bin/env python3
"""Aggregate confirmatory breadth, spatial-range, and directionality results."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Embed TrueType outlines in publication PDFs rather than Matplotlib's Type 3 fonts.
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.spatial.transform import Rotation  # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

HERE = Path(__file__).resolve().parent
COMMON = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(1, str(COMMON))
from common import sha256  # noqa: E402
from evaluate_couplings import METHODS, RESULT_ROOT  # noqa: E402

BOOTSTRAP_SEED = 2026081902
N_BOOTSTRAP = 100_000
DIPOLE_METHODS = (
    "exact_dipole_base_com",
    "exact_dipole_solute_com",
    "exact_dipole_nuclear_center",
)
LABELS = {
    "qm": "QM response",
    "glider": "GLIDER",
    "mace_polar_l": "MACE-POLAR-1-L",
    "exact_dipole_base_com": "Exact dipole (base COM)",
    "exact_dipole_solute_com": "Exact dipole (solute COM)",
    "exact_dipole_nuclear_center": "Exact dipole (nuclear centre)",
    "zero": "Zero response",
    "full_frozen_electrostatic": "Full frozen electrostatic",
}
COLORS = {
    "glider": "#1769AA",
    "mace_polar_l": "#D97706",
    "exact_dipole_base_com": "#218739",
    "exact_dipole_solute_com": "#218739",
    "exact_dipole_nuclear_center": "#218739",
    "zero": "#777777",
}


def safe_correlations(reference: np.ndarray, predicted: np.ndarray) -> tuple[float, float]:
    reference = np.asarray(reference, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if len(reference) < 2 or np.std(reference) == 0 or np.std(predicted) == 0:
        return float("nan"), float("nan")
    return float(pearsonr(reference, predicted).statistic), float(
        spearmanr(reference, predicted).statistic
    )


def ci(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, (0.025, 0.975))
    return float(low), float(high)


def error_metrics(reference: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    error = np.asarray(prediction) - np.asarray(reference)
    pearson, spearman = safe_correlations(reference, prediction)
    return {
        "MAE": float(np.mean(np.abs(error))),
        "RMSE": float(np.sqrt(np.mean(error**2))),
        "signed_bias": float(np.mean(error)),
        "Pearson_r": pearson,
        "Spearman_rho": spearman,
    }


def best_dipole_origin(energy: pd.DataFrame) -> str:
    reference = energy["qm_energy_kcal_mol"].to_numpy()
    maes = {
        method: np.mean(np.abs(energy[f"{method}_energy_kcal_mol"].to_numpy() - reference))
        for method in DIPOLE_METHODS
    }
    return min(maes, key=lambda method: (maes[method], method))


def best_vector_dipole_origin(table: pd.DataFrame) -> str:
    values = {
        method: float(
            np.sqrt(np.mean(table[table.method == method].vector_error_norm**2))
        )
        for method in DIPOLE_METHODS
    }
    return min(values, key=lambda method: (values[method], method))


def best_orientation_dipole_origin(
    table: pd.DataFrame, metric: str
) -> str:
    values = {
        method: float(table[table.method == method][metric].mean())
        for method in DIPOLE_METHODS
    }
    return min(values, key=lambda method: (values[method], method))


def bootstrap_indices(n: int) -> np.ndarray:
    rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    return rng.integers(0, n, size=(N_BOOTSTRAP, n))


def breadth_energy(
    table: pd.DataFrame, output: Path, best_dipole: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reference = table["qm_energy_kcal_mol"].to_numpy()
    rows = []
    per_rows = []
    for method in METHODS:
        prediction = table[f"{method}_energy_kcal_mol"].to_numpy()
        metrics = error_metrics(reference, prediction)
        rows.append({"endpoint": "W4_energy", "method": method, "n_solutes": len(table), **metrics})
        for row, error in zip(table.itertuples(), prediction - reference, strict=True):
            per_rows.append(
                {
                    "molecule_id": row.molecule_id,
                    "family": row.family,
                    "method": method,
                    "error_kcal_mol": error,
                    "absolute_error_kcal_mol": abs(error),
                    "squared_error": error**2,
                }
            )
    aggregate = pd.DataFrame(rows)
    per_solute = pd.DataFrame(per_rows)
    aggregate.to_csv(output / "breadth_energy_statistics.csv", index=False)
    per_solute.to_csv(output / "breadth_energy_errors_per_solute.csv", index=False)

    sampled = bootstrap_indices(len(table))
    boot = {}
    boot_rows = []
    for method in METHODS:
        values = per_solute[per_solute.method == method]
        abs_error = values.absolute_error_kcal_mol.to_numpy()
        squared = values.squared_error.to_numpy()
        boot[method] = {
            "MAE": np.mean(abs_error[sampled], axis=1),
            "RMSE": np.sqrt(np.mean(squared[sampled], axis=1)),
        }
        direct = aggregate[aggregate.method == method].iloc[0]
        for metric, distribution in boot[method].items():
            low, high = ci(distribution)
            boot_rows.append(
                {
                    "endpoint": "W4_energy",
                    "quantity": "method_metric",
                    "method": method,
                    "comparator": "",
                    "metric": metric,
                    "point_estimate": direct[metric],
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    for comparator in ("mace_polar_l", best_dipole):
        for metric in ("MAE", "RMSE"):
            delta = boot["glider"][metric] - boot[comparator][metric]
            low, high = ci(delta)
            point = (
                aggregate.loc[aggregate.method == "glider", metric].iloc[0]
                - aggregate.loc[aggregate.method == comparator, metric].iloc[0]
            )
            boot_rows.append(
                {
                    "endpoint": "W4_energy",
                    "quantity": "paired_glider_minus_comparator",
                    "method": "glider",
                    "comparator": comparator,
                    "metric": metric,
                    "point_estimate": point,
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    bootstrap = pd.DataFrame(boot_rows)
    return aggregate, bootstrap


def vector_statistics(
    table: pd.DataFrame, endpoint: str, angular_threshold: float, output: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    for method, group in table.groupby("method", sort=False):
        meaningful = group[
            (group.reference_magnitude >= angular_threshold)
            & np.isfinite(group.angular_error_degrees)
        ]
        rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "n_solutes": len(group),
                "vector_RMSE": float(np.sqrt(np.mean(group.vector_error_norm**2))),
                "magnitude_MAE": float(np.mean(np.abs(group.magnitude_error))),
                "magnitude_RMSE": float(np.sqrt(np.mean(group.magnitude_error**2))),
                "mean_cosine_similarity": float(meaningful.cosine_similarity.mean()),
                "median_angular_error_degrees": float(meaningful.angular_error_degrees.median()),
                "n_directionally_meaningful": len(meaningful),
                "reference_RMS_magnitude": float(np.sqrt(np.mean(group.reference_magnitude**2))),
            }
        )
    aggregate = pd.DataFrame(rows)
    aggregate.to_csv(output / f"breadth_{endpoint}_statistics.csv", index=False)

    solutes = sorted(table.molecule_id.unique())
    sampled = bootstrap_indices(len(solutes))
    method_boot = {}
    boot_rows = []
    for method in METHODS:
        values = table[table.method == method].set_index("molecule_id").loc[solutes]
        method_boot[method] = {
            "vector_RMSE": np.sqrt(np.mean(values.vector_error_norm.to_numpy()[sampled] ** 2, axis=1)),
            "magnitude_MAE": np.mean(np.abs(values.magnitude_error.to_numpy()[sampled]), axis=1),
        }
        direct = aggregate[aggregate.method == method].iloc[0]
        for metric, distribution in method_boot[method].items():
            low, high = ci(distribution)
            boot_rows.append(
                {
                    "endpoint": endpoint,
                    "quantity": "method_metric",
                    "method": method,
                    "comparator": "",
                    "metric": metric,
                    "point_estimate": direct[metric],
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    bootstrap = pd.DataFrame(boot_rows)
    return aggregate, bootstrap, table


def orientation_statistics(
    table: pd.DataFrame, output: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    for probe_id, case_group in table.groupby("probe_id", sort=True):
        reference_table = case_group[case_group.method == "qm"].sort_values("rotation_index")
        reference = reference_table.energy_kcal_mol.to_numpy()
        centered_reference = reference - np.mean(reference)
        reference_amplitude = float(np.ptp(reference))
        rotations = {
            row.rotation_id: Rotation.from_quat(
                [row.quaternion_x, row.quaternion_y, row.quaternion_z, row.quaternion_w]
            )
            for row in reference_table.itertuples()
        }
        reference_minimum = reference_table.iloc[int(np.argmin(reference))].rotation_id
        reference_maximum = reference_table.iloc[int(np.argmax(reference))].rotation_id
        for method, method_group in case_group.groupby("method", sort=False):
            method_group = method_group.sort_values("rotation_index")
            predicted = method_group.energy_kcal_mol.to_numpy()
            error = predicted - reference
            centered_error = predicted - np.mean(predicted) - centered_reference
            pearson, spearman = safe_correlations(reference, predicted)
            predicted_minimum = method_group.iloc[int(np.argmin(predicted))].rotation_id
            predicted_maximum = method_group.iloc[int(np.argmax(predicted))].rotation_id
            if reference_amplitude >= 0.05 and method != "full_frozen_electrostatic":
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
                    "probe_id": probe_id,
                    "case_id": method_group.iloc[0].case_id,
                    "molecule_id": method_group.iloc[0].molecule_id,
                    "family": method_group.iloc[0].family,
                    "method": method,
                    "profile_RMSE_kcal_mol": float(np.sqrt(np.mean(error**2))),
                    "centered_profile_RMSE_kcal_mol": float(np.sqrt(np.mean(centered_error**2))),
                    "profile_MAE_kcal_mol": float(np.mean(np.abs(error))),
                    "profile_bias_kcal_mol": float(np.mean(error)),
                    "Pearson_r": pearson,
                    "Spearman_rho": spearman,
                    "reference_amplitude_kcal_mol": reference_amplitude,
                    "predicted_amplitude_kcal_mol": float(np.ptp(predicted)),
                    "amplitude_error_kcal_mol": float(np.ptp(predicted) - reference_amplitude),
                    "minimum_geodesic_error_degrees": minimum_angle,
                    "maximum_geodesic_error_degrees": maximum_angle,
                }
            )
    per_solute = pd.DataFrame(rows)
    per_solute.to_csv(output / "orientation_statistics_per_solute.csv", index=False)
    aggregate = (
        per_solute.groupby("method", sort=False)
        .agg(
            n_solutes=("probe_id", "size"),
            mean_profile_RMSE_kcal_mol=("profile_RMSE_kcal_mol", "mean"),
            mean_centered_profile_RMSE_kcal_mol=("centered_profile_RMSE_kcal_mol", "mean"),
            mean_absolute_amplitude_error_kcal_mol=("amplitude_error_kcal_mol", lambda x: float(np.mean(np.abs(x)))),
            median_Pearson_r=("Pearson_r", "median"),
            median_Spearman_rho=("Spearman_rho", "median"),
            median_minimum_geodesic_error_degrees=("minimum_geodesic_error_degrees", "median"),
        )
        .reset_index()
    )
    aggregate.to_csv(output / "orientation_aggregate_statistics.csv", index=False)

    solutes = sorted(per_solute.molecule_id.unique())
    sampled = bootstrap_indices(len(solutes))
    boot_rows = []
    for metric_column, metric_label, aggregate_column in (
        (
            "profile_RMSE_kcal_mol",
            "mean_profile_RMSE_kcal_mol",
            "mean_profile_RMSE_kcal_mol",
        ),
        (
            "centered_profile_RMSE_kcal_mol",
            "mean_centered_profile_RMSE_kcal_mol",
            "mean_centered_profile_RMSE_kcal_mol",
        ),
    ):
        for method in METHODS:
            values = (
                per_solute[per_solute.method == method]
                .set_index("molecule_id")
                .loc[solutes, metric_column]
                .to_numpy()
            )
            distribution = np.mean(values[sampled], axis=1)
            low, high = ci(distribution)
            direct = aggregate.loc[aggregate.method == method, aggregate_column].iloc[0]
            boot_rows.append(
                {
                    "endpoint": "orientation",
                    "quantity": "method_metric",
                    "method": method,
                    "comparator": "",
                    "metric": metric_label,
                    "point_estimate": direct,
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    bootstrap = pd.DataFrame(boot_rows)
    return aggregate, bootstrap, per_solute


def add_paired_bootstrap(
    bootstrap: pd.DataFrame,
    per_solute: pd.DataFrame,
    endpoint: str,
    metric_column: str,
    metric_label: str,
    best_dipole: str,
) -> pd.DataFrame:
    solutes = sorted(per_solute.molecule_id.unique())
    sampled = bootstrap_indices(len(solutes))
    distributions = {}
    for method in METHODS:
        values = (
            per_solute[per_solute.method == method]
            .set_index("molecule_id")
            .loc[solutes, metric_column]
            .to_numpy()
        )
        distributions[method] = np.mean(values[sampled], axis=1)
    rows = []
    for comparator in ("mace_polar_l", best_dipole):
        delta = distributions["glider"] - distributions[comparator]
        low, high = ci(delta)
        rows.append(
            {
                "endpoint": endpoint,
                "quantity": "paired_glider_minus_comparator",
                "method": "glider",
                "comparator": comparator,
                "metric": metric_label,
                "point_estimate": float(np.mean(distributions["glider"]) - np.mean(distributions[comparator])),
                "ci95_low": low,
                "ci95_high": high,
            }
        )
    return pd.concat([bootstrap, pd.DataFrame(rows)], ignore_index=True)


def win_fractions(
    energy_errors: pd.DataFrame,
    torque: pd.DataFrame,
    orientation: pd.DataFrame,
    dipole_comparators: dict[str, str],
) -> pd.DataFrame:
    rows = []
    endpoints = {
        "W4_energy_absolute_error": (energy_errors, "absolute_error_kcal_mol"),
        "W4_torque_vector_error": (torque, "vector_error_norm"),
        "W4_orientation_absolute_profile_RMSE": (orientation, "profile_RMSE_kcal_mol"),
        "W4_orientation_centered_profile_RMSE": (orientation, "centered_profile_RMSE_kcal_mol"),
    }
    for endpoint, (table, metric) in endpoints.items():
        pivot = table.pivot(index="molecule_id", columns="method", values=metric)
        for comparator in ("mace_polar_l", dipole_comparators[endpoint]):
            difference = pivot.glider - pivot[comparator]
            rows.append(
                {
                    "endpoint": endpoint,
                    "comparator": comparator,
                    "n_solutes": len(pivot),
                    "n_glider_improved": int(np.sum(difference < 0)),
                    "fraction_glider_improved": float(np.mean(difference < 0)),
                    "median_paired_glider_minus_comparator": float(np.median(difference)),
                }
            )
    return pd.DataFrame(rows)


def range_and_cumulative(
    energy: pd.DataFrame, output: Path, best_dipole: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rank_rows = []
    rank_boot_rows = []
    for rank, group in energy.groupby("water_rank", sort=True):
        group = group.sort_values("molecule_id")
        reference = group.qm_energy_kcal_mol.to_numpy()
        sampled = bootstrap_indices(len(group))
        distributions = {}
        for method in METHODS:
            error = group[f"{method}_energy_kcal_mol"].to_numpy() - reference
            metrics = error_metrics(reference, group[f"{method}_energy_kcal_mol"].to_numpy())
            rank_rows.append(
                {"water_rank": rank, "method": method, "n_solutes": len(group), **metrics}
            )
            distributions[method] = np.mean(np.abs(error)[sampled], axis=1)
            low, high = ci(distributions[method])
            rank_boot_rows.append(
                {
                    "water_rank": rank,
                    "quantity": "method_MAE",
                    "method": method,
                    "comparator": "",
                    "point_estimate": metrics["MAE"],
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
        for comparator in ("mace_polar_l", best_dipole):
            delta = distributions["glider"] - distributions[comparator]
            low, high = ci(delta)
            rank_boot_rows.append(
                {
                    "water_rank": rank,
                    "quantity": "paired_glider_minus_comparator_MAE",
                    "method": "glider",
                    "comparator": comparator,
                    "point_estimate": float(np.mean(delta)),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    rank = pd.DataFrame(rank_rows)
    rank.to_csv(output / "rank_statistics.csv", index=False)
    pd.DataFrame(rank_boot_rows).to_csv(output / "rank_bootstrap_intervals.csv", index=False)

    lower = np.floor(energy.oxygen_distance_A.min() * 2) / 2
    upper = np.ceil(energy.oxygen_distance_A.max() * 2) / 2
    edges = np.arange(lower, upper + 0.5001, 0.5)
    energy = energy.copy()
    energy["distance_bin"] = pd.cut(energy.oxygen_distance_A, edges, right=False, include_lowest=True)
    distance_rows = []
    distance_boot_rows = []
    for distance_bin, group in energy.groupby("distance_bin", observed=True, sort=True):
        reference = group.qm_energy_kcal_mol.to_numpy()
        method_per_solute = {}
        for method in METHODS:
            prediction = group[f"{method}_energy_kcal_mol"].to_numpy()
            error = prediction - reference
            per_solute = pd.DataFrame(
                {"molecule_id": group.molecule_id, "absolute_error": np.abs(error), "squared_error": error**2}
            ).groupby("molecule_id", as_index=False).mean()
            distance_rows.append(
                {
                    "distance_bin": str(distance_bin),
                    "distance_low_A": float(distance_bin.left),
                    "distance_high_A": float(distance_bin.right),
                    "method": method,
                    "n_probes": len(group),
                    "n_solutes": group.molecule_id.nunique(),
                    "equal_solute_MAE": float(per_solute.absolute_error.mean()),
                    "equal_solute_RMSE": float(np.sqrt(per_solute.squared_error.mean())),
                }
            )
            method_per_solute[method] = per_solute.set_index("molecule_id").absolute_error
        common_solutes = sorted(
            set.intersection(*(set(values.index) for values in method_per_solute.values()))
        )
        sampled = bootstrap_indices(len(common_solutes))
        distributions = {}
        for method, values in method_per_solute.items():
            array = values.loc[common_solutes].to_numpy()
            distributions[method] = np.mean(array[sampled], axis=1)
            low, high = ci(distributions[method])
            distance_boot_rows.append(
                {
                    "distance_bin": str(distance_bin),
                    "quantity": "method_equal_solute_MAE",
                    "method": method,
                    "comparator": "",
                    "n_solutes": len(common_solutes),
                    "point_estimate": float(np.mean(array)),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
        for comparator in ("mace_polar_l", best_dipole):
            delta = distributions["glider"] - distributions[comparator]
            low, high = ci(delta)
            distance_boot_rows.append(
                {
                    "distance_bin": str(distance_bin),
                    "quantity": "paired_glider_minus_comparator_MAE",
                    "method": "glider",
                    "comparator": comparator,
                    "n_solutes": len(common_solutes),
                    "point_estimate": float(np.mean(delta)),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    distance = pd.DataFrame(distance_rows)
    distance.to_csv(output / "distance_bin_statistics.csv", index=False)
    pd.DataFrame(distance_boot_rows).to_csv(
        output / "distance_bin_bootstrap_intervals.csv", index=False
    )

    cumulative_rows = []
    for molecule_id, group in energy.groupby("molecule_id", sort=True):
        group = group.sort_values("water_rank")
        for method in METHODS:
            values = group[f"{method}_energy_kcal_mol"].to_numpy()
            for water_rank_value, cumulative_value in zip(
                group.water_rank, np.cumsum(values), strict=True
            ):
                cumulative_rows.append(
                    {
                        "molecule_id": molecule_id,
                        "family": group.family.iloc[0],
                        "through_water_rank": water_rank_value,
                        "method": method,
                        "cumulative_energy_kcal_mol": cumulative_value,
                    }
                )
    cumulative = pd.DataFrame(cumulative_rows)
    reference = cumulative[cumulative.method == "qm"][["molecule_id", "through_water_rank", "cumulative_energy_kcal_mol"]].rename(
        columns={"cumulative_energy_kcal_mol": "qm_cumulative_energy_kcal_mol"}
    )
    cumulative = cumulative.merge(reference, on=["molecule_id", "through_water_rank"])
    cumulative["error_kcal_mol"] = cumulative.cumulative_energy_kcal_mol - cumulative.qm_cumulative_energy_kcal_mol
    cumulative.to_csv(output / "cumulative_outer_environment.csv", index=False)
    cumulative_stats = (
        cumulative.groupby(["through_water_rank", "method"], sort=True)
        .agg(
            n_solutes=("molecule_id", "size"),
            MAE=("error_kcal_mol", lambda x: float(np.mean(np.abs(x)))),
            RMSE=("error_kcal_mol", lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))),
            signed_bias=("error_kcal_mol", "mean"),
        )
        .reset_index()
    )
    cumulative_stats.to_csv(output / "cumulative_statistics.csv", index=False)
    cumulative_boot_rows = []
    for rank_value, group in cumulative.groupby("through_water_rank", sort=True):
        distributions = {}
        sampled = bootstrap_indices(group.molecule_id.nunique())
        for method in METHODS:
            values = group[group.method == method].sort_values("molecule_id")
            errors = values.error_kcal_mol.to_numpy()
            distributions[method] = np.mean(np.abs(errors)[sampled], axis=1)
            low, high = ci(distributions[method])
            cumulative_boot_rows.append(
                {
                    "through_water_rank": rank_value,
                    "quantity": "method_MAE",
                    "method": method,
                    "comparator": "",
                    "point_estimate": float(np.mean(np.abs(errors))),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
        for comparator in ("mace_polar_l", best_dipole):
            delta = distributions["glider"] - distributions[comparator]
            low, high = ci(delta)
            cumulative_boot_rows.append(
                {
                    "through_water_rank": rank_value,
                    "quantity": "paired_glider_minus_comparator_MAE",
                    "method": "glider",
                    "comparator": comparator,
                    "point_estimate": float(np.mean(delta)),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    pd.DataFrame(cumulative_boot_rows).to_csv(
        output / "cumulative_bootstrap_intervals.csv", index=False
    )
    return rank, distance, cumulative_stats


def torque_range_statistics(
    torque: pd.DataFrame, output: Path, best_dipole: str
) -> pd.DataFrame:
    rows = []
    bootstrap_rows = []
    for rank, group in torque.groupby("water_rank", sort=True):
        distributions = {}
        point_estimates = {}
        sampled = bootstrap_indices(group.molecule_id.nunique())
        for method in METHODS:
            values = group[group.method == method].sort_values("molecule_id")
            errors = values.vector_error_norm.to_numpy()
            meaningful = values.reference_magnitude.to_numpy() >= 0.01
            angular = values.angular_error_degrees.to_numpy()[meaningful]
            rows.append(
                {
                    "water_rank": rank,
                    "method": method,
                    "n_solutes": len(values),
                    "vector_RMSE": float(np.sqrt(np.mean(errors**2))),
                    "magnitude_MAE": float(np.mean(np.abs(values.magnitude_error))),
                    "median_angular_error_degrees": float(np.median(angular))
                    if len(angular)
                    else float("nan"),
                    "n_directionally_meaningful": int(np.sum(meaningful)),
                    "reference_RMS_magnitude": float(
                        np.sqrt(np.mean(values.reference_magnitude.to_numpy() ** 2))
                    ),
                }
            )
            point_estimates[method] = float(np.sqrt(np.mean(errors**2)))
            distributions[method] = np.sqrt(np.mean(errors[sampled] ** 2, axis=1))
            low, high = ci(distributions[method])
            bootstrap_rows.append(
                {
                    "water_rank": rank,
                    "quantity": "method_vector_RMSE",
                    "method": method,
                    "comparator": "",
                    "point_estimate": float(np.sqrt(np.mean(errors**2))),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
        for comparator in ("mace_polar_l", best_dipole):
            delta = distributions["glider"] - distributions[comparator]
            low, high = ci(delta)
            bootstrap_rows.append(
                {
                    "water_rank": rank,
                    "quantity": "paired_glider_minus_comparator_vector_RMSE",
                    "method": "glider",
                    "comparator": comparator,
                    "point_estimate": point_estimates["glider"]
                    - point_estimates[comparator],
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output / "torque_rank_statistics.csv", index=False)
    pd.DataFrame(bootstrap_rows).to_csv(
        output / "torque_rank_bootstrap_intervals.csv", index=False
    )
    return result


def directionality_fractions(
    energy: pd.DataFrame, torque: pd.DataFrame, orientation: pd.DataFrame, output: Path
) -> pd.DataFrame:
    rows = []
    for row in energy.itertuples():
        denominator = abs(row.full_frozen_electrostatic_kcal_mol)
        numerator = abs(row.qm_energy_kcal_mol)
        rows.append(
            {
                "molecule_id": row.molecule_id,
                "family": row.family,
                "water_rank": row.water_rank,
                "oxygen_distance_A": row.oxygen_distance_A,
                "quantity": "absolute_energy",
                "response_value": numerator,
                "full_frozen_value": denominator,
                "denominator_threshold": 0.1,
                "ratio": numerator / denominator if denominator >= 0.1 else np.nan,
            }
        )
    torque_pivot = torque.pivot_table(
        index=["molecule_id", "family", "water_rank", "oxygen_distance_A"],
        columns="method",
        values="magnitude",
    ).reset_index()
    for row in torque_pivot.itertuples():
        denominator = row.full_frozen_electrostatic
        numerator = row.qm
        rows.append(
            {
                "molecule_id": row.molecule_id,
                "family": row.family,
                "water_rank": row.water_rank,
                "oxygen_distance_A": row.oxygen_distance_A,
                "quantity": "torque_magnitude",
                "response_value": numerator,
                "full_frozen_value": denominator,
                "denominator_threshold": 0.01,
                "ratio": numerator / denominator if denominator >= 0.01 else np.nan,
            }
        )
    for molecule_id, group in orientation.groupby("molecule_id", sort=True):
        amplitudes = group.groupby("method").energy_kcal_mol.agg(lambda x: float(np.ptp(x)))
        denominator = amplitudes["full_frozen_electrostatic"]
        numerator = amplitudes["qm"]
        rows.append(
            {
                "molecule_id": molecule_id,
                "family": group.family.iloc[0],
                "water_rank": 4,
                "oxygen_distance_A": group.oxygen_distance_A.iloc[0],
                "quantity": "orientation_modulation",
                "response_value": numerator,
                "full_frozen_value": denominator,
                "denominator_threshold": 0.05,
                "ratio": numerator / denominator if denominator >= 0.05 else np.nan,
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(output / "directionality_fractions.csv", index=False)
    summary = (
        result.groupby(["quantity", "water_rank"], dropna=False, sort=True)
        .agg(
            n_total=("molecule_id", "size"),
            n_defined=("ratio", "count"),
            median_ratio=("ratio", "median"),
            mean_ratio=("ratio", "mean"),
            median_response_value=("response_value", "median"),
            median_full_frozen_value=("full_frozen_value", "median"),
        )
        .reset_index()
    )
    summary.to_csv(output / "directionality_fraction_statistics.csv", index=False)
    return result


def sensitivity_tables(
    energy_errors: pd.DataFrame,
    torque: pd.DataFrame,
    orientation: pd.DataFrame,
    output: Path,
) -> None:
    endpoint_tables = {
        "W4_energy_MAE": (energy_errors, "absolute_error_kcal_mol"),
        "W4_torque_vector_RMSE": (torque.assign(metric=torque.vector_error_norm**2), "metric"),
        "W4_orientation_centered_RMSE": (orientation.assign(metric=orientation.centered_profile_RMSE_kcal_mol**2), "metric"),
    }
    family_rows = []
    loo_rows = []
    for endpoint, (table, metric) in endpoint_tables.items():
        for method, method_table in table.groupby("method", sort=False):
            family_values = method_table.groupby("family")[metric].mean()
            value = float(family_values.mean())
            if "RMSE" in endpoint:
                value = float(np.sqrt(value))
            family_rows.append(
                {
                    "endpoint": endpoint,
                    "method": method,
                    "n_families": len(family_values),
                    "equal_family_value": value,
                }
            )
        for held_out in sorted(table.molecule_id.unique()):
            subset = table[table.molecule_id != held_out]
            for method, method_table in subset.groupby("method", sort=False):
                value = float(method_table[metric].mean())
                if "RMSE" in endpoint:
                    value = float(np.sqrt(value))
                loo_rows.append(
                    {
                        "endpoint": endpoint,
                        "held_out_solute": held_out,
                        "method": method,
                        "n_solutes": subset.molecule_id.nunique(),
                        "value": value,
                    }
                )
    pd.DataFrame(family_rows).to_csv(output / "family_sensitivity.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(output / "leave_one_solute_out.csv", index=False)


def save_figure(figure: plt.Figure, stem: Path) -> None:
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def make_figures(
    breadth_energy: pd.DataFrame,
    range_energy: pd.DataFrame,
    energy_errors: pd.DataFrame,
    torque: pd.DataFrame,
    orientation_stats: pd.DataFrame,
    rank: pd.DataFrame,
    cumulative: pd.DataFrame,
    torque_rank: pd.DataFrame,
    fractions: pd.DataFrame,
    selected_origins: dict[str, str],
    figures: Path,
) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    energy_methods = (
        "glider",
        "mace_polar_l",
        selected_origins["W4_energy_MAE"],
    )
    torque_methods = (
        "glider",
        "mace_polar_l",
        selected_origins["W4_torque_vector_RMSE"],
    )
    absolute_orientation_methods = (
        "glider",
        "mace_polar_l",
        selected_origins["W4_orientation_absolute_profile_RMSE"],
    )
    centered_orientation_methods = (
        "glider",
        "mace_polar_l",
        selected_origins["W4_orientation_centered_profile_RMSE"],
    )
    reference = breadth_energy.qm_energy_kcal_mol.to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1))
    values = [reference]
    for method, marker in zip(energy_methods, ("o", "s", "^"), strict=True):
        prediction = breadth_energy[f"{method}_energy_kcal_mol"].to_numpy()
        values.append(prediction)
        axes[0].scatter(reference, prediction, label=LABELS[method], color=COLORS[method], marker=marker, s=43)
    low, high = min(map(np.min, values)), max(map(np.max, values))
    pad = 0.08 * (high - low or 1)
    axes[0].plot([low - pad, high + pad], [low - pad, high + pad], color="0.4", lw=1)
    axes[0].set(xlabel="QM response coupling (kcal mol$^{-1}$)", ylabel="Predicted response coupling (kcal mol$^{-1}$)")
    axes[0].legend(frameon=False, fontsize=8)
    pivot = energy_errors[energy_errors.method.isin(energy_methods)].pivot(index="molecule_id", columns="method", values="absolute_error_kcal_mol")
    for _, row in pivot.iterrows():
        axes[1].plot(range(3), [row[m] for m in energy_methods], color="0.75", lw=0.8)
    for index, method in enumerate(energy_methods):
        axes[1].scatter(np.full(len(pivot), index), pivot[method], color=COLORS[method], s=35, zorder=3)
    axes[1].set_xticks(range(3), ["GLIDER", "MACE-L", "Exact dipole"], rotation=15)
    axes[1].set_ylabel("Absolute error (kcal mol$^{-1}$)")
    fig.tight_layout()
    save_figure(fig, figures / "confirmatory_breadth_energy")

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.1))
    torque_pivot = torque[torque.method.isin(torque_methods)].pivot(index="molecule_id", columns="method", values="vector_error_norm")
    absolute_orient_pivot = orientation_stats[
        orientation_stats.method.isin(absolute_orientation_methods)
    ].pivot(index="molecule_id", columns="method", values="profile_RMSE_kcal_mol")
    centered_orient_pivot = orientation_stats[
        orientation_stats.method.isin(centered_orientation_methods)
    ].pivot(index="molecule_id", columns="method", values="centered_profile_RMSE_kcal_mol")
    for axis, pivot_table, methods, ylabel in (
        (axes[0], torque_pivot, torque_methods, "Torque vector error (kcal mol$^{-1}$)"),
        (axes[1], absolute_orient_pivot, absolute_orientation_methods, "Absolute orientation-profile RMSE (kcal mol$^{-1}$)"),
        (axes[2], centered_orient_pivot, centered_orientation_methods, "Centred orientation-profile RMSE (kcal mol$^{-1}$)"),
    ):
        for _, row in pivot_table.iterrows():
            axis.plot(range(3), [row[m] for m in methods], color="0.75", lw=0.8)
        for index, method in enumerate(methods):
            axis.scatter(np.full(len(pivot_table), index), pivot_table[method], color=COLORS[method], s=35, zorder=3)
        axis.set_xticks(range(3), ["GLIDER", "MACE-L", "Exact dipole"], rotation=15)
        axis.set_ylabel(ylabel)
    fig.tight_layout()
    save_figure(fig, figures / "confirmatory_torque_orientation")

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.1))
    for method in energy_methods:
        values = rank[rank.method == method]
        axes[0].plot(values.water_rank, values.MAE, marker="o", color=COLORS[method], label=LABELS[method])
        raw_error = np.abs(
            range_energy[f"{method}_energy_kcal_mol"]
            - range_energy.qm_energy_kcal_mol
        )
        axes[1].scatter(
            range_energy.oxygen_distance_A,
            raw_error,
            color=COLORS[method],
            alpha=0.45,
            s=18,
            label=LABELS[method],
        )
    axes[0].set(xlabel="External-water rank", ylabel="Energy MAE (kcal mol$^{-1}$)", xticks=range(4, 13))
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set(xlabel="O distance to nearest solute atom (Å)", ylabel="Absolute energy error (kcal mol$^{-1}$)")
    fig.tight_layout()
    save_figure(fig, figures / "spatial_range_rank_distance")

    fig, axis = plt.subplots(figsize=(5.3, 4.1))
    for method in torque_methods:
        values = torque_rank[torque_rank.method == method]
        axis.plot(
            values.water_rank,
            values.vector_RMSE,
            marker="o",
            color=COLORS[method],
            label=LABELS[method],
        )
    axis.set(
        xlabel="External-water rank",
        ylabel="Torque vector RMSE (kcal mol$^{-1}$)",
        xticks=range(4, 13),
    )
    axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save_figure(fig, figures / "spatial_range_torque")

    fig, axis = plt.subplots(figsize=(5.3, 4.1))
    for method in energy_methods:
        values = cumulative[cumulative.method == method]
        axis.plot(values.through_water_rank, values.MAE, marker="o", color=COLORS[method], label=LABELS[method])
    axis.set(xlabel="Cumulative probes W4 through W$k$", ylabel="Cumulative-coupling MAE (kcal mol$^{-1}$)", xticks=range(4, 13))
    axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save_figure(fig, figures / "cumulative_outer_environment")

    w4 = fractions[fractions.water_rank == 4]
    quantities = ["absolute_energy", "torque_magnitude", "orientation_modulation"]
    data = [w4[w4.quantity == quantity].ratio.dropna().to_numpy() for quantity in quantities]
    fig, axis = plt.subplots(figsize=(5.5, 4.2))
    axis.boxplot(data, tick_labels=["Energy offset", "Torque", "Orientation\nmodulation"], showfliers=False)
    for index, values in enumerate(data, 1):
        axis.scatter(np.full(len(values), index), values, color="#1769AA", alpha=0.65, s=28)
    axis.set_ylabel("Response / full frozen electrostatic magnitude")
    fig.tight_layout()
    save_figure(fig, figures / "directionality_fraction")

    # Compact main-text summary. Endpoint-specific dipole origins follow the
    # frozen panel-wide selection rule; no origin is selected case by case.
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 8.0))
    values = [reference]
    for method, marker in zip(energy_methods, ("o", "s", "^"), strict=True):
        prediction = breadth_energy[f"{method}_energy_kcal_mol"].to_numpy()
        values.append(prediction)
        axes[0, 0].scatter(
            reference,
            prediction,
            label=LABELS[method],
            color=COLORS[method],
            marker=marker,
            s=40,
        )
    low, high = min(map(np.min, values)), max(map(np.max, values))
    pad = 0.08 * (high - low or 1)
    axes[0, 0].plot(
        [low - pad, high + pad], [low - pad, high + pad], color="0.4", lw=1
    )
    axes[0, 0].set(
        xlabel="QM response coupling (kcal mol$^{-1}$)",
        ylabel="Predicted response coupling (kcal mol$^{-1}$)",
    )
    axes[0, 0].legend(frameon=False, fontsize=7)

    for axis, pivot_table, methods, ylabel in (
        (axes[0, 1], torque_pivot, torque_methods, "Torque vector error (kcal mol$^{-1}$)"),
        (axes[1, 0], absolute_orient_pivot, absolute_orientation_methods, "Absolute orientation-profile RMSE (kcal mol$^{-1}$)"),
    ):
        for _, row in pivot_table.iterrows():
            axis.plot(range(3), [row[m] for m in methods], color="0.75", lw=0.8)
        for index, method in enumerate(methods):
            axis.scatter(
                np.full(len(pivot_table), index),
                pivot_table[method],
                color=COLORS[method],
                s=32,
                zorder=3,
            )
        axis.set_xticks(
            range(3), ["GLIDER", "MACE-L", "Exact dipole"], rotation=15
        )
        axis.set_ylabel(ylabel)

    for method in energy_methods:
        values = rank[rank.method == method]
        axes[1, 1].plot(
            values.water_rank,
            values.MAE,
            marker="o",
            color=COLORS[method],
            label=LABELS[method],
        )
    axes[1, 1].set(
        xlabel="External-water rank",
        ylabel="Energy MAE (kcal mol$^{-1}$)",
        xticks=range(4, 13),
    )
    for label, axis in zip("abcd", axes.flat, strict=True):
        axis.text(
            -0.15,
            1.05,
            label,
            transform=axis.transAxes,
            fontweight="bold",
            fontsize=12,
            va="top",
        )
    fig.tight_layout()
    save_figure(fig, figures / "downstream_physics_summary")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluated", type=Path, default=RESULT_ROOT / "evaluated")
    parser.add_argument("--output", type=Path, default=RESULT_ROOT / "analysis")
    parser.add_argument("--figures", type=Path, default=RESULT_ROOT / "figures")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    energy = pd.read_csv(args.evaluated / "range_energy_per_probe.csv")
    breadth = pd.read_csv(args.evaluated / "breadth_energy_per_solute.csv")
    force = pd.read_csv(args.evaluated / "breadth_force_per_solute.csv")
    torque = pd.read_csv(args.evaluated / "breadth_torque_per_solute.csv")
    range_torque = pd.read_csv(args.evaluated / "range_torque_per_probe.csv")
    orientation = pd.read_csv(args.evaluated / "orientation_profiles.csv")
    best_energy_dipole = best_dipole_origin(breadth)

    energy_aggregate, energy_boot = breadth_energy(
        breadth, args.output, best_energy_dipole
    )
    energy_errors = pd.read_csv(args.output / "breadth_energy_errors_per_solute.csv")
    force_aggregate, force_boot, force_per = vector_statistics(force, "force", 0.05, args.output)
    torque_aggregate, torque_boot, torque_per = vector_statistics(torque, "torque", 0.01, args.output)
    orientation_aggregate, orientation_boot, orientation_per = orientation_statistics(orientation, args.output)
    best_force_dipole = best_vector_dipole_origin(force_per)
    best_torque_dipole = best_vector_dipole_origin(torque_per)
    best_orientation_absolute_dipole = best_orientation_dipole_origin(
        orientation_per, "profile_RMSE_kcal_mol"
    )
    best_orientation_centered_dipole = best_orientation_dipole_origin(
        orientation_per, "centered_profile_RMSE_kcal_mol"
    )
    selected_origins = {
        "W4_energy_MAE": best_energy_dipole,
        "W4_force_vector_RMSE": best_force_dipole,
        "W4_torque_vector_RMSE": best_torque_dipole,
        "W4_orientation_absolute_profile_RMSE": best_orientation_absolute_dipole,
        "W4_orientation_centered_profile_RMSE": best_orientation_centered_dipole,
    }
    (args.output / "selected_fixed_dipole_origin.json").write_text(
        json.dumps(
            {
                "selection_rule": "one panel-wide fixed origin per endpoint, lowest equal-solute aggregate primary error; never casewise",
                "selected_origins": selected_origins,
                "all_origins_retained": True,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    torque_boot = add_paired_bootstrap(
        torque_boot,
        torque_per,
        "torque",
        "vector_error_norm",
        "vector_RMSE",
        best_torque_dipole,
    )
    orientation_boot = add_paired_bootstrap(
        orientation_boot,
        orientation_per,
        "orientation",
        "profile_RMSE_kcal_mol",
        "mean_profile_RMSE_kcal_mol",
        best_orientation_absolute_dipole,
    )
    orientation_boot = add_paired_bootstrap(
        orientation_boot,
        orientation_per,
        "orientation",
        "centered_profile_RMSE_kcal_mol",
        "mean_centered_profile_RMSE_kcal_mol",
        best_orientation_centered_dipole,
    )
    pd.concat([energy_boot, force_boot, torque_boot, orientation_boot], ignore_index=True).to_csv(
        args.output / "blocked_bootstrap_intervals.csv", index=False
    )
    wins = win_fractions(
        energy_errors,
        torque_per,
        orientation_per,
        {
            "W4_energy_absolute_error": best_energy_dipole,
            "W4_torque_vector_error": best_torque_dipole,
            "W4_orientation_absolute_profile_RMSE": best_orientation_absolute_dipole,
            "W4_orientation_centered_profile_RMSE": best_orientation_centered_dipole,
        },
    )
    wins.to_csv(args.output / "cross_solute_win_fractions.csv", index=False)
    rank, distance, cumulative = range_and_cumulative(
        energy, args.output, best_energy_dipole
    )
    torque_rank = torque_range_statistics(
        range_torque, args.output, best_torque_dipole
    )
    fractions = directionality_fractions(energy, range_torque, orientation, args.output)
    sensitivity_tables(energy_errors, torque_per, orientation_per, args.output)
    make_figures(
        breadth,
        energy,
        energy_errors,
        torque_per,
        orientation_per,
        rank,
        cumulative,
        torque_rank,
        fractions,
        selected_origins,
        args.figures,
    )

    response = breadth.qm_energy_kcal_mol.to_numpy()
    physical = {
        "n_confirmatory_solutes": len(breadth),
        "n_source_families": int(breadth.family.nunique()),
        "response_energy_RMS_kcal_mol": float(np.sqrt(np.mean(response**2))),
        "response_energy_median_absolute_kcal_mol": float(np.median(np.abs(response))),
        "response_energy_min_kcal_mol": float(np.min(response)),
        "response_energy_max_kcal_mol": float(np.max(response)),
        "kBT_298.15K_kcal_mol": 0.592484949,
        "response_RMS_over_kBT": float(np.sqrt(np.mean(response**2)) / 0.592484949),
        "best_fixed_dipole_origin_for_W4_energy": best_energy_dipole,
    }
    (args.output / "physical_scale.json").write_text(
        json.dumps(physical, indent=2, sort_keys=True) + "\n"
    )
    manifest = {
        "protocol": "results/posthoc_downstream_response_confirmatory/protocol_freeze.json",
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_resamples": N_BOOTSTRAP,
        "selected_fixed_dipole_origins": selected_origins,
        "outputs": {
            path.name: sha256(path)
            for path in sorted(args.output.iterdir())
            if path.is_file() and path.name != "analysis_manifest.json"
        },
        "figures": {
            path.name: sha256(path) for path in sorted(args.figures.iterdir()) if path.is_file()
        },
        "prospective_claim": False,
    }
    (args.output / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(physical, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
