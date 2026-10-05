#!/usr/bin/env python3
"""Analyze the frozen exact-QM global-multipole hierarchy control."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Embed TrueType outlines in publication PDFs rather than Matplotlib's Type 3 fonts.
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from multipole_control import CONTROL_ROOT, RESULT_ROOT  # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

BOOTSTRAP_SEED = 2026081902
N_BOOTSTRAP = 100_000
HIERARCHIES = (
    "exact_dipole",
    "exact_dipole_quadrupole",
    "exact_dipole_quadrupole_octupole",
)
EXISTING_METHODS = ("glider", "mace_polar_l", "zero")
LABELS = {
    "glider": "GLIDER",
    "mace_polar_l": "MACE-POLAR-1-L",
    "zero": "Zero response",
    "exact_dipole": "Exact dipole",
    "exact_dipole_quadrupole": "Exact dipole + quadrupole",
    "exact_dipole_quadrupole_octupole": "Exact dipole + quadrupole + octupole",
}
COLORS = {
    "glider": "#1769AA",
    "mace_polar_l": "#D97706",
    "zero": "#777777",
    "exact_dipole": "#218739",
    "exact_dipole_quadrupole": "#8B5A2B",
    "exact_dipole_quadrupole_octupole": "#7A3E9D",
}


def correlations(reference: np.ndarray, prediction: np.ndarray) -> tuple[float, float]:
    if np.std(reference) == 0 or np.std(prediction) == 0:
        return float("nan"), float("nan")
    return (
        float(pearsonr(reference, prediction).statistic),
        float(spearmanr(reference, prediction).statistic),
    )


def ci(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, (0.025, 0.975))
    return float(low), float(high)


def bootstrap_indices(n: int) -> np.ndarray:
    rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    return rng.integers(0, n, size=(N_BOOTSTRAP, n))


def choose_origins(
    energy: pd.DataFrame, torque: pd.DataFrame, orientation: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for hierarchy in HIERARCHIES:
        values = energy[energy.hierarchy == hierarchy]
        scores = values.groupby("origin").absolute_error_kcal_mol.mean()
        rows.append(
            {
                "endpoint": "W4_energy_MAE",
                "hierarchy": hierarchy,
                "selected_origin": scores.idxmin(),
                "selection_score": scores.min(),
            }
        )
        values = torque[torque.hierarchy == hierarchy]
        scores = values.groupby("origin").vector_error_norm.apply(
            lambda x: float(np.sqrt(np.mean(x**2)))
        )
        rows.append(
            {
                "endpoint": "W4_torque_vector_RMSE",
                "hierarchy": hierarchy,
                "selected_origin": scores.idxmin(),
                "selection_score": scores.min(),
            }
        )
        values = orientation[orientation.hierarchy == hierarchy]
        per_solute = (
            values.groupby(["origin", "molecule_id"])
            .profile_squared_error.mean()
            .pow(0.5)
        )
        scores = per_solute.groupby("origin").mean()
        rows.append(
            {
                "endpoint": "W4_orientation_absolute_profile_RMSE",
                "hierarchy": hierarchy,
                "selected_origin": scores.idxmin(),
                "selection_score": scores.min(),
            }
        )
    return pd.DataFrame(rows)


def selected_origin(
    selection: pd.DataFrame, endpoint: str, hierarchy: str
) -> str:
    return str(
        selection[
            (selection.endpoint == endpoint) & (selection.hierarchy == hierarchy)
        ].selected_origin.iloc[0]
    )


def assemble_energy(
    multipoles: pd.DataFrame, existing: pd.DataFrame, selection: pd.DataFrame
) -> pd.DataFrame:
    metadata = [
        "probe_id",
        "case_id",
        "molecule_id",
        "family",
        "water_rank",
        "source_water_index",
        "oxygen_distance_A",
    ]
    rows = []
    for row in existing.itertuples():
        for method in EXISTING_METHODS:
            prediction = float(getattr(row, f"{method}_energy_kcal_mol"))
            reference = float(row.qm_energy_kcal_mol)
            rows.append(
                {
                    **{column: getattr(row, column) for column in metadata},
                    "method": method,
                    "origin": "",
                    "reference_energy_kcal_mol": reference,
                    "predicted_energy_kcal_mol": prediction,
                    "error_kcal_mol": prediction - reference,
                    "absolute_error_kcal_mol": abs(prediction - reference),
                }
            )
    references = existing.set_index("probe_id").qm_energy_kcal_mol
    for hierarchy in HIERARCHIES:
        origin = selected_origin(selection, "W4_energy_MAE", hierarchy)
        values = multipoles[
            (multipoles.hierarchy == hierarchy) & (multipoles.origin == origin)
        ]
        for row in values.itertuples():
            reference = float(references.loc[row.probe_id])
            prediction = float(row.energy_kcal_mol)
            rows.append(
                {
                    **{column: getattr(row, column) for column in metadata},
                    "method": hierarchy,
                    "origin": origin,
                    "reference_energy_kcal_mol": reference,
                    "predicted_energy_kcal_mol": prediction,
                    "error_kcal_mol": prediction - reference,
                    "absolute_error_kcal_mol": abs(prediction - reference),
                }
            )
    return pd.DataFrame(rows)


def energy_statistics(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, group in table.groupby("method", sort=False):
        reference = group.reference_energy_kcal_mol.to_numpy()
        prediction = group.predicted_energy_kcal_mol.to_numpy()
        error = prediction - reference
        pearson, spearman = correlations(reference, prediction)
        rows.append(
            {
                "method": method,
                "n_solutes": group.molecule_id.nunique(),
                "MAE_kcal_mol": float(np.mean(np.abs(error))),
                "RMSE_kcal_mol": float(np.sqrt(np.mean(error**2))),
                "signed_bias_kcal_mol": float(np.mean(error)),
                "Pearson_r": pearson,
                "Spearman_rho": spearman,
            }
        )
    return pd.DataFrame(rows)


def assemble_torque(
    multipoles: pd.DataFrame, existing: pd.DataFrame, selection: pd.DataFrame
) -> pd.DataFrame:
    keep = existing[existing.method.isin(EXISTING_METHODS)].copy()
    rows = [keep]
    for hierarchy in HIERARCHIES:
        origin = selected_origin(selection, "W4_torque_vector_RMSE", hierarchy)
        values = multipoles[
            (multipoles.hierarchy == hierarchy) & (multipoles.origin == origin)
        ].copy()
        values["method"] = hierarchy
        rows.append(values[keep.columns])
    return pd.concat(rows, ignore_index=True)


def torque_statistics(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, group in table.groupby("method", sort=False):
        meaningful = group[
            (group.reference_magnitude >= 0.01)
            & np.isfinite(group.angular_error_degrees)
        ]
        rows.append(
            {
                "method": method,
                "n_solutes": group.molecule_id.nunique(),
                "vector_RMSE_kcal_mol": float(
                    np.sqrt(np.mean(group.vector_error_norm**2))
                ),
                "magnitude_MAE_kcal_mol": float(
                    np.mean(np.abs(group.magnitude_error))
                ),
                "median_angular_error_degrees": float(
                    meaningful.angular_error_degrees.median()
                ),
                "n_directionally_meaningful": len(meaningful),
            }
        )
    return pd.DataFrame(rows)


def orientation_per_solute(
    multipoles: pd.DataFrame, existing: pd.DataFrame, selection: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reference = existing[existing.method == "qm"][
        ["probe_id", "rotation_id", "energy_kcal_mol"]
    ].rename(columns={"energy_kcal_mol": "reference_energy_kcal_mol"})
    profiles = [existing[existing.method.isin(EXISTING_METHODS)].merge(reference)]
    for hierarchy in HIERARCHIES:
        origin = selected_origin(
            selection, "W4_orientation_absolute_profile_RMSE", hierarchy
        )
        values = multipoles[
            (multipoles.hierarchy == hierarchy) & (multipoles.origin == origin)
        ].copy()
        values["method"] = hierarchy
        profiles.append(values)
    profile_table = pd.concat(profiles, ignore_index=True, sort=False)
    rows = []
    for (method, molecule_id), group in profile_table.groupby(
        ["method", "molecule_id"], sort=False
    ):
        prediction = group.energy_kcal_mol.to_numpy()
        reference_values = group.reference_energy_kcal_mol.to_numpy()
        error = prediction - reference_values
        centered_error = (prediction - prediction.mean()) - (
            reference_values - reference_values.mean()
        )
        pearson, spearman = correlations(reference_values, prediction)
        rows.append(
            {
                "method": method,
                "molecule_id": molecule_id,
                "family": group.family.iloc[0],
                "profile_RMSE_kcal_mol": float(np.sqrt(np.mean(error**2))),
                "centered_profile_RMSE_kcal_mol": float(
                    np.sqrt(np.mean(centered_error**2))
                ),
                "profile_MAE_kcal_mol": float(np.mean(np.abs(error))),
                "profile_bias_kcal_mol": float(np.mean(error)),
                "Pearson_r": pearson,
                "Spearman_rho": spearman,
                "reference_amplitude_kcal_mol": float(np.ptp(reference_values)),
                "predicted_amplitude_kcal_mol": float(np.ptp(prediction)),
                "amplitude_error_kcal_mol": float(
                    np.ptp(prediction) - np.ptp(reference_values)
                ),
            }
        )
    return pd.DataFrame(rows), profile_table


def orientation_statistics(table: pd.DataFrame) -> pd.DataFrame:
    return (
        table.groupby("method", sort=False)
        .agg(
            n_solutes=("molecule_id", "size"),
            mean_profile_RMSE_kcal_mol=("profile_RMSE_kcal_mol", "mean"),
            mean_centered_profile_RMSE_kcal_mol=(
                "centered_profile_RMSE_kcal_mol",
                "mean",
            ),
            mean_absolute_amplitude_error_kcal_mol=(
                "amplitude_error_kcal_mol",
                lambda x: float(np.mean(np.abs(x))),
            ),
            median_Pearson_r=("Pearson_r", "median"),
            median_Spearman_rho=("Spearman_rho", "median"),
        )
        .reset_index()
    )


def blocked_bootstrap(
    energy: pd.DataFrame, torque: pd.DataFrame, orientation: pd.DataFrame
) -> pd.DataFrame:
    methods = EXISTING_METHODS + HIERARCHIES
    solutes = sorted(energy.molecule_id.unique())
    sampled = bootstrap_indices(len(solutes))
    endpoint_values: dict[str, dict[str, np.ndarray]] = {
        "W4_energy_MAE": {},
        "W4_energy_RMSE": {},
        "W4_torque_vector_RMSE": {},
        "W4_orientation_absolute_profile_RMSE": {},
        "W4_orientation_centered_profile_RMSE": {},
    }
    point: dict[str, dict[str, float]] = {key: {} for key in endpoint_values}
    for method in methods:
        e = energy[energy.method == method].set_index("molecule_id").loc[solutes]
        absolute = e.absolute_error_kcal_mol.to_numpy()
        squared = e.error_kcal_mol.to_numpy() ** 2
        endpoint_values["W4_energy_MAE"][method] = np.mean(
            absolute[sampled], axis=1
        )
        endpoint_values["W4_energy_RMSE"][method] = np.sqrt(
            np.mean(squared[sampled], axis=1)
        )
        point["W4_energy_MAE"][method] = float(np.mean(absolute))
        point["W4_energy_RMSE"][method] = float(np.sqrt(np.mean(squared)))
        t = torque[torque.method == method].set_index("molecule_id").loc[solutes]
        squared = t.vector_error_norm.to_numpy() ** 2
        endpoint_values["W4_torque_vector_RMSE"][method] = np.sqrt(
            np.mean(squared[sampled], axis=1)
        )
        point["W4_torque_vector_RMSE"][method] = float(np.sqrt(np.mean(squared)))
        o = orientation[orientation.method == method].set_index("molecule_id").loc[
            solutes
        ]
        for endpoint, column in (
            ("W4_orientation_absolute_profile_RMSE", "profile_RMSE_kcal_mol"),
            (
                "W4_orientation_centered_profile_RMSE",
                "centered_profile_RMSE_kcal_mol",
            ),
        ):
            values = o[column].to_numpy()
            endpoint_values[endpoint][method] = np.mean(values[sampled], axis=1)
            point[endpoint][method] = float(np.mean(values))
    rows = []
    for endpoint, method_values in endpoint_values.items():
        for method, distribution in method_values.items():
            low, high = ci(distribution)
            rows.append(
                {
                    "endpoint": endpoint,
                    "quantity": "method_metric",
                    "method": method,
                    "comparator": "",
                    "point_estimate": point[endpoint][method],
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
        for comparator in (
            "mace_polar_l",
            "exact_dipole",
            "exact_dipole_quadrupole",
            "exact_dipole_quadrupole_octupole",
        ):
            delta = method_values["glider"] - method_values[comparator]
            low, high = ci(delta)
            rows.append(
                {
                    "endpoint": endpoint,
                    "quantity": "paired_glider_minus_comparator",
                    "method": "glider",
                    "comparator": comparator,
                    "point_estimate": point[endpoint]["glider"]
                    - point[endpoint][comparator],
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    return pd.DataFrame(rows)


def wins(
    energy: pd.DataFrame, torque: pd.DataFrame, orientation: pd.DataFrame
) -> pd.DataFrame:
    endpoint_data = {
        "W4_energy_absolute_error": (energy, "absolute_error_kcal_mol"),
        "W4_torque_vector_error": (torque, "vector_error_norm"),
        "W4_orientation_absolute_profile_RMSE": (
            orientation,
            "profile_RMSE_kcal_mol",
        ),
        "W4_orientation_centered_profile_RMSE": (
            orientation,
            "centered_profile_RMSE_kcal_mol",
        ),
    }
    rows = []
    for endpoint, (table, metric) in endpoint_data.items():
        pivot = table.pivot(index="molecule_id", columns="method", values=metric)
        for comparator in (
            "mace_polar_l",
            "exact_dipole",
            "exact_dipole_quadrupole",
            "exact_dipole_quadrupole_octupole",
        ):
            difference = pivot.glider - pivot[comparator]
            rows.append(
                {
                    "endpoint": endpoint,
                    "comparator": comparator,
                    "n_solutes": len(pivot),
                    "n_glider_improved": int(np.sum(difference < 0)),
                    "fraction_glider_improved": float(np.mean(difference < 0)),
                    "median_paired_glider_minus_comparator": float(
                        np.median(difference)
                    ),
                }
            )
    return pd.DataFrame(rows)


def sensitivity(
    energy: pd.DataFrame, torque: pd.DataFrame, orientation: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    endpoints = {
        "W4_energy_MAE": (energy, "absolute_error_kcal_mol", False),
        "W4_torque_vector_RMSE": (torque, "vector_error_norm", True),
        "W4_orientation_absolute_profile_RMSE": (
            orientation,
            "profile_RMSE_kcal_mol",
            False,
        ),
        "W4_orientation_centered_profile_RMSE": (
            orientation,
            "centered_profile_RMSE_kcal_mol",
            False,
        ),
    }
    family_rows = []
    loo_rows = []
    for endpoint, (table, column, square_first) in endpoints.items():
        values = table.copy()
        values["metric"] = values[column] ** 2 if square_first else values[column]
        for method, group in values.groupby("method", sort=False):
            family_value = group.groupby("family").metric.mean().mean()
            family_rows.append(
                {
                    "endpoint": endpoint,
                    "method": method,
                    "n_families": group.family.nunique(),
                    "equal_family_value": float(
                        np.sqrt(family_value) if square_first else family_value
                    ),
                }
            )
        for held_out in sorted(values.molecule_id.unique()):
            subset = values[values.molecule_id != held_out]
            for method, group in subset.groupby("method", sort=False):
                value = group.metric.mean()
                loo_rows.append(
                    {
                        "endpoint": endpoint,
                        "held_out_solute": held_out,
                        "method": method,
                        "n_solutes": group.molecule_id.nunique(),
                        "value": float(np.sqrt(value) if square_first else value),
                    }
                )
    return pd.DataFrame(family_rows), pd.DataFrame(loo_rows)


def range_statistics(
    multipoles: pd.DataFrame, existing: pd.DataFrame, selection: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    table = assemble_energy(multipoles, existing, selection)
    rows = []
    boot_rows = []
    for rank, rank_table in table.groupby("water_rank", sort=True):
        solutes = sorted(rank_table.molecule_id.unique())
        sampled = bootstrap_indices(len(solutes))
        for method, group in rank_table.groupby("method", sort=False):
            values = group.set_index("molecule_id").loc[solutes]
            absolute = values.absolute_error_kcal_mol.to_numpy()
            distribution = np.mean(absolute[sampled], axis=1)
            low, high = ci(distribution)
            rows.append(
                {
                    "water_rank": rank,
                    "method": method,
                    "n_solutes": len(solutes),
                    "MAE_kcal_mol": float(np.mean(absolute)),
                    "RMSE_kcal_mol": float(
                        np.sqrt(np.mean(values.error_kcal_mol.to_numpy() ** 2))
                    ),
                }
            )
            boot_rows.append(
                {
                    "water_rank": rank,
                    "method": method,
                    "point_estimate_MAE_kcal_mol": float(np.mean(absolute)),
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    return table, pd.DataFrame(rows), pd.DataFrame(boot_rows)


def extent_statistics(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (rank, origin), group in table.groupby(["water_rank", "origin"]):
        rows.append(
            {
                "water_rank": rank,
                "origin": origin,
                "n_probes": len(group),
                "fraction_oxygen_outside_max_nuclear_radius": float(
                    np.mean(group.oxygen_over_max_base_nuclear_radius > 1)
                ),
                "fraction_oxygen_outside_r90": float(
                    np.mean(group.oxygen_over_r90 > 1)
                ),
                "fraction_oxygen_outside_r95": float(
                    np.mean(group.oxygen_over_r95 > 1)
                ),
                "fraction_oxygen_outside_r99": float(
                    np.mean(group.oxygen_over_r99 > 1)
                ),
                "median_oxygen_over_r99": float(group.oxygen_over_r99.median()),
                "minimum_oxygen_over_r99": float(group.oxygen_over_r99.min()),
                "maximum_oxygen_over_r99": float(group.oxygen_over_r99.max()),
            }
        )
    return pd.DataFrame(rows)


def save_figure(
    energy: pd.DataFrame,
    torque: pd.DataFrame,
    orientation: pd.DataFrame,
    rank: pd.DataFrame,
    output: Path,
) -> None:
    methods = (
        "glider",
        "mace_polar_l",
        "exact_dipole",
        "exact_dipole_quadrupole",
        "exact_dipole_quadrupole_octupole",
    )
    energy_values = energy_statistics(energy).set_index("method")
    torque_values = torque_statistics(torque).set_index("method")
    orientation_values = orientation_statistics(orientation).set_index("method")
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.8))
    panels = (
        (
            axes[0, 0],
            [energy_values.loc[method, "MAE_kcal_mol"] for method in methods],
            "W4 energy MAE\n(kcal mol$^{-1}$)",
        ),
        (
            axes[0, 1],
            [torque_values.loc[method, "vector_RMSE_kcal_mol"] for method in methods],
            "W4 torque vector RMSE\n(kcal mol$^{-1}$)",
        ),
        (
            axes[1, 0],
            [
                orientation_values.loc[
                    method, "mean_profile_RMSE_kcal_mol"
                ]
                for method in methods
            ],
            "W4 orientation-profile RMSE\n(kcal mol$^{-1}$)",
        ),
    )
    short_labels = [
        "GLIDER",
        "MACE-L",
        r"$\mu$",
        r"$\mu$+$\Theta$",
        r"$\mu$+$\Theta$+$\Omega$",
    ]
    for axis, values, ylabel in panels:
        axis.bar(
            range(len(methods)),
            values,
            color=[COLORS[method] for method in methods],
        )
        axis.set_xticks(range(len(methods)), short_labels)
        axis.set_ylabel(ylabel)
    for method in methods:
        values = rank[rank.method == method]
        axes[1, 1].plot(
            values.water_rank,
            values.MAE_kcal_mol,
            marker="o",
            color=COLORS[method],
            label=LABELS[method],
        )
    axes[1, 1].set(
        xlabel="External-water rank",
        ylabel="Energy MAE (kcal mol$^{-1}$)",
        xticks=range(4, 13),
    )
    axes[1, 1].legend(frameon=False, fontsize=7)
    for label, axis in zip("abcd", axes.flat, strict=True):
        axis.text(
            -0.14,
            1.04,
            label,
            transform=axis.transAxes,
            fontweight="bold",
            fontsize=12,
            va="top",
        )
    fig.tight_layout()
    fig.savefig(output / "multipole_hierarchy_control.pdf", bbox_inches="tight")
    fig.savefig(output / "multipole_hierarchy_control.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=CONTROL_ROOT)
    parser.add_argument("--output", type=Path, default=CONTROL_ROOT / "analysis")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    existing_energy = pd.read_csv(RESULT_ROOT / "evaluated/breadth_energy_per_solute.csv")
    range_energy = pd.read_csv(RESULT_ROOT / "evaluated/range_energy_per_probe.csv")
    existing_torque = pd.read_csv(RESULT_ROOT / "evaluated/breadth_torque_per_solute.csv")
    existing_orientation = pd.read_csv(RESULT_ROOT / "evaluated/orientation_profiles.csv")
    raw_energy = pd.read_csv(args.input / "W4_energy_all_origins.csv")
    references = existing_energy.set_index("probe_id").qm_energy_kcal_mol
    raw_energy["reference_energy_kcal_mol"] = raw_energy.probe_id.map(references)
    raw_energy["absolute_error_kcal_mol"] = np.abs(
        raw_energy.energy_kcal_mol - raw_energy.reference_energy_kcal_mol
    )
    raw_torque = pd.read_csv(args.input / "W4_torque_all_origins.csv")
    raw_orientation = pd.read_csv(
        args.input / "W4_orientation_profiles_all_origins.csv"
    )
    raw_orientation["profile_squared_error"] = (
        raw_orientation.energy_kcal_mol
        - raw_orientation.reference_energy_kcal_mol
    ) ** 2
    selection = choose_origins(raw_energy, raw_torque, raw_orientation)
    selection.to_csv(args.output / "selected_panel_wide_origins.csv", index=False)

    energy = assemble_energy(raw_energy, existing_energy, selection)
    energy.to_csv(args.output / "W4_energy_per_solute.csv", index=False)
    energy_stats = energy_statistics(energy)
    energy_stats.to_csv(args.output / "W4_energy_statistics.csv", index=False)
    torque = assemble_torque(raw_torque, existing_torque, selection)
    torque.to_csv(args.output / "W4_torque_per_solute.csv", index=False)
    torque_stats = torque_statistics(torque)
    torque_stats.to_csv(args.output / "W4_torque_statistics.csv", index=False)
    orientation, orientation_profiles = orientation_per_solute(
        raw_orientation, existing_orientation, selection
    )
    orientation.to_csv(args.output / "W4_orientation_per_solute.csv", index=False)
    orientation_profiles.to_csv(
        args.output / "W4_orientation_profiles_selected.csv", index=False
    )
    orientation_stats = orientation_statistics(orientation)
    orientation_stats.to_csv(
        args.output / "W4_orientation_statistics.csv", index=False
    )

    bootstrap = blocked_bootstrap(energy, torque, orientation)
    bootstrap.to_csv(args.output / "blocked_bootstrap_intervals.csv", index=False)
    wins(energy, torque, orientation).to_csv(
        args.output / "cross_solute_win_fractions.csv", index=False
    )
    family, loo = sensitivity(energy, torque, orientation)
    family.to_csv(args.output / "equal_family_sensitivity.csv", index=False)
    loo.to_csv(args.output / "leave_one_solute_out.csv", index=False)

    raw_range = pd.read_csv(args.input / "W4_W12_energy_all_origins.csv")
    range_table, rank_stats, rank_bootstrap = range_statistics(
        raw_range, range_energy, selection
    )
    range_table.to_csv(args.output / "W4_W12_energy_per_probe.csv", index=False)
    rank_stats.to_csv(args.output / "W4_W12_rank_statistics.csv", index=False)
    rank_bootstrap.to_csv(
        args.output / "W4_W12_rank_bootstrap_intervals.csv", index=False
    )

    extent = pd.read_csv(args.input / "source_extent_geometry.csv")
    extent_stats = extent_statistics(extent)
    extent_stats.to_csv(args.output / "source_extent_statistics.csv", index=False)
    save_figure(energy, torque, orientation, rank_stats, args.output)

    summary = {
        "bootstrap": {"seed": BOOTSTRAP_SEED, "resamples": N_BOOTSTRAP},
        "selected_panel_wide_origins": selection.to_dict(orient="records"),
        "W4_energy": energy_stats.to_dict(orient="records"),
        "W4_torque": torque_stats.to_dict(orient="records"),
        "W4_orientation": orientation_stats.to_dict(orient="records"),
        "W4_source_extent": extent_stats[extent_stats.water_rank == 4].to_dict(
            orient="records"
        ),
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
