#!/usr/bin/env python3
"""Check and visualize the post hoc cyclic-carbamate QM separation references."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DISTANCES = (3, 4, 5, 6, 8, 10, 15, 20, 50, 100)
SYSTEMS = ("single_water", "whole_environment")
BLUE = "#0072B2"
AMBER = "#B66B21"
INK = "#202124"
MUTED = "#687078"
GRID = "#D8DDE2"


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def validate(directory: Path) -> list[dict[str, float | str]]:
    with (directory / "summary.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    expected = {(system, distance) for system in SYSTEMS for distance in DISTANCES}
    actual = {(row["system"], int(row["distance_A"])) for row in rows}
    if len(rows) != 20 or actual != expected:
        raise ValueError("Expected the complete, unique 20-condition scan")
    for row in rows:
        config = (f"dissociation__dev_cyclic_carbamate__{row['system']}__"
                  f"{row['distance_A']}A")
        record = json.loads((directory / f"{config}.json").read_text())
        if record["config_id"] != config or record["n_probe_points"] != 512:
            raise ValueError(f"Incorrect record for {config}")
        arrays = directory / record["arrays_file"]
        if digest(arrays) != record["arrays_sha256"]:
            raise ValueError(f"Array digest mismatch: {config}")
        source = ROOT / record["geometry_source"]
        for path, key in ((source, "geometry_source_sha256"),
                          (source.parent / "solute_probe_points.npz", "probe_source_sha256"),
                          (source.parent / "predicted_probe_potentials.npz",
                           "prediction_source_sha256")):
            if digest(path) != record[key]:
                raise ValueError(f"Input digest mismatch: {config}, {path}")
        for component in ("full_complex", "solute_with_neighbour_ghosts",
                          "neighbour_with_solute_ghosts"):
            attempts = record["attempts"][component]
            if not attempts or not attempts[-1]["converged"]:
                raise ValueError(f"Unconverged {component}: {config}")
        for key in ("qm_response_rms_mEh_per_e", "glider_rms_mEh_per_e",
                    "prior_rms_mEh_per_e", "glider_error_rms_mEh_per_e",
                    "prior_error_rms_mEh_per_e"):
            row[key] = float(row[key])
            if not np.isclose(row[key], record[key], rtol=0, atol=1e-12):
                raise ValueError(f"Summary mismatch: {config}, {key}")
        row["distance_A"] = int(row["distance_A"])
    return rows


def plot(rows: list[dict], output: Path) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.4,
        "axes.labelsize": 8.5, "axes.titlesize": 9,
        "xtick.labelsize": 7.7, "ytick.labelsize": 7.7,
        "axes.edgecolor": INK, "axes.labelcolor": INK,
        "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "svg.hashsalt": "cyclic_carbamate_separation_qm_v1",
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.8, 5.2), facecolor="white",
                             constrained_layout=True)
    top_series = (
        ("QM reference", "qm_response_rms_mEh_per_e", INK),
        ("GLIDER", "glider_rms_mEh_per_e", BLUE),
        ("Frozen prior", "prior_rms_mEh_per_e", AMBER),
    )
    bottom_series = (
        ("GLIDER error", "glider_error_rms_mEh_per_e", BLUE),
        ("Prior error", "prior_error_rms_mEh_per_e", AMBER),
    )
    for column, system in enumerate(SYSTEMS):
        group = sorted((row for row in rows if row["system"] == system),
                       key=lambda row: row["distance_A"])
        x = [row["distance_A"] for row in group]
        top, bottom = axes[:, column]
        for label, key, color in top_series:
            top.plot(x, [row[key] for row in group], marker="o", markersize=3.5,
                     linewidth=1.55, color=color, label=label)
        for label, key, color in bottom_series:
            bottom.plot(x, [row[key] for row in group], marker="o", markersize=3.5,
                        linewidth=1.55, color=color, label=label)
        for ax in (top, bottom):
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(2.7, 112)
            ax.set_xticks([3, 5, 10, 20, 50, 100],
                          labels=["3", "5", "10", "20", "50", "100"])
            ax.minorticks_off()
            ax.grid(axis="y", color=GRID, linewidth=0.55)
            ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False)
        top.set_title(("a" if column == 0 else "b") + "  " +
                      ("One water moved" if column == 0 else "Four-water environment moved"),
                      loc="left", weight="bold")
        bottom.set_title(("c" if column == 0 else "d") +
                         "  Error against QM", loc="left", weight="bold")
        bottom.set_xlabel("Minimum interfragment distance (Å)")
        if column == 0:
            top.set_ylabel("Response-potential RMS (mEh/e)")
            bottom.set_ylabel("Error RMS (mEh/e)")
        if column == 0:
            top.legend(frameon=False, fontsize=7.3, ncol=1, loc="lower left")
        bottom.set_ylim(0.04, 10)
    fig.savefig(output / "diagnostic.svg", facecolor="white")
    fig.savefig(output / "diagnostic.pdf", facecolor="white")
    fig.savefig(output / "diagnostic.png", dpi=240, facecolor="white")
    plt.close(fig)


def cpu_comparison(output: Path, cpu_directories: list[Path]) -> str | None:
    comparisons = {}
    for directory in cpu_directories:
        for path in directory.glob("dissociation__*.npz"):
            gpu_path = output / path.name
            if not gpu_path.exists():
                continue
            with np.load(path, allow_pickle=False) as cpu, np.load(
                gpu_path, allow_pickle=False
            ) as gpu:
                difference = (cpu["qm_response_esp_hartree_per_e"] -
                              gpu["qm_response_esp_hartree_per_e"])
                comparisons[path.name] = (
                    float(np.max(np.abs(difference))) * 1000,
                    float(np.sqrt(np.mean(difference**2))) * 1000,
                )
    if not comparisons:
        return None
    values = list(comparisons.values())
    return (f"Independent CPU calculations cross-check {len(values)} distinct "
            f"conditions of the "
            f"20 GPU references. Across their 512-point fields, the largest "
            f"pointwise CPU–GPU difference is {max(x[0] for x in values):.5f} "
            f"mEh/e and the largest per-case RMS difference is "
            f"{max(x[1] for x in values):.5f} mEh/e.")


def note(rows: list[dict], output: Path, cpu_directories: list[Path]) -> None:
    def get(system: str, distance: int) -> dict:
        return next(row for row in rows if row["system"] == system
                    and row["distance_A"] == distance)

    lines = [
        "# QM separation check: dev_cyclic_carbamate",
        "",
        "This post hoc diagnostic evaluates the archived geometries at 3, 4, 5, 6, 8, "
        "10, 15, 20, 50 and 100 Å. One series moves a single water; the other "
        "moves the intact four-water environment. The solute and its 512-point "
        "probe surface stay fixed.",
        "",
        "The reference is counterpoise-consistent complex-minus-fragments "
        "DF-RKS ωB97X-D3(BJ)/def2-TZVPD, grid level 4, SCF tolerance 1e-10. "
        "All 60 component SCFs converged. No predicted values were refitted.",
        "",
        "| System | Distance | QM response RMS | GLIDER RMS | Prior RMS | "
        "GLIDER error | Prior error |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for system in SYSTEMS:
        for distance in (3, 20, 50, 100):
            row = get(system, distance)
            lines.append("| {} | {} Å | {:.4f} | {:.4f} | {:.4f} | {:.4f} | {:.4f} |".format(
                "One water" if system == "single_water" else "Four waters",
                distance, row["qm_response_rms_mEh_per_e"],
                row["glider_rms_mEh_per_e"], row["prior_rms_mEh_per_e"],
                row["glider_error_rms_mEh_per_e"],
                row["prior_error_rms_mEh_per_e"]))
    lines += [
        "",
        "All potentials and potential errors in the table are RMS over the "
        "same 512 solute-centred probes, in mEh/e. These are one-solute "
        "diagnostics, not mean benchmark errors across solutes.",
        "",
        "![Separation diagnostic](diagnostic.svg)",
        "",
        "The QM response tends to zero as the fragments separate. At 100 Å, "
        "both predictors retain a nonzero field; on these geometries GLIDER's "
        "residual exceeds the prior's. For the four-water system specifically, "
        "GLIDER's predicted RMS grows from 0.6965 at 20 Å to 1.5563 mEh/e "
        "at 100 Å as QM falls from 0.0065 to 0.00009 mEh/e. This directly "
        "confirms an out-of-range limitation in the frozen model. The "
        "short-range errors for this solute "
        "also vary sharply with distance, so this scan should not be used as a "
        "contact-range performance summary. The same method on other solutes "
        "would be needed to characterize the distribution of separation behavior.",
        "",
        "As an independent physical check, the counterpoise interaction energies "
        "at 100 Å are below 0.0001 kcal/mol in magnitude in both systems.",
        "",
    ]
    comparison = cpu_comparison(output, cpu_directories)
    if comparison:
        lines += [comparison, ""]
    lines += ["[Complete numeric table](summary.csv) · [Vector figure](diagnostic.svg) · "
              "[PDF figure](diagnostic.pdf) · [PNG figure](diagnostic.png)", ""]
    (output / "FINDINGS.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--cpu-reference", type=Path, nargs="*", default=[])
    args = parser.parse_args()
    rows = validate(args.directory)
    plot(rows, args.directory)
    note(rows, args.directory, args.cpu_reference)
    print(f"Validated {len(rows)} conditions and wrote figure and findings")


if __name__ == "__main__":
    main()
