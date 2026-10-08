#!/usr/bin/env python3
"""Draw the post hoc contact audit of the original non-water panel.

Input is produced by scripts/reproduce/audit_nonwater_contacts.py.  The plot
never filters the frozen panel or replaces its originally reported score.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BLUE = "#0072B2"
INK = "#202124"
MUTED = "#687078"
GRID = "#D8DDE2"
RED = "#B33A3A"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "experiments/nonwater/contact_audit.csv"
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "experiments/nonwater")
    args = parser.parse_args()
    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 72:
        raise ValueError("Expected all 72 frozen non-water configurations")

    groups = [
        (r"$E_{\mathrm{int}}<0$", lambda x: x < 0),
        (r"$0\leq E_{\mathrm{int}}\leq 10$", lambda x: 0 <= x <= 10),
        (r"$E_{\mathrm{int}}>10$", lambda x: x > 10),
    ]
    grouped = [
        [row for row in rows if criterion(float(row["cp_interaction_energy_kcal_mol"]))]
        for _, criterion in groups
    ]
    if [len(group) for group in grouped] != [10, 10, 52]:
        raise ValueError("Energy-bin counts differ from the validated contact audit")

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.2,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 7.6,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "glider_nonwater_audit_v1",
        }
    )
    fig = plt.figure(figsize=(7.1, 3.0), facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.55], wspace=0.38)
    ax_count = fig.add_subplot(grid[0, 0])
    ax_effect = fig.add_subplot(grid[0, 1])

    count = [len(group) for group in grouped]
    colors = ["#B8CBD9", "#739EBB", BLUE]
    ax_count.barh(range(3), count, color=colors, height=0.54, edgecolor="none")
    ax_count.set_yticks(range(3), [r"$<0$", r"$0$ to $10$", r"$>10$"])
    ax_count.invert_yaxis()
    ax_count.set_xlim(0, 60)
    ax_count.set_xlabel("Configurations")
    ax_count.set_title("a  CP interaction energy (kcal mol$^{-1}$)", loc="left", weight="bold")
    for index, value in enumerate(count):
        ax_count.text(value + 1, index, str(value), va="center", fontsize=8, color=INK)
    ax_count.spines[["top", "right", "left"]].set_visible(False)
    ax_count.tick_params(axis="y", length=0)
    ax_count.grid(axis="x", color=GRID, linewidth=0.45)
    ax_count.set_axisbelow(True)

    for index, group in enumerate(grouped):
        for row in group:
            tag = hashlib.sha256(str(row["config_id"]).encode()).digest()
            jitter = (int.from_bytes(tag[:2], "big") / 65535 - 0.5) * 0.48
            ax_effect.scatter(
                index + jitter,
                float(row["glider_minus_mace_l_nrmse"]),
                s=18,
                facecolors=BLUE,
                edgecolors="white",
                linewidths=0.35,
                alpha=0.78,
                zorder=3,
            )
    ax_effect.axhline(0, color=RED, linewidth=0.8, zorder=2)
    ax_effect.set_xlim(-0.55, 2.55)
    ax_effect.set_ylim(-0.78, 0.18)
    ax_effect.set_xticks(
        range(3),
        ["$<0$\n2 solutes", "$0$ to $10$\n5 solutes", "$>10$\n11 solutes"],
    )
    ax_effect.set_title("b  GLIDER − MACE-L response-ESP NRMSE", loc="left", weight="bold")
    ax_effect.grid(axis="y", color=GRID, linewidth=0.45)
    ax_effect.set_axisbelow(True)
    ax_effect.spines[["top", "right"]].set_visible(False)
    ax_effect.text(
        0.99,
        0.96,
        "Below zero favours GLIDER",
        transform=ax_effect.transAxes,
        ha="right",
        va="top",
        color=MUTED,
        fontsize=7,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.87, bottom=0.22, wspace=0.31)
    for suffix in ("svg", "pdf"):
        metadata = {"Creator": "GLIDER non-water contact audit"}
        metadata["Date" if suffix == "svg" else "CreationDate"] = None
        fig.savefig(
            args.output_dir / f"contact_audit.{suffix}",
            facecolor="white",
            transparent=False,
            metadata=metadata,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
