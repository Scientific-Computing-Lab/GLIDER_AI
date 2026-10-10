#!/usr/bin/env python3
"""Render the journal Fig. 9 reading aid from the audited ARROW release data."""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/arrow_comparison"
METHODS = ("glider", "arrow_smeared", "mace_polar_l")
COLORS = {"glider": "#0072B2", "arrow_smeared": "#7B4AB2", "mace_polar_l": "#D55E00"}
LABELS = {"glider": "GLIDER", "arrow_smeared": "ARROW", "mace_polar_l": "MACE-L"}
MARKERS = {"glider": "o", "arrow_smeared": "^", "mace_polar_l": "s"}


def main() -> None:
    audit = json.loads((BASE / "audit.json").read_text())
    sets = ("panel_1", "panel_2", "panel_3", "liquid", "shell_size")
    set_labels = ("Panel I", "Panel II", "Panel III", "Liquid", "Shell")
    bias = defaultdict(list)
    with (BASE / "shell_size_signed_bias_per_config.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["method"] in METHODS:
                bias[(int(row["n_waters"]), row["method"])].append(
                    float(row["signed_error_solute_surface_mEh"])
                )

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "axes.edgecolor": "#687078",
                         "text.color": "#202124", "axes.labelcolor": "#202124"})
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.4, 4.0), layout="constrained")
    fig.patch.set_facecolor("#FFFFFF")
    for ax in (left, right):
        ax.set_facecolor("#FFFFFF")
        ax.grid(axis="y", color="#D8DDE2", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_axisbelow(True)

    for method, offset in zip(METHODS, (-0.16, 0.0, 0.16)):
        values = [audit["sets"][name]["nrmse"][method] for name in sets]
        left.scatter([i + offset for i in range(len(sets))], values,
                     color=COLORS[method], marker=MARKERS[method], s=48,
                     label=LABELS[method], zorder=3)
    left.set_xticks(range(len(sets)), set_labels)
    left.set_ylim(0, 0.87)
    left.set_ylabel("Response-ESP NRMSE")
    left.set_title("a   Accuracy on shared chemistry", loc="left", weight="bold", pad=15)
    left.legend(frameon=False, loc="upper left", ncol=3, fontsize=9)

    waters = (1, 3, 6, 12)
    for method in METHODS:
        means = [statistics.mean(bias[(size, method)]) for size in waters]
        if any(len(bias[(size, method)]) != 6 for size in waters):
            raise ValueError(f"Incomplete signed-bias series: {method}")
        right.plot(waters, means, color=COLORS[method], marker=MARKERS[method],
                   markersize=6, linewidth=1.8, label=LABELS[method])
    right.axhline(0, color="#687078", linewidth=0.9, linestyle=":")
    right.set_xticks(waters)
    right.set_ylim(-2.35, 2.35)
    right.set_xlabel("Number of waters")
    right.set_ylabel("Solute-surface signed error (mEh/e)")
    right.set_title("b   Signed error as water shells grow", loc="left", weight="bold", pad=15)
    right.legend(frameon=False, loc="lower left", ncol=3, fontsize=9)

    output = BASE / "figure_09.svg"
    fig.savefig(output, format="svg", facecolor="#FFFFFF",
                metadata={"Title": "Independent ARROW comparison on shared GLIDER configurations",
                          "Description": "Panel a compares response-ESP NRMSE in five shared sets. Panel b shows signed solute-surface error from one to 12 waters."})
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
