#!/usr/bin/env python3
"""Compare the original contacts with the separately frozen follow-up."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
BLUE = "#0072B2"
TEAL = "#009E73"
AMBER = "#D55E00"
INK = "#202124"
MUTED = "#687078"
GRID = "#D8DDE2"
SPECIES = {"NH3": (BLUE, "NH$_3$"), "CH3OH": (TEAL, "CH$_3$OH"), "CH3CN": (AMBER, "CH$_3$CN")}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def jitter(identifier: str, amplitude: float = 0.18) -> float:
    digest = hashlib.sha256(identifier.encode()).digest()
    return (int.from_bytes(digest[:2], "big") / 65535 - 0.5) * 2 * amplitude


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, default=ROOT / "experiments/nonwater/contact_audit.csv")
    parser.add_argument("--followup", type=Path, default=ROOT / "experiments/nonwater_contact")
    args = parser.parse_args()
    original = read_csv(args.original)
    new = read_csv(args.followup / "results.csv")
    if len(original) != 72 or len(new) != 72:
        raise ValueError("Expected 72 original geometries and 36 two-method follow-up scores")
    by_case: dict[str, dict[str, dict[str, str]]] = {}
    for row in new:
        by_case.setdefault(row["config_id"], {})[row["method"]] = row
    if len(by_case) != 36 or any(set(item) != {"glider", "mace_polar_l"} for item in by_case.values()):
        raise ValueError("Follow-up methods or case coverage differ from the freeze")

    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.3, "axes.labelsize": 8.5,
        "axes.titlesize": 9, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "pdf.fonttype": 42,
        "svg.fonttype": "none", "svg.hashsalt": "glider_nonwater_contact_v1",
    })
    figure, (energy_axis, score_axis) = plt.subplots(
        1, 2, figsize=(7.1, 3.25), gridspec_kw={"width_ratios": [0.9, 1.2], "wspace": 0.4},
        facecolor="white",
    )
    for group, index, color in ((original, 0, MUTED), (by_case, 1, BLUE)):
        items = group if isinstance(group, list) else [methods["glider"] for methods in group.values()]
        for row in items:
            energy_axis.scatter(
                index + jitter(row["config_id"]),
                float(row["cp_interaction_energy_kcal_mol"]),
                s=16, color=color, alpha=0.7, edgecolor="white", linewidth=0.25,
                zorder=3,
            )
    energy_axis.axhline(0, color=INK, linewidth=0.7, zorder=2)
    energy_axis.set_yscale("symlog", linthresh=10, linscale=0.9)
    energy_axis.set_ylim(-15, 210)
    energy_axis.set_xlim(-0.5, 1.5)
    energy_axis.set_xticks([0, 1], ["Original\n72", "Follow-up\n36"])
    energy_axis.set_yticks([-10, -5, 0, 5, 10, 50, 100])
    energy_axis.set_yticklabels(["−10", "−5", "0", "5", "10", "50", "100"])
    energy_axis.set_ylabel("CP interaction energy (kcal mol$^{-1}$)")
    energy_axis.set_title("a  Fixed-geometry contact check", loc="left", weight="bold")
    energy_axis.grid(axis="y", color=GRID, linewidth=0.45)
    energy_axis.spines[["top", "right"]].set_visible(False)
    energy_axis.set_axisbelow(True)
    old_negative = sum(float(row["cp_interaction_energy_kcal_mol"]) < 0 for row in original)
    new_negative = sum(float(methods["glider"]["cp_interaction_energy_kcal_mol"]) < 0 for methods in by_case.values())
    energy_axis.text(0, -12.8, f"{old_negative}/72 < 0", ha="center", va="bottom", fontsize=7, color=MUTED)
    energy_axis.text(1, -12.8, f"{new_negative}/36 < 0", ha="center", va="bottom", fontsize=7, color=BLUE)

    maximum = 0.0
    shown_species: set[str] = set()
    for config_id, methods in sorted(by_case.items()):
        glider = methods["glider"]
        baseline = methods["mace_polar_l"]
        species = glider["perturbant"]
        color, label = SPECIES[species]
        x = float(baseline["esp_nrmse"])
        y = float(glider["esp_nrmse"])
        maximum = max(maximum, x, y)
        score_axis.scatter(x, y, s=24, color=color, alpha=0.8,
                           edgecolor="white", linewidth=0.35, zorder=3,
                           label=label if species not in shown_species else None)
        shown_species.add(species)
    limit = max(0.8, maximum * 1.08)
    score_axis.plot([0, limit], [0, limit], color=MUTED, linewidth=0.8, linestyle="--", zorder=1)
    score_axis.set_xlim(0, limit)
    score_axis.set_ylim(0, limit)
    score_axis.set_aspect("equal", adjustable="box")
    score_axis.set_xlabel("MACE-POLAR-1-L response-ESP NRMSE")
    score_axis.set_ylabel("GLIDER response-ESP NRMSE")
    score_axis.set_title("b  All 36 frozen follow-up scores", loc="left", weight="bold")
    score_axis.grid(color=GRID, linewidth=0.45)
    score_axis.spines[["top", "right"]].set_visible(False)
    score_axis.set_axisbelow(True)
    handles, labels = score_axis.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    score_axis.legend(unique.values(), unique.keys(), frameon=False, loc="upper left", fontsize=7)

    figure.subplots_adjust(left=0.1, right=0.98, top=0.87, bottom=0.2)
    for extension in ("svg", "pdf"):
        metadata = {"Creator": "GLIDER non-water contact follow-up"}
        metadata["Date" if extension == "svg" else "CreationDate"] = None
        figure.savefig(args.followup / f"contact_followup.{extension}",
                       facecolor="white", transparent=False,
                       metadata=metadata)
    plt.close(figure)


if __name__ == "__main__":
    main()
