#!/usr/bin/env python3
"""Render the repository gallery from the released numerical tables.

The banner and architecture are explicitly conceptual. Result plots read the
unaltered CSVs in data/figure_data; no benchmark scores are recomputed here.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
DATA = ROOT / "data" / "figure_data"

NAVY = "#14263d"
INK = "#24364b"
MUTED = "#617388"
BLUE = "#1678b7"
ORANGE = "#df843d"
MAUVE = "#aa5a8e"
TEAL = "#36a398"
GREEN = "#6ca389"
PALE = "#f6f9fb"
LINE = "#dce6ec"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "svg.fonttype": "none",
        "svg.hashsalt": "glider-ai-paper-companion",
        "axes.edgecolor": LINE,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def rows(path: str) -> list[dict[str, str]]:
    with (DATA / path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def save(fig: plt.Figure, name: str) -> None:
    ASSETS.mkdir(exist_ok=True)
    fig.savefig(ASSETS / name, format="svg", metadata={"Date": None})
    plt.close(fig)


def box(ax, x, y, w, h, face, edge="none", radius=0.035, lw=1):
    shape = FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=lw,
    )
    ax.add_patch(shape)
    return shape


def banner() -> None:
    fig, ax = plt.subplots(figsize=(12, 4.35))
    fig.subplots_adjust(0, 0, 1, 1)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 4.35)
    ax.axis("off")
    box(ax, 0.02, 0.02, 11.96, 4.31, NAVY, radius=0.2)

    # Decorative, conceptual spatial field; it is not a measured potential.
    for x, y, w, h, c, a in [
        (9.08, 2.08, 4.5, 2.7, BLUE, 0.14),
        (10.40, 2.03, 3.4, 2.5, MAUVE, 0.13),
        (9.45, 1.54, 4.3, 2.1, TEAL, 0.13),
        (10.70, 1.25, 3.2, 1.9, ORANGE, 0.10),
    ]:
        ax.add_patch(Ellipse((x, y), w, h, facecolor=c, edgecolor="none", alpha=a))
    bonds = [(8.10, 2.10, 8.70, 2.62), (8.70, 2.62, 9.34, 2.11),
             (9.34, 2.11, 10.10, 2.36), (9.34, 2.11, 9.72, 1.43),
             (10.10, 2.36, 10.80, 2.83), (10.10, 2.36, 10.85, 1.84)]
    for x1, y1, x2, y2 in bonds:
        ax.plot([x1, x2], [y1, y2], color="#a3c5d4", alpha=0.55, lw=2.7, zorder=3)
    for x, y, radius, c in [
        (8.10, 2.10, .13, "#bed8df"), (8.70, 2.62, .20, "#f6f9fb"),
        (9.34, 2.11, .24, "#85b8d4"), (10.10, 2.36, .21, "#f6f9fb"),
        (9.72, 1.43, .13, "#d7e7ea"), (10.80, 2.83, .13, "#d7e7ea"),
        (10.85, 1.84, .13, "#d7e7ea"),
    ]:
        ax.add_patch(Circle((x, y), radius, facecolor=c, edgecolor=NAVY, lw=1.2, zorder=4))
    for radius in (1.12, 1.62, 2.14):
        ax.add_patch(Circle((9.59, 2.12), radius, fill=False,
                            edgecolor="#6caecb", linestyle=(0, (2, 7)),
                            lw=1, alpha=.32, zorder=2))

    ax.text(.60, 3.68, "GLIDER  /  PAPER COMPANION", color="#8ed3df", fontsize=11,
            fontweight="bold", va="center")
    ax.text(.60, 2.73, "Molecular response,", color="white", fontsize=31,
            fontweight="bold", va="center")
    ax.text(.60, 2.16, "kept spatial.", color="white", fontsize=31,
            fontweight="bold", va="center")
    ax.text(.63, 1.45, "From frozen polar pretraining to a reusable",
            color="#c3d3df", fontsize=13, va="center")
    ax.text(.63, 1.12, "interaction-induced electrostatic field.",
            color="#c3d3df", fontsize=13, va="center")
    box(ax, .61, .40, 1.67, .39, "#234a69", radius=.16)
    box(ax, 2.41, .40, 1.95, .39, "#234a69", radius=.16)
    box(ax, 4.49, .40, 2.09, .39, "#234a69", radius=.16)
    ax.text(1.45, .595, "48 labelled environments", color="white", fontsize=8.3,
            ha="center", va="center")
    ax.text(3.39, .595, "56 response-test solutes", color="white", fontsize=8.3,
            ha="center", va="center")
    ax.text(5.54, .595, "3 new neighbour species", color="white", fontsize=8.3,
            ha="center", va="center")
    save(fig, "glider-banner.svg")


def architecture() -> None:
    fig, ax = plt.subplots(figsize=(12, 4.45))
    fig.subplots_adjust(.02, .03, .98, .97)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 4.3)
    ax.axis("off")
    box(ax, .05, .06, 11.9, 4.18, PALE, radius=.18)
    ax.text(.48, 3.77, "How the field is built", fontsize=18, fontweight="bold", color=NAVY)
    ax.text(.48, 3.43, "Frozen inputs vary with geometry; only the response correction is learned from 48 labelled environments.",
            fontsize=9.9, color=MUTED)

    for x, y, w, h, fc, eyebrow, main, sub in [
        (.48, 1.97, 2.42, 1.13, "#e1f3f0", "FROZEN FEATURES", "MACE-POLAR-1-M", "geometry representation"),
        (.48, .65, 2.42, 1.13, "#e7eef8", "FROZEN PREDICTIONS", "MACE-M + MACE-L", "averaged starting response"),
        (4.07, 1.18, 2.80, 1.47, "#e5f1f8", "SPARSE RESPONSE LABELS", "GLIDER corrections", "local sites + global dipole"),
        (8.01, 1.18, 3.35, 1.47, "#f3eaf2", "CONSTRAINED OUTPUT", "Spatial response ESP", "net charge = 0  ·  dipole reconciled"),
    ]:
        box(ax, x, y, w, h, fc, edge=LINE, radius=.14)
        ax.text(x+.19, y+h-.24, eyebrow, color=MUTED, fontsize=7.9, fontweight="bold", va="center")
        ax.text(x+.19, y+h-.60, main, color=NAVY, fontsize=12.5, fontweight="bold", va="center")
        ax.text(x+.19, y+.20, sub, color=MUTED, fontsize=8.9, va="center")

    for a, b in [((2.90, 2.48), (4.07, 2.12)), ((2.90, 1.20), (4.07, 1.75)),
                 ((6.87, 1.92), (8.01, 1.92))]:
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=15,
                                     color="#557c96", linewidth=1.7))
    ax.text(5.74, .30, "48 three-water environments  ·  14 chemistries  ·  five frozen ensemble members",
            ha="center", fontsize=9.3, color=MUTED)
    save(fig, "architecture.svg")


def transfer() -> None:
    glider = rows("Fig4_glider.csv")
    baseline = rows("Fig4_mace_polar_l.csv")
    paired = rows("Fig4_paired.csv")
    fig, ax = plt.subplots(figsize=(11.6, 4.25))
    fig.subplots_adjust(left=.19, right=.95, top=.79, bottom=.22)
    fig.suptitle("Spatial response transfers to unseen-to-response solutes", x=.065,
                 y=.96, ha="left", fontsize=17, fontweight="bold", color=NAVY)
    fig.text(.065, .865, "Equal-solute response-ESP NRMSE; lower is better. Error bars show method-wise 95% intervals.",
             fontsize=10, color=MUTED)
    panel_names=[f"Panel {s}  ·  {n} solutes" for s,n in zip(('I','II','III'),(12,24,20))]
    y=[2,1,0]
    for i,(g,b,p) in enumerate(zip(glider,baseline,paired)):
        yy=y[i]
        for row,offset,color,label in [(b,.14,ORANGE,"Polar baseline"),(g,-.14,BLUE,"GLIDER")]:
            center=float(row["esp_nrmse"])
            low=float(row["esp_ci95_low"]);high=float(row["esp_ci95_high"])
            ax.plot([low,high],[yy+offset]*2,color=color,lw=2.3,alpha=.58,zorder=2)
            ax.plot([low,low],[yy+offset-.043,yy+offset+.043],color=color,lw=1.5,zorder=2)
            ax.plot([high,high],[yy+offset-.043,yy+offset+.043],color=color,lw=1.5,zorder=2)
            ax.scatter([center],[yy+offset],s=92,color=color,edgecolor="white",lw=1.6,zorder=4,
                       label=label if i==0 else None)
        ax.text(.84,yy,f'{float(p["relative_error_reduction_pct"]):.1f}% lower  ·  '
                f'{p["glider_molecule_wins"]}/{p["n_molecules"]} solutes',
                fontsize=9.2,color=INK,va="center")
    ax.set_yticks(y,panel_names)
    ax.set_xlim(.16,1.12)
    ax.set_ylim(-.48,2.48)
    ax.set_xticks([.2,.4,.6,.8,1.0])
    ax.set_xlabel("Response-ESP NRMSE",labelpad=10)
    ax.grid(axis="x",color=LINE,lw=.9)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.tick_params(axis="y",length=0,pad=11)
    ax.tick_params(axis="x",length=0,pad=7)
    ax.legend(frameon=False,ncol=2,loc="upper right",bbox_to_anchor=(1.0,1.27),fontsize=9.5)
    fig.text(.065,.052,"All predictions preceded the QM references. Foundation-pretraining exposure is unknown.",
             fontsize=8.8,color=MUTED)
    save(fig,"prospective-transfer.svg")


def coupling() -> None:
    rank=rows("downstream/energy_rank.csv")
    fig,(left,right)=plt.subplots(1,2,figsize=(11.6,4.2),gridspec_kw={"width_ratios":[.96,1.25]})
    fig.subplots_adjust(left=.17,right=.97,top=.68,bottom=.23,wspace=.32)
    fig.suptitle("Reuse one frozen response for a held-out water",x=.064,y=.96,
                 ha="left",fontsize=17,fontweight="bold",color=NAVY)
    fig.text(.064,.865,"Ten additional solutes. One-way electrostatic coupling, not total interaction energy.",
             fontsize=10,color=MUTED)
    w4=[("GLIDER",.0487,BLUE),("Zero response",.0853,"#798b99"),
        ("Polar baseline",.114,ORANGE),("Exact response dipole",.149,MAUVE)]
    ys=[3,2,1,0]
    for (label,value,color),y in zip(w4,ys):
        left.plot([0,value],[y,y],color=color,lw=3,alpha=.55)
        left.scatter(value,y,s=105,color=color,zorder=3)
        left.text(value+.007,y,f"{value:.3f}",va="center",fontsize=9.3,color=INK)
    left.set_yticks(ys,[x[0] for x in w4])
    left.set_xlim(0,.19);left.set_ylim(-.5,3.55)
    left.set_xticks([0,.05,.10,.15])
    left.set_xlabel("W4 energy MAE (kcal mol$^{-1}$)",labelpad=10)
    left.set_title("Single held-out probe",loc="left",fontsize=11.2,fontweight="bold",color=NAVY,pad=14)
    for name,color in [("glider",BLUE),("mace_polar_l",ORANGE),("exact_dipole",MAUVE)]:
        label={"glider":"GLIDER","mace_polar_l":"Polar baseline","exact_dipole":"Exact response dipole"}[name]
        right.plot([int(x["rank"]) for x in rank],[float(x[name]) for x in rank],
                   color=color,marker="o",markersize=5.4,lw=2,label=label)
    right.set_xticks(range(4,13))
    right.set_ylim(0,.24)
    right.set_xlabel("External-water rank (W4-W12)",labelpad=10)
    right.set_ylabel("Energy MAE (kcal mol$^{-1}$)")
    right.set_title("Same base, independent outer waters",loc="left",fontsize=11.2,
                    fontweight="bold",color=NAVY,pad=14)
    right.legend(frameon=False,fontsize=8.8,loc="lower left",ncol=3,
                 bbox_to_anchor=(0,1.32),borderaxespad=0,
                 handlelength=1.6,columnspacing=.9)
    for ax in (left,right):
        ax.grid(axis="x" if ax is left else "y",color=LINE,lw=.9)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():spine.set_visible(False)
        ax.tick_params(length=0,pad=7)
    fig.text(.064,.052,"The exact dipole uses the best single panel-wide origin for each endpoint among three predefined origins.",
             fontsize=8.8,color=MUTED)
    save(fig,"frozen-coupling.svg")


def distance() -> None:
    values=rows("FigS15_dense_energy.csv")
    fig,ax=plt.subplots(figsize=(11.6,4.2))
    fig.subplots_adjust(left=.11,right=.95,top=.74,bottom=.25)
    fig.suptitle("A local advantage, not a universal multipole verdict",x=.065,y=.96,
                 ha="left",fontsize=17,fontweight="bold",color=NAVY)
    fig.text(.065,.86,"Fixed-charge probe, ten frozen bases. Higher moments recover in the far field.",
             fontsize=10,color=MUTED)
    x=[float(r["clearance_A"]) for r in values]
    for col,label,color,marker in [
        ("glider","GLIDER",BLUE,"o"),
        ("exact_dipole","Exact dipole",MAUVE,"^"),
        ("exact_dipole_quadrupole","Exact D+Q",ORANGE,"D"),
        ("exact_dipole_quadrupole_octupole","Exact D+Q+O",GREEN,"s"),
    ]:
        ax.plot(x,[float(r[col]) for r in values],color=color,lw=2,marker=marker,
                markersize=5,label=label)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xlim(.6,15.3)
    ax.set_xlabel("Oxygen clearance from base van der Waals surface (Å)",labelpad=10)
    ax.set_ylabel("Coupling MAE (kcal mol$^{-1}$)")
    ax.grid(axis="y",which="major",color=LINE,lw=.9)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.tick_params(length=0,pad=7)
    ax.legend(frameon=False,ncol=4,loc="upper right",fontsize=8.8)
    save(fig,"distance-and-moments.svg")


def main() -> None:
    banner()
    architecture()
    transfer()
    coupling()
    distance()
    print("Wrote five SVG figures to assets/ from unchanged release data.")


if __name__ == "__main__":
    main()
