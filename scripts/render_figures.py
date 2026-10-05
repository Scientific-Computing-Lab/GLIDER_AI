#!/usr/bin/env python3
"""Render the repository's vector reading aids from released numerical tables.

The architecture is conceptual. Result plots read the unaltered CSVs in
figures/; no benchmark scores are recomputed here.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
DATA = ROOT / "figures"

NAVY = "#14263d"
INK = "#24364b"
MUTED = "#617388"
BLUE = "#1678b7"
ORANGE = "#df843d"
MAUVE = "#aa5a8e"
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
    path = ASSETS / name
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


def box(ax, x, y, w, h, face, edge="none", radius=0.035, lw=1):
    shape = FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=lw,
    )
    ax.add_patch(shape)
    return shape


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
    ax.text(5.74, .30, "48 three-/four-water environments  ·  14 chemistries  ·  five frozen ensemble members",
            ha="center", fontsize=9.3, color=MUTED)
    save(fig, "architecture.svg")


def transfer() -> None:
    glider = rows("figure_03/glider.csv")
    baseline = rows("figure_03/mace_polar_l.csv")
    paired = rows("figure_03/paired.csv")
    fig, ax = plt.subplots(figsize=(11.6, 4.25))
    fig.subplots_adjust(left=.19, right=.95, top=.79, bottom=.22)
    fig.suptitle("Spatial response transfers beyond response supervision", x=.065,
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
    fig.text(.065,.052,"GLIDER predictions preceded response QM. Some public baselines were evaluated later.",
             fontsize=8.8,color=MUTED)
    save(fig,"prospective-transfer.svg")


def coupling() -> None:
    rank=rows("figure_04/energy_rank.csv")
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
    values=rows("figure_05/dense_energy.csv")
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
    architecture()
    transfer()
    coupling()
    distance()
    print("Wrote five SVG figures to assets/ from unchanged release data.")


def separation():
    import pandas as pd
    from matplotlib.ticker import ScalarFormatter
    data=pd.read_csv(DATA/'figure_S07/separation.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,3.8))
    for ax,system,title,upper in zip(
        axes,
        ['single_water','whole_environment'],
        ['Solute + one water','Solute + intact water cluster'],
        [3,200],
    ):
        for method,color,marker,label in [
            ('glider',BLUE,'o','GLIDER'),
            ('averaged_prior',ORANGE,'s','Frozen averaged prior'),
        ]:
            q=data[(data.system==system)&(data.method==method)].sort_values('distance_A')
            ax.plot(q.distance_A,q.esp_rms_mEh_per_e,color=color,lw=1.7,marker=marker,ms=4,label=label)
        ax.set(xscale='log',yscale='log',xlim=(2.8,110),ylim=(.09,upper),
               xlabel='Minimum fragment distance (Å)',
               ylabel='Predicted response-potential RMS (mEh/e)',title=title)
        ax.set_xticks([3,5,10,20,50,100])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.spines[['top','right']].set_visible(False)
        ax.grid(alpha=.2)
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('The predicted response grows beyond 20 Å',fontsize=13,fontweight='bold',y=.98)
    fig.tight_layout(rect=(0,0,1,.91));save(fig,'separation-control.svg')

if __name__ == "__main__":
    separation()
    main()
