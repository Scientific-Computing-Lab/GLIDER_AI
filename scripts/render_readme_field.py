#!/usr/bin/env python3
"""Compose the paper's exact Fig. 1c QM layers for the README.

The four transparent PNGs in assets/figure1_source are unchanged scientific
layers. This script only overlays the matching molecule and field layers and
adds labels, frames, and a qualitative colour key. It does not generate a new
electronic-structure result or rescale the two fields independently.
"""

from __future__ import annotations

from pathlib import Path

from matplotlib import font_manager
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "assets" / "figure1_source"
OUTPUT = ROOT / "assets" / "spatial-response.png"

INK = "#172C40"
MUTED = "#587083"
LINE = "#DEE7EA"
BLUE = "#245B97"
RED = "#C96558"
BACKGROUND = "#F7F9F8"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    family = "DejaVu Sans"
    path = font_manager.findfont(font_manager.FontProperties(family=family, weight="bold" if bold else "regular"))
    return ImageFont.truetype(path, size)


def field_panel(canvas: Image.Image, name: str, left: int, label: str, moment: str) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((left, 140, left + 715, 686), radius=24, fill="white", outline=LINE, width=2)
    draw.text((left + 35, 174), label, font=font(29, True), fill=INK)
    draw.text((left + 667, 178), moment, font=font(26), fill=MUTED, anchor="ra")

    field = Image.open(SOURCE / f"field_{name}.png").convert("RGBA")
    molecule = Image.open(SOURCE / f"molecule_{name}.png").convert("RGBA")
    assert field.size == molecule.size
    # Both scientific layers receive the same transform, preserving their alignment.
    scale = min(650 / field.width, 438 / field.height)
    size = (round(field.width * scale), round(field.height * scale))
    field = field.resize(size, Image.Resampling.LANCZOS)
    molecule = molecule.resize(size, Image.Resampling.LANCZOS)
    plot = Image.new("RGBA", size)
    plot.alpha_composite(field)
    plot.alpha_composite(molecule)
    x = left + (715 - size[0]) // 2
    y = 226 + (430 - size[1]) // 2
    canvas.alpha_composite(plot, (x, y))


def main() -> None:
    canvas = Image.new("RGBA", (1600, 804), BACKGROUND)
    draw = ImageDraw.Draw(canvas)
    draw.text((68, 48), "Similar magnitudes. Different spatial response.", font=font(51, True), fill=INK)
    draw.text((70, 111), "TWO WATER ARRANGEMENTS  /  EXACT QM RESPONSE ESP  /  PAPER FIG. 1c", font=font(21, True), fill=MUTED)
    field_panel(canvas, "equilibrium", 65, "Equilibrium", "|Δμ| = 0.464 D")
    field_panel(canvas, "cooperative", 820, "Cooperative", "|Δμ| = 0.468 D")

    # A qualitative key, matching the sign convention in the paper. The source
    # figure supplies the quantitative colour scale and full caption.
    draw.ellipse((71, 728, 89, 746), fill=BLUE)
    draw.text((100, 722), "negative response potential", font=font(22), fill=MUTED)
    draw.ellipse((441, 728, 459, 746), fill=RED)
    draw.text((470, 722), "positive response potential", font=font(22), fill=MUTED)
    draw.text((1530, 724), "QM reference, not model output", font=font(21), fill=MUTED, anchor="ra")
    canvas.convert("RGB").save(OUTPUT, optimize=True)


if __name__ == "__main__":
    main()
