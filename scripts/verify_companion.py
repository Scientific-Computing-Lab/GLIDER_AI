#!/usr/bin/env python3
"""Check that this curated paper companion is complete and internally coherent."""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IGNORED_DIRS = {".git", ".qa", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache",
                "build", "dist", "example_prediction"}


def ignored(path: Path) -> bool:
    parts = path.relative_to(ROOT).parts
    return any(part in IGNORED_DIRS or part.endswith(".egg-info") for part in parts) or \
        parts[:2] == ("third_party", "checkpoints")


def table(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def check_links() -> None:
    broken = []
    pattern = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    html_images = re.compile(r'<img\s+[^>]*src="([^"]+)"')
    for markdown in ROOT.rglob("*.md"):
        if ignored(markdown) or "provenance" in markdown.relative_to(ROOT).parts:
            continue
        content = markdown.read_text()
        for target in pattern.findall(content) + html_images.findall(content):
            target = target.split("#", 1)[0]
            if not target or ":" in target.split("/", 1)[0] or target.startswith("#"):
                continue
            if not (markdown.parent / target).exists():
                broken.append(f"{markdown.relative_to(ROOT)} -> {target}")
    assert not broken, "Broken local Markdown links:\n" + "\n".join(broken)


def check_headlines() -> None:
    rows = table("figures/figure_03/paired.csv")
    assert len(rows) == 3
    assert sum(int(r["n_molecules"]) for r in rows) == 56
    assert sum(int(r["glider_molecule_wins"]) for r in rows) == 54
    expected = [
        (0.321799354439422, 0.524476079542344, 38.6436546886518),
        (0.342174873318955, 0.571007635391217, 40.0752543204577),
        (0.426580605812233, 0.629706400723695, 32.2572225211651),
    ]
    for row, (glider, baseline, reduction) in zip(rows, expected):
        assert abs(float(row["glider_value"]) - glider) < 1e-12
        assert abs(float(row["comparator_value"]) - baseline) < 1e-12
        assert abs(float(row["relative_error_reduction_pct"]) - reduction) < 1e-10
        assert float(row["ci95_high"]) < 0
    nonwater = table("figures/figure_S15/nonwater_effects.csv")
    overall = next(r for r in nonwater if r["subset"] == "Overall")
    assert int(overall["n_solutes"]) == 12
    assert abs(float(overall["esp_glider"]) - 0.450371053714298) < 1e-12
    assert abs(float(overall["esp_comparator"]) - 0.617170202219818) < 1e-12
    ranks = table("figures/figure_04/energy_rank.csv")
    assert [int(r["rank"]) for r in ranks] == list(range(4, 13))
    w4 = ranks[0]
    assert abs(float(w4["glider"]) - 0.04868219341822043) < 1e-12
    assert abs(float(w4["mace_polar_l"]) - 0.11401221689887638) < 1e-12
    separation = table("figures/figure_S07/separation.csv")
    assert len(separation) == 40
    for system, glider_20, glider_100, prior_20, prior_100 in (
        ("single_water", 0.563891174, 1.540125107, 0.195369430, 0.113099723),
        ("whole_environment", 0.924588953, 2.842730048, 0.300614666, 0.239373936),
    ):
        def amplitude(method: str, distance: float) -> float:
            row = next(r for r in separation if r["system"] == system
                       and r["method"] == method and float(r["distance_A"]) == distance)
            return float(row["esp_rms_mEh_per_e"])
        assert abs(amplitude("glider", 20) - glider_20) < 1e-8
        assert abs(amplitude("glider", 100) - glider_100) < 1e-8
        assert abs(amplitude("averaged_prior", 20) - prior_20) < 1e-8
        assert abs(amplitude("averaged_prior", 100) - prior_100) < 1e-8
        assert glider_100 > glider_20 and prior_100 < prior_20
    components = table("experiments/dissociation_extended/posthoc_components.csv")
    assert len(components) == 280


def main():
    check_headlines()
    check_links()
    assert not (ROOT/".github/README.md").exists()
    for name in ['training','panel_1','panel_2','panel_3','nonwater','nonwater_contact','liquid','shell_size','heldout_water','heldout_water_pilot','distance_sweep','global_branch','dissociation','dissociation_extended']:
        assert (ROOT/'experiments'/name/'README.md').is_file(),name
    print('PASS: experiment guides, active local links, headline values and visual gallery')
if __name__=='__main__':main()
