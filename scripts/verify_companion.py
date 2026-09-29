#!/usr/bin/env python3
"""Check that this curated paper companion is complete and internally coherent."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "PROVENANCE.json"
ASSET_NAMES = {
    "glider-banner.svg",
    "architecture.svg",
    "prospective-transfer.svg",
    "frozen-coupling.svg",
    "distance-and-moments.svg",
}
IGNORED_DIRS = {".git", ".qa", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache",
                "build", "dist", "example_prediction"}


def ignored(path: Path) -> bool:
    parts = path.relative_to(ROOT).parts
    return any(part in IGNORED_DIRS or part.endswith(".egg-info") for part in parts) or \
        parts[:2] == ("third_party", "checkpoints")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_digest(directory: Path) -> str:
    """Hash immutable archive files, excluding newly written directory guides."""
    h = hashlib.sha256()
    for path in sorted(p for p in directory.rglob("*") if p.is_file() and p.name != "README.md"):
        h.update(str(path.relative_to(directory)).encode() + b"\0")
        h.update(digest(path).encode() + b"\n")
    return h.hexdigest()


def table(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def check_links() -> None:
    broken = []
    pattern = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    for markdown in ROOT.rglob("*.md"):
        if ignored(markdown):
            continue
        for target in pattern.findall(markdown.read_text()):
            target = target.split("#", 1)[0]
            if not target or ":" in target.split("/", 1)[0] or target.startswith("#"):
                continue
            if not (markdown.parent / target).exists():
                broken.append(f"{markdown.relative_to(ROOT)} -> {target}")
    assert not broken, "Broken local Markdown links:\n" + "\n".join(broken)


def check_headlines() -> None:
    rows = table("data/figure_data/Fig4_paired.csv")
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
    nonwater = table("data/figure_data/Fig5_nonwater_effects.csv")
    overall = next(r for r in nonwater if r["subset"] == "Overall")
    assert int(overall["n_solutes"]) == 12
    assert abs(float(overall["esp_glider"]) - 0.450371053714298) < 1e-12
    assert abs(float(overall["esp_comparator"]) - 0.617170202219818) < 1e-12
    ranks = table("data/figure_data/downstream/energy_rank.csv")
    assert [int(r["rank"]) for r in ranks] == list(range(4, 13))
    w4 = ranks[0]
    assert abs(float(w4["glider"]) - 0.04868219341822043) < 1e-12
    assert abs(float(w4["mace_polar_l"]) - 0.11401221689887638) < 1e-12


def main() -> None:
    provenance = json.loads(MANIFEST.read_text())
    count = 0
    for group in (
        "code_and_checkpoint_sha256",
        "camera_ready_scientific_data_sha256",
        "derived_plotting_data_sha256",
    ):
        for relative, expected in provenance[group].items():
            path = ROOT / relative
            assert path.is_file(), relative
            assert digest(path) == expected, f"Content hash changed: {relative}"
            count += 1
    for relative, expected in provenance["benchmark_archive_tree_sha256"].items():
        assert tree_digest(ROOT / relative) == expected, f"Archive changed: {relative}"
    check_headlines()
    check_links()
    directories = [
        p for p in ROOT.rglob("*")
        if p.is_dir() and not ignored(p)
    ]
    missing_guides = [str(p.relative_to(ROOT)) for p in directories if not (p / "README.md").is_file()]
    assert not missing_guides, f"Directories without README: {missing_guides}"
    forbidden = {".pdf", ".tex", ".zip", ".docx", ".pptx"}
    copied_manuscripts = [
        str(p.relative_to(ROOT)) for p in ROOT.rglob("*")
        if p.is_file() and p.suffix.lower() in forbidden and ".git" not in p.parts
    ]
    assert not copied_manuscripts, f"Manuscript/build files included: {copied_manuscripts}"
    assert {p.name for p in (ROOT / "assets").glob("*.svg")} == ASSET_NAMES
    assert (ROOT / "examples/one_response_geometry.extxyz").read_text().splitlines()[0] == "28"
    print(f"Companion verified: {count} exact source files, {len(directories)} directory guides, "
          "headline values, links and SVG gallery.")


if __name__ == "__main__":
    main()
