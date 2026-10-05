#!/usr/bin/env python3
"""Build a byte-size and SHA-256 inventory of the downstream result package."""

from __future__ import annotations

import csv
from pathlib import Path

from common import sha256

ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response"
OUTPUT = RESULT_ROOT / "artifact_manifest.csv"


def main() -> None:
    paths = sorted(
        path
        for path in RESULT_ROOT.rglob("*")
        if path.is_file() and path != OUTPUT
    )
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("relative_path", "bytes", "sha256")
        )
        writer.writeheader()
        for path in paths:
            writer.writerow(
                {
                    "relative_path": path.relative_to(RESULT_ROOT).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    print(f"wrote {len(paths)} entries to {OUTPUT}")


if __name__ == "__main__":
    main()
