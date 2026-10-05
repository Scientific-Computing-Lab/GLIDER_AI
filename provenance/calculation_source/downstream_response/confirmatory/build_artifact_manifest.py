#!/usr/bin/env python3
"""Build the byte-size and SHA-256 inventory for the confirmatory package."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

COMMON = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(COMMON))
from common import sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[4]
RESULT_ROOT = ROOT / "results/posthoc_downstream_response_confirmatory"
OUTPUT = RESULT_ROOT / "artifact_manifest.csv"


def main() -> None:
    paths = sorted(
        path for path in RESULT_ROOT.rglob("*") if path.is_file() and path != OUTPUT
    )
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("relative_path", "bytes", "sha256"),
            lineterminator="\n",
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
