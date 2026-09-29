#!/usr/bin/env python3
"""Download official MACE-POLAR-1 checkpoints and verify publication hashes."""

from __future__ import annotations

import argparse
import hashlib
import urllib.request
from pathlib import Path

BASE = "https://github.com/ACEsuit/mace-foundations/releases/download/mace_polar_1"
CHECKPOINTS = {
    "S": (
        "MACE-POLAR-1-S.model",
        "e4495612037b3b3312633182882a38a694ecac9ea0be2b9889ac0b2a84a99510",
    ),
    "M": (
        "MACE-POLAR-1-M.model",
        "fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a",
    ),
    "L": (
        "MACE-POLAR-1-L.model",
        "9f65f8dc6ddaff1d631e299cb531376a7da5e68d1bef04f34a2d5073d5ef114b",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sizes", nargs="+", choices=sorted(CHECKPOINTS), default=["M", "L"])
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for size in args.sizes:
        name, expected = CHECKPOINTS[size]
        destination = args.output_dir / name
        if not destination.exists():
            urllib.request.urlretrieve(f"{BASE}/{name}", destination)
        observed = sha256(destination)
        if observed != expected:
            raise RuntimeError(f"Hash mismatch for {name}: {observed}")
        print(f"{name}: verified {observed}")


if __name__ == "__main__":
    main()
