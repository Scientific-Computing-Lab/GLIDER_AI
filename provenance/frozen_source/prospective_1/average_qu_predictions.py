#!/usr/bin/env python3
"""Create the preregistered equal-weight MACE-POLAR M/L response control."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from response_learning import sha256


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    first, second = np.load(args.first), np.load(args.second)
    if not np.array_equal(
        first["config_ids"], second["config_ids"]
    ) or not np.array_equal(first["offsets"], second["offsets"]):
        raise RuntimeError("Prediction archives are not aligned")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        config_ids=first["config_ids"],
        offsets=first["offsets"],
        charges_e=0.5 * (first["charges_e"] + second["charges_e"]),
        dipoles_e_bohr=0.5 * (first["dipoles_e_bohr"] + second["dipoles_e_bohr"]),
    )
    manifest = {
        "method": "unfitted equal-weight MACE-POLAR-1-M/L q/u average",
        "first_sha256": sha256(args.first),
        "second_sha256": sha256(args.second),
        "output_sha256": sha256(args.output),
        "reference_labels_accessed": False,
    }
    args.output.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
