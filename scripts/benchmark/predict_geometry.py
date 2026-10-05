#!/usr/bin/env python3
"""Geometry-to-response GLIDER inference using official MACE-POLAR checkpoints."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FROZEN = ROOT / "provenance/frozen_source/prospective_1"


def run(*arguments: str) -> None:
    subprocess.run([sys.executable, *arguments], check=True, cwd=ROOT)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, action="append", required=True)
    parser.add_argument("--mace-polar-m", type=Path, required=True)
    parser.add_argument("--mace-polar-l", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    configuration_args = [
        value
        for path in args.configurations
        for value in ("--configurations", str(path.resolve()))
    ]
    with tempfile.TemporaryDirectory(prefix="glider-inference-") as temporary:
        work = Path(temporary)
        features = work / "features_m.npz"
        prior_m = work / "prior_m.npz"
        prior_l = work / "prior_l.npz"
        prior_ml = work / "prior_ml.npz"
        run(
            str(FROZEN / "extract_mace_features.py"),
            *configuration_args,
            "--checkpoint",
            str(args.mace_polar_m.resolve()),
            "--output",
            str(features),
            "--device",
            args.device,
        )
        for checkpoint, output in (
            (args.mace_polar_m, prior_m),
            (args.mace_polar_l, prior_l),
        ):
            run(
                str(FROZEN / "predict_mace_polar_cp_response.py"),
                *configuration_args,
                "--checkpoint",
                str(checkpoint.resolve()),
                "--output",
                str(output),
                "--device",
                args.device,
            )
        run(
            str(FROZEN / "average_qu_predictions.py"),
            "--first",
            str(prior_m),
            "--second",
            str(prior_l),
            "--output",
            str(prior_ml),
        )
        run(
            str(ROOT / "scripts/benchmark/predict_from_features.py"),
            *configuration_args,
            "--features",
            str(features),
            "--checkpoint",
            str(ROOT / "checkpoints/glider_site_response_ensemble.pt"),
            "--base-predictions",
            str(prior_ml),
            "--output",
            str(args.output.resolve()),
        )


if __name__ == "__main__":
    main()
