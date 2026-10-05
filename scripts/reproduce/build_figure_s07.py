#!/usr/bin/env python3
"""Join the archived separation summaries into the Figure S7 plot table."""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    ROOT / "experiments/dissociation/summary.csv",
    ROOT / "experiments/dissociation_extended/summary.csv",
)
OUTPUT = ROOT / "figures/figure_S07/separation.csv"


def main() -> None:
    rows = {}
    fieldnames = None
    for source in SOURCES:
        with source.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"Incompatible columns: {source}")
            for row in reader:
                key = (row["system"], row["method"], float(row["distance_A"]))
                if key in rows:
                    raise ValueError(f"Duplicate condition: {key}")
                rows[key] = row

    if len(rows) != 40:
        raise ValueError(f"Expected 40 method-by-distance rows, found {len(rows)}")
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for key in sorted(rows):
            writer.writerow(rows[key])
    print(f"Wrote {len(rows)} rows to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
