#!/usr/bin/env python3
"""Create an identity-only FreeSolv table without loading property columns."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path

EXPECTED_SHA256 = "2d13f095713bc39b85f85dd7b4e5483fbb12fc694bf253bb1d92a4c4d484f260"
URL = "https://raw.githubusercontent.com/MobleyLab/FreeSolv/6c7d19b4b565537365ffd22006aa2cd4643200c6/database.txt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    if not args.raw_cache.exists():
        args.raw_cache.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(URL, args.raw_cache)
    if sha256(args.raw_cache) != EXPECTED_SHA256:
        raise RuntimeError("Pinned FreeSolv source checksum mismatch")

    records: list[tuple[str, str]] = []
    # Deliberately do not name, parse, store, or expose fields after index 1.
    with args.raw_cache.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split(";")
            records.append((fields[0].strip(), fields[1].strip()))
    if len(records) != 642 or len({key for key, _ in records}) != 642:
        raise RuntimeError("Unexpected pinned FreeSolv identity count")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["freesolv_id", "smiles"])
        writer.writerows(records)
    manifest = {
        "upstream_url": URL,
        "upstream_commit": "6c7d19b4b565537365ffd22006aa2cd4643200c6",
        "upstream_sha256": sha256(args.raw_cache),
        "extraction_script_sha256": sha256(Path(__file__)),
        "identity_only_sha256": sha256(args.output),
        "records": len(records),
        "fields_accessed": ["semicolon field 0: FreeSolv ID", "semicolon field 1: SMILES"],
        "fields_after_index_1_parsed_or_stored": False,
        "hydration_targets_accessed": False,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
