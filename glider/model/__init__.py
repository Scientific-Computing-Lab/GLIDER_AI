"""Frozen-model provenance helpers.

The byte-exact training/inference implementation is retained under
``evidence/frozen_source/prospective_1``. It remains the authoritative path for
historical predictions; this package deliberately does not relabel refactored
code as the source of the frozen result.
"""

from pathlib import Path


def authoritative_source(root: Path | None = None) -> Path:
    base = root or Path(__file__).resolve().parents[2]
    return base / "evidence/frozen_source/prospective_1/response_learning.py"


__all__ = ["authoritative_source"]
