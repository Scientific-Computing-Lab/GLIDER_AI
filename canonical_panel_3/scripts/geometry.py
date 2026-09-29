"""Byte-stable Panel-III wrapper around the frozen panel-II geometry rules."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(WORKSPACE / "complete_interaction_hamiltonian/src"))
sys.path.insert(0, str(WORKSPACE / "manybody_completion/src"))
from build_configurations import build_cluster, conformers, make_atoms  # noqa: E402
from build_clusters import cooperative_chain, reorient  # noqa: E402


def build_molecule_environments(source_id: str, smiles: str, stratum: str, protocol_hash: str):
    seed = int(hashlib.sha256(f"{protocol_hash}|{source_id}|{smiles}".encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    molecule, coordinates = conformers(smiles, seed % (2**31 - 1))
    solute = coordinates[0]
    heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    hetero = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() in (7, 8, 9, 15, 16, 17, 35, 53)]
    anchors = hetero + [index for index in heavy if index not in hetero]
    equilibrium = build_cluster(solute, anchors, 3, rng, "equilibrium")
    specifications = [
        ("equilibrium", equilibrium),
        ("cooperative", cooperative_chain(solute, anchors, 3, rng)),
        ("orientation", reorient(equilibrium, rng)),
        ("compressed", cooperative_chain(solute, anchors, 3, rng, compressed=True)),
    ]
    result = []
    for regime, waters in specifications:
        config_id = f"canonical_p3__{source_id}__cluster_3water__{regime}"
        metadata = {
            "config_id": config_id,
            "molecule_id": source_id,
            "smiles": smiles,
            "split": "canonical_prospective_3",
            "chemical_role": stratum,
            "category": "cluster_3water",
            "motif": regime,
            "regime": regime,
            "conformer_index": 0,
        }
        result.append(make_atoms(molecule, solute, waters, metadata))
    return result
