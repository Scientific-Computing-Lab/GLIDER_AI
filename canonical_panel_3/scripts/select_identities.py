#!/usr/bin/env python3
"""Frozen, label-blind deterministic FreeSolv identity selector."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors

from geometry import build_molecule_environments

ROOT = Path(__file__).resolve().parents[2]
ALLOWED = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
SEED = "GLIDER_CANONICAL_PANEL_3_FREESOLV_V052_20260814"
DESCRIPTORS = ["molecular_weight", "heavy_atoms", "hbd", "hba", "tpsa", "logp", "rings", "rotatable_bonds", "aromatic_fraction", "heteroatom_count"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(source_id: str, smiles: str, enforce_panel_eligibility: bool = True) -> dict | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    bonds = {Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC}
    heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if not (len(Chem.GetMolFrags(mol)) == 1 and sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 0):
        return None
    if enforce_panel_eligibility and not 4 <= heavy <= 15:
        return None
    if not {a.GetSymbol() for a in mol.GetAtoms()} <= ALLOWED:
        return None
    if any(a.GetNumRadicalElectrons() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    if any(b.GetBondType() not in bonds for b in mol.GetBonds()):
        return None
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    connectivity_key = Chem.MolToInchiKey(mol).split("-")[0]
    aromatic = sum(a.GetIsAromatic() for a in mol.GetAtoms())
    counts = Counter(a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms())
    symbols = {a.GetSymbol() for a in mol.GetAtoms()}
    if symbols & {"S", "P"}:
        stratum = 4
    elif symbols & {"F", "Cl", "Br", "I"}:
        stratum = 5
    elif "N" in symbols or any(a.GetIsAromatic() and a.GetSymbol() in {"N", "O", "S"} for a in mol.GetAtoms()):
        stratum = 3
    elif "O" in symbols:
        stratum = 2
    else:
        stratum = 1
    return {
        "freesolv_id": source_id,
        "input_smiles": smiles,
        "canonical_smiles": canonical,
        "connectivity_inchikey": connectivity_key,
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "elemental_composition": ";".join(f"{k}:{counts[k]}" for k in sorted(counts)),
        "molecular_weight": Descriptors.MolWt(mol),
        "heavy_atoms": heavy,
        "formal_charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "logp": Crippen.MolLogP(mol),
        "rings": Lipinski.RingCount(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "aromatic_fraction": aromatic / heavy,
        "heteroatom_count": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "stratum": stratum,
        "_mol": mol,
        "_fp": FPGEN.GetFingerprint(mol),
    }


def historical() -> list[dict]:
    result = []
    specifications = [
        (ROOT / "data/development/configuration_registry.csv", "ORIGINAL_DEVELOPMENT"),
        (ROOT / "data/prospective_1/configuration_registry.csv", None),
        (ROOT / "data/prospective_2/configuration_registry.csv", None),
    ]
    for path, pool_status in specifications:
        table = pd.read_csv(path)
        if pool_status is not None:
            table = table[table.pool_status == pool_status]
        for row in table.drop_duplicates("molecule_id").itertuples():
            parsed = record(str(row.molecule_id), str(row.smiles), enforce_panel_eligibility=False)
            if parsed is None:
                raise RuntimeError(f"Historical identity did not parse: {row.molecule_id}")
            result.append(parsed)
    if len(result) != 50 or len({row["canonical_smiles"] for row in result}) != 50:
        raise RuntimeError("Historical record is not exactly 50 unique identities")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--identities", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    source = pd.read_csv(args.identities)
    if list(source.columns) != ["freesolv_id", "smiles"]:
        raise RuntimeError("Selector accepts the identity-only schema exclusively")
    history = historical()
    historical_smiles = {row["canonical_smiles"] for row in history}
    historical_keys = {row["connectivity_inchikey"] for row in history}
    candidates, exclusions = [], []
    for row in source.itertuples(index=False):
        parsed = record(str(row.freesolv_id), str(row.smiles))
        if parsed is None:
            exclusions.append({"freesolv_id": row.freesolv_id, "reason": "frozen_eligibility"})
            continue
        similarities = DataStructs.BulkTanimotoSimilarity(parsed["_fp"], [x["_fp"] for x in history])
        parsed["maximum_historical_morgan_similarity"] = max(similarities)
        parsed["nearest_historical_id"] = history[int(np.argmax(similarities))]["freesolv_id"]
        if parsed["canonical_smiles"] in historical_smiles or parsed["connectivity_inchikey"] in historical_keys:
            exclusions.append({"freesolv_id": row.freesolv_id, "reason": "exact_historical_overlap"})
        elif parsed["maximum_historical_morgan_similarity"] >= 0.70:
            exclusions.append({"freesolv_id": row.freesolv_id, "reason": "historical_similarity_ge_0.70"})
        else:
            candidates.append(parsed)
    if not candidates:
        raise RuntimeError("No eligible candidates")
    matrix = np.asarray([[float(row[key]) for key in DESCRIPTORS] for row in candidates])
    median = np.median(matrix, axis=0)
    q25, q75 = np.percentile(matrix, [25, 75], axis=0)
    scale = np.where(q75 > q25, q75 - q25, 1.0)
    for row, vector in zip(candidates, (matrix - median) / scale):
        row["_z"] = vector
    hist_matrix = np.asarray([[float(row[key]) for key in DESCRIPTORS] for row in history])
    for row, vector in zip(history, (hist_matrix - median) / scale):
        row["_z"] = vector

    def distance(left: dict, right: dict) -> float:
        fingerprint = 1.0 - DataStructs.TanimotoSimilarity(left["_fp"], right["_fp"])
        descriptor = min(float(np.linalg.norm(left["_z"] - right["_z"])) / np.sqrt(10.0), 1.0)
        return 0.70 * fingerprint + 0.30 * descriptor

    chosen: list[dict] = []
    failures: list[dict] = []
    protocol_hash = sha256(args.protocol)
    requests = ["aliphatic", "aromatic", "aliphatic", "aromatic"]
    for round_index in range(4):
        for stratum in range(1, 6):
            chosen_ids = {row["freesolv_id"] for row in chosen}
            pool = [row for row in candidates if row["stratum"] == stratum and row["freesolv_id"] not in chosen_ids and all(f["freesolv_id"] != row["freesolv_id"] for f in failures)]
            if stratum == 1:
                aromatic = requests[round_index] == "aromatic"
                pool = [row for row in pool if (row["aromatic_fraction"] > 0) == aromatic]
            ranking = []
            reference = history + chosen
            for row in pool:
                score = min(distance(row, other) for other in reference)
                tie = hashlib.sha256(f"{row['freesolv_id']}|{SEED}".encode()).hexdigest()
                ranking.append((-score, tie, row, score))
            ranking.sort(key=lambda item: (item[0], item[1]))
            accepted = None
            for _, _, row, score in ranking:
                try:
                    frames = build_molecule_environments(row["freesolv_id"], row["canonical_smiles"], f"stratum_{stratum}", protocol_hash)
                    if len(frames) != 4:
                        raise RuntimeError("not four environments")
                except Exception as error:
                    failures.append({"freesolv_id": row["freesolv_id"], "stratum": stratum, "reason": type(error).__name__, "detail": str(error)})
                    continue
                row["selection_round"] = round_index + 1
                row["maxmin_score_at_selection"] = score
                row["tie_break_sha256"] = hashlib.sha256(f"{row['freesolv_id']}|{SEED}".encode()).hexdigest()
                accepted = row
                chosen.append(row)
                break
            if accepted is None:
                raise RuntimeError(f"No constructable candidate for stratum {stratum}, round {round_index + 1}")
    if len(chosen) != 20 or Counter(row["stratum"] for row in chosen) != Counter({1: 4, 2: 4, 3: 4, 4: 4, 5: 4}):
        raise RuntimeError("Frozen 4x5 selection failed")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    private = {"_mol", "_fp", "_z"}
    pd.DataFrame([{key: value for key, value in row.items() if key not in private} for row in candidates]).sort_values(["stratum", "freesolv_id"]).to_csv(args.output_dir / "eligible_pool.csv", index=False)
    pd.DataFrame([{key: value for key, value in row.items() if key not in private} for row in chosen]).sort_values(["stratum", "selection_round"]).to_csv(args.output_dir / "selected_molecules.csv", index=False)
    pd.DataFrame(exclusions + failures).to_csv(args.output_dir / "exclusions_and_construction_failures.csv", index=False)
    manifest = {
        "identity_only_source_sha256": sha256(args.identities),
        "selector_sha256": sha256(Path(__file__)),
        "geometry_preflight_sha256": sha256(Path(__file__).with_name("geometry.py")),
        "protocol_sha256": protocol_hash,
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "eligible_count": len(candidates),
        "selected_count": len(chosen),
        "stratum_counts": {str(key): value for key, value in sorted(Counter(row["stratum"] for row in chosen).items())},
        "construction_failures": failures,
        "hydration_targets_accessed": False,
        "model_predictions_accessed": False,
        "response_labels_accessed": False,
    }
    manifest["eligible_pool_sha256"] = sha256(args.output_dir / "eligible_pool.csv")
    manifest["selected_molecules_sha256"] = sha256(args.output_dir / "selected_molecules.csv")
    (args.output_dir / "SELECTION_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
