#!/usr/bin/env python3
"""Build frozen chemistry-disjoint solute--water interaction configurations."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import write
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign


ROOT = Path(__file__).resolve().parents[1]
SEED = 20260812
WATER = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.9572 * np.cos(np.deg2rad(52.26)), 0.9572 * np.sin(np.deg2rad(52.26)), 0.0],
        [0.9572 * np.cos(np.deg2rad(52.26)), -0.9572 * np.sin(np.deg2rad(52.26)), 0.0],
    ]
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def rotation_from_x(target: np.ndarray) -> np.ndarray:
    target = np.asarray(target, float) / np.linalg.norm(target)
    source = np.array([1.0, 0.0, 0.0])
    cross = np.cross(source, target)
    sine = np.linalg.norm(cross)
    cosine = float(source @ target)
    if sine < 1e-12:
        return np.eye(3) if cosine > 0 else np.diag([-1.0, -1.0, 1.0])
    skew = np.array(
        [[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]]
    )
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / sine**2)


def conformers(smiles: str, seed: int) -> tuple[Chem.Mol, list[np.ndarray]]:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ids = list(
        AllChem.EmbedMultipleConfs(
            mol,
            numConfs=16,
            randomSeed=seed,
            pruneRmsThresh=0.12,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
        )
    )
    if not ids:
        raise RuntimeError(f"Embedding failed for {smiles}")
    results = AllChem.MMFFOptimizeMoleculeConfs(
        mol, mmffVariant="MMFF94s", maxIters=1000
    )
    order = sorted(range(len(ids)), key=lambda i: (float(results[i][1]), ids[i]))
    selected = [order[0]]
    for idx in order[1:]:
        if (
            rdMolAlign.GetBestRMS(mol, mol, prbId=ids[idx], refId=ids[selected[0]])
            >= 0.35
        ):
            selected.append(idx)
            break
    coords = []
    for idx in selected:
        conf = mol.GetConformer(ids[idx])
        coords.append(
            np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        )
    return mol, coords


def direction_outward(
    coords: np.ndarray, anchor: int, rng: np.random.Generator
) -> np.ndarray:
    direction = coords[anchor] - coords.mean(axis=0) + 0.35 * rng.normal(size=3)
    if np.linalg.norm(direction) < 0.2:
        direction = rng.normal(size=3)
    return direction / np.linalg.norm(direction)


def orient_water(radial: np.ndarray, mode: str, rng: np.random.Generator) -> np.ndarray:
    radial = radial / np.linalg.norm(radial)
    if mode == "water_acceptor":
        # Water oxygen approaches a solute donor: hydrogens point away.
        dipole = radial
        rotation = rotation_from_x(dipole)
    elif mode == "water_donor":
        dipole = -radial
        rotation = rotation_from_x(dipole)
    elif mode == "tangential":
        trial = np.cross(radial, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(trial) < 1e-8:
            trial = np.cross(radial, np.array([0.0, 1.0, 0.0]))
        rotation = rotation_from_x(trial)
    elif mode == "random":
        rotation = random_rotation(rng)
    else:
        raise KeyError(mode)
    return WATER @ rotation.T


def place_at_minimum_distance(
    solute: np.ndarray,
    anchor: int,
    radial: np.ndarray,
    target: float,
    mode: str,
    rng: np.random.Generator,
) -> np.ndarray:
    radial = radial / np.linalg.norm(radial)
    oriented = orient_water(radial, mode, rng)

    def placed(offset: float) -> np.ndarray:
        return oriented + solute[anchor] + offset * radial

    def residual(offset: float) -> float:
        wat = placed(offset)
        return float(
            np.linalg.norm(solute[:, None] - wat[None, :], axis=2).min() - target
        )

    high = 18.0
    low = high
    while low > 0 and residual(low) > 0:
        low -= 0.025
    if residual(low) > 0:
        raise RuntimeError("Could not bracket water placement")
    for _ in range(70):
        mid = 0.5 * (low + high)
        if residual(mid) > 0:
            high = mid
        else:
            low = mid
    wat = placed(0.5 * (low + high))
    achieved = np.linalg.norm(solute[:, None] - wat[None, :], axis=2).min()
    if not np.isclose(achieved, target, atol=1e-7):
        raise AssertionError((achieved, target))
    return wat


def build_cluster(
    solute: np.ndarray,
    anchors: list[int],
    n_waters: int,
    rng: np.random.Generator,
    regime: str,
) -> list[np.ndarray]:
    waters: list[np.ndarray] = []
    for wi in range(n_waters):
        accepted = None
        for _ in range(200):
            anchor = anchors[(wi + int(rng.integers(len(anchors)))) % len(anchors)]
            radial = direction_outward(solute, anchor, rng)
            if regime == "equilibrium":
                distance = float(rng.uniform(1.85, 2.75))
            else:
                distance = float(rng.uniform(1.45, 1.75))
            candidate = place_at_minimum_distance(
                solute, anchor, radial, distance, "random", rng
            )
            if waters:
                oo = [np.linalg.norm(candidate[0] - water[0]) for water in waters]
                all_dist = [
                    np.linalg.norm(candidate[:, None] - water[None, :], axis=2).min()
                    for water in waters
                ]
                if min(oo) < 2.25 or min(all_dist) < 1.25:
                    continue
            accepted = candidate
            break
        if accepted is None:
            raise RuntimeError(f"Could not place water {wi} in {regime} cluster")
        waters.append(accepted)
    return waters


def make_atoms(
    mol: Chem.Mol,
    solute: np.ndarray,
    waters: list[np.ndarray],
    metadata: dict,
) -> Atoms:
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = np.vstack([solute] + waters) if waters else solute.copy()
    atoms = Atoms(
        symbols + [x for _ in waters for x in ("O", "H", "H")], positions=positions
    )
    atoms.info.update(metadata)
    atoms.info["n_solute_atoms"] = len(symbols)
    atoms.info["n_waters"] = len(waters)
    atoms.arrays["fragment_id"] = np.array(
        [0] * len(symbols) + [i + 1 for i in range(len(waters)) for _ in range(3)]
    )
    return atoms


def main() -> None:
    output = ROOT / "data/configurations"
    output.mkdir(parents=True, exist_ok=True)
    molecules = pd.read_csv(ROOT / "config/molecules.csv")
    all_atoms: list[Atoms] = []
    rows: list[dict] = []

    for mi, row in molecules.iterrows():
        rng = np.random.default_rng(SEED + mi * 1009)
        mol, confs = conformers(row.smiles, SEED + mi)
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        hetero = [
            a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (7, 8, 9, 16, 17)
        ]
        anchors = hetero + [idx for idx in heavy if idx not in hetero]
        base = confs[0]

        def append(
            category: str,
            suffix: str,
            solute: np.ndarray,
            waters: list[np.ndarray],
            **extra,
        ):
            config_id = f"{row.split}__{row.molecule_id}__{category}__{suffix}"
            metadata = {
                "config_id": config_id,
                "molecule_id": row.molecule_id,
                "smiles": row.smiles,
                "split": row.split,
                "chemical_role": row.chemical_role,
                "category": category,
                **extra,
            }
            atoms = make_atoms(mol, solute, waters, metadata)
            nsolute = len(base)
            mind = float(
                np.linalg.norm(
                    atoms.positions[:nsolute, None] - atoms.positions[None, nsolute:],
                    axis=2,
                ).min()
            )
            metadata["minimum_interfragment_distance_angstrom"] = mind
            all_atoms.append(atoms)
            rows.append(
                {
                    **metadata,
                    "n_atoms": len(atoms),
                    "n_solute_atoms": nsolute,
                    "n_waters": len(waters),
                }
            )

        scan_modes = ["water_donor"] + (["tangential"] if row.split == "test" else [])
        for mode in scan_modes:
            anchor = anchors[0]
            radial = direction_outward(base, anchor, rng)
            for di, distance in enumerate((1.45, 1.70, 2.00, 2.40, 3.00, 4.00, 6.00)):
                wat = place_at_minimum_distance(
                    base, anchor, radial, distance, mode, rng
                )
                regime = (
                    "repulsive"
                    if distance < 1.8
                    else ("equilibrium" if distance <= 3.0 else "asymptotic")
                )
                append(
                    "scan",
                    f"{mode}_r{di}",
                    base,
                    [wat],
                    scan_mode=mode,
                    scan_distance_angstrom=distance,
                    regime=regime,
                    conformer_index=0,
                )

        for di in range(3):
            coords = confs[di % len(confs)]
            anchor = anchors[di % len(anchors)]
            radial = direction_outward(coords, anchor, rng)
            distance = float((2.05, 2.45, 2.85)[di])
            mode = ("water_donor", "water_acceptor", "random")[di]
            wat = place_at_minimum_distance(coords, anchor, radial, distance, mode, rng)
            append(
                "dimer",
                f"eq{di}",
                coords,
                [wat],
                scan_mode=mode,
                regime="equilibrium",
                conformer_index=di % len(confs),
            )

        for nw in (2, 4):
            for ci in range(2):
                coords = confs[ci % len(confs)]
                waters = build_cluster(coords, anchors, nw, rng, "equilibrium")
                append(
                    f"cluster_{nw}water",
                    f"env{ci}",
                    coords,
                    waters,
                    regime="equilibrium",
                    conformer_index=ci % len(confs),
                )

    all_atoms = sorted(all_atoms, key=lambda x: x.info["config_id"])
    registry = pd.DataFrame(rows).sort_values("config_id").reset_index(drop=True)
    assert registry.config_id.tolist() == [x.info["config_id"] for x in all_atoms]
    registry.to_csv(output / "configuration_registry.csv", index=False)
    write(output / "configurations_all.extxyz", all_atoms, format="extxyz")
    for split in ("train", "validation", "test"):
        write(
            output / f"configurations_{split}.extxyz",
            [x for x in all_atoms if x.info["split"] == split],
            format="extxyz",
        )
    manifest = {
        "seed": SEED,
        "n_configurations": len(registry),
        "counts_by_split": registry.groupby("split").size().astype(int).to_dict(),
        "counts_by_category": registry.groupby("category").size().astype(int).to_dict(),
        "molecules_by_split": {
            k: list(v)
            for k, v in registry.groupby("split").molecule_id.unique().items()
        },
        "files": {p.name: sha256(p) for p in sorted(output.glob("*")) if p.is_file()},
        "experimental_hydration_targets_accessed": False,
    }
    (output / "configuration_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
