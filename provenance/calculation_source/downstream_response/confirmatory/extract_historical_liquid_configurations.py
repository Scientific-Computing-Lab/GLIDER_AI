#!/usr/bin/env python3
"""Materialize the target-free confirmatory liquid configurations from Git history.

Only geometry/topology data are read.  No GLIDER, comparator, or QM response
prediction is opened by this program.  Each source is a completed OpenMM
equilibration state already present in the Transformato repository history.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import subprocess
import xml.etree.ElementTree as ET
from collections import deque
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parents[4]
WORKSPACE = ROOT.parent
SOURCE_REPOSITORY = WORKSPACE / "foundation/sources/repos/transformato"


@dataclass(frozen=True)
class Source:
    molecule_id: str
    family: str
    canonical_smiles: str
    commit: str
    topology_path: str
    state_path: str
    solute_sdf_path: str


# Frozen by a target-free workspace/Git-history audit.  The tautomeric duplicate
# of 2OJ9 and test/alias copies are deliberately not counted as new solutes.
SOURCES = (
    Source(
        "2OJ9-original",
        "2OJ9",
        "Cc1cc(-n2ccnc2)cc2nc(-c3c(NCc4ccccn4)cc[nH]c3=O)[nH]c12",
        "3c952fa3607e410feb0453e9bb6b7cdff6c45884",
        "data/2OJ9-original/waterbox/openmm/step3_input.pdb",
        "data/2OJ9-original/waterbox/openmm/step4_equilibration.rst",
        "data/2OJ9-original/waterbox/bmi/BMI.sdf",
    ),
    Source(
        "cdk2-1h1q",
        "cdk2",
        "c1ccc(Nc2nc(OCC3CCCCC3)c3nc[nH]c3n2)cc1",
        "c2c447ebe57226083ccbe4a2b4073ba36adcac95",
        "data/cdk2-1h1q/waterbox/openmm/step3_input.pdb",
        "data/cdk2-1h1q/waterbox/openmm/step4_equilibration.rst",
        "data/cdk2-1h1q/waterbox/unk/UNK.sdf",
    ),
    Source(
        "cdk2-28",
        "cdk2",
        "CNS(=O)(=O)c1ccc(Nc2nc(OCC3CCCCC3)c3nc[nH]c3n2)cc1",
        "c2c447ebe57226083ccbe4a2b4073ba36adcac95",
        "data/cdk2-28/waterbox/openmm/step3_input.pdb",
        "data/cdk2-28/waterbox/openmm/step4_equilibration.rst",
        "data/cdk2-28/28.sdf",
    ),
    Source(
        "ethane",
        "small_aliphatic",
        "CC",
        "9e33a91f55a15bd0d25790d42e91681bfcb66dfb",
        "data/ethane/waterbox/openmm/step3_charmm2omm.pdb",
        "data/ethane/waterbox/openmm/step4_equilibration.rst",
        "data/ethane/waterbox/lig/lig.sdf",
    ),
    Source(
        "ethanol",
        "small_aliphatic",
        "CCO",
        "9e33a91f55a15bd0d25790d42e91681bfcb66dfb",
        "data/ethanol/waterbox/openmm/step3_charmm2omm.pdb",
        "data/ethanol/waterbox/openmm/step4_equilibration.rst",
        "data/ethanol/waterbox/unl/UNL.sdf",
    ),
    Source(
        "jnk1-17124",
        "jnk1",
        "CCOc1nc(NC(=O)Cc2cc(OC)c(Br)cc2OC)cc(N)c1C#N",
        "cf83b1ce0d272f37ee42b09ef3e5b3ef987a0a08",
        "data/jnk1-17124/waterbox/openmm/step3_input.pdb",
        "data/jnk1-17124/waterbox/openmm/step4_equilibration.rst",
        "data/jnk1-17124/ligand_17124.sdf",
    ),
    Source(
        "jnk1-18631",
        "jnk1",
        "CCOc1nc(NC(=O)Cc2ccccc2OC)cc(N)c1C#N",
        "cf83b1ce0d272f37ee42b09ef3e5b3ef987a0a08",
        "data/jnk1-18631/waterbox/openmm/step3_input.pdb",
        "data/jnk1-18631/waterbox/openmm/step4_equilibration.rst",
        "data/jnk1-18631/ligand_18631.sdf",
    ),
    Source(
        "methane",
        "small_aliphatic",
        "C",
        "34bf648301496897e676e7a285f756e30280d247",
        "data/methane/waterbox/openmm/step3_charmm2omm.pdb",
        "data/methane/waterbox/openmm/step4_equilibration.rst",
        "data/methane/methane.sdf",
    ),
    Source(
        "tyk2-ejm_42",
        "tyk2",
        "CCC(=O)Nc1cc(NC(=O)c2c(Cl)cccc2Cl)ccn1",
        "2e0666db35b4dec53175ea55fa6a5aacd579bf93",
        "data/tyk2-ejm_42/waterbox/openmm/step3_input.pdb",
        "data/tyk2-ejm_42/waterbox/openmm/step4_equilibration.rst",
        "data/tyk2-ejm_42/ejm_42.sdf",
    ),
    Source(
        "tyk2-ejm_45",
        "tyk2",
        "O=C(CC1CC1)Nc1cc(NC(=O)c2c(Cl)cccc2Cl)ccn1",
        "2e0666db35b4dec53175ea55fa6a5aacd579bf93",
        "data/tyk2-ejm_45/waterbox/openmm/step3_input.pdb",
        "data/tyk2-ejm_45/waterbox/openmm/step4_equilibration.rst",
        "data/tyk2-ejm_45/ejm_45.sdf",
    ),
)

ALLOWED_ELEMENTS = {"H", "C", "N", "O", "F", "S", "Cl", "Br", "I"}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.check_output(
        ["git", "-C", str(SOURCE_REPOSITORY), "show", f"{commit}:{path}"]
    )


def parse_sdf(value: bytes) -> tuple[list[str], list[tuple[int, int]]]:
    lines = value.decode().splitlines()
    if len(lines) < 4:
        raise RuntimeError("truncated SDF")
    atom_count = int(lines[3][0:3])
    bond_count = int(lines[3][3:6])
    elements = [lines[4 + index][31:34].strip() for index in range(atom_count)]
    bonds = []
    for line in lines[4 + atom_count : 4 + atom_count + bond_count]:
        bonds.append((int(line[0:3]) - 1, int(line[3:6]) - 1))
    return elements, bonds


def parse_pdb(value: bytes) -> list[dict[str, object]]:
    atoms = []
    for line in io.StringIO(value.decode()):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            {
                "atom_name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22],
                "resid": line[22:26].strip(),
                "element": line[76:78].strip(),
            }
        )
    return atoms


def parse_state(value: bytes) -> tuple[np.ndarray, np.ndarray, float, str]:
    root = ET.fromstring(value)
    positions = np.asarray(
        [
            [float(item.attrib[axis]) for axis in ("x", "y", "z")]
            for item in root.find("Positions")
        ]
    )
    box_node = root.find("PeriodicBoxVectors")
    box = np.asarray(
        [
            [float(item.attrib[axis]) for axis in ("x", "y", "z")]
            for item in box_node
        ]
    )
    return positions * 10.0, box * 10.0, float(root.attrib["time"]), root.attrib["openmmVersion"]


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    fractional = np.linalg.solve(box.T, np.asarray(delta).T).T
    fractional -= np.rint(fractional)
    return fractional @ box


def unwrap_solute(
    positions: np.ndarray, bonds: list[tuple[int, int]], box: np.ndarray
) -> np.ndarray:
    adjacency: list[list[int]] = [[] for _ in positions]
    for first, second in bonds:
        adjacency[first].append(second)
        adjacency[second].append(first)
    result = np.full_like(positions, np.nan)
    result[0] = positions[0]
    queue = deque([0])
    while queue:
        first = queue.popleft()
        for second in adjacency[first]:
            if np.all(np.isfinite(result[second])):
                continue
            result[second] = result[first] + minimum_image(
                positions[second] - positions[first], box
            )
            queue.append(second)
    if not np.all(np.isfinite(result)):
        raise RuntimeError("SDF solute graph is disconnected")
    return result


def water_groups(atoms: list[dict[str, object]]) -> list[tuple[int, list[int]]]:
    grouped: dict[tuple[str, str], list[int]] = {}
    for index, atom in enumerate(atoms):
        if atom["resname"] != "TIP":
            continue
        key = (str(atom["chain"]), str(atom["resid"]))
        grouped.setdefault(key, []).append(index)
    result = []
    for source_index, indices in enumerate(grouped.values()):
        if len(indices) != 3:
            raise RuntimeError("expected every TIP residue to contain O-H-H")
        names = [str(atoms[index]["atom_name"]) for index in indices]
        oxygen = names.index("OH2") if "OH2" in names else 0
        ordered = [indices[oxygen], *[index for i, index in enumerate(indices) if i != oxygen]]
        result.append((source_index, ordered))
    return result


def materialize(source: Source):
    topology_bytes = git_bytes(source.commit, source.topology_path)
    state_bytes = git_bytes(source.commit, source.state_path)
    sdf_bytes = git_bytes(source.commit, source.solute_sdf_path)
    topology = parse_pdb(topology_bytes)
    positions, box, time_ps, openmm_version = parse_state(state_bytes)
    if len(topology) != len(positions):
        raise RuntimeError(f"topology/state atom mismatch: {source.molecule_id}")
    sdf_elements, bonds = parse_sdf(sdf_bytes)

    solute_topology_indices = []
    solute_elements = []
    for index, atom in enumerate(topology):
        if atom["resname"] in {"TIP", "POT", "CLA"}:
            continue
        if str(atom["element"]).upper() == "DU" or str(atom["atom_name"]).startswith("LP"):
            continue
        solute_topology_indices.append(index)
        solute_elements.append(str(atom["element"]))
    for index, element in enumerate(solute_elements):
        if not element:
            name = str(topology[solute_topology_indices[index]]["atom_name"])
            solute_elements[index] = "Cl" if name.upper().startswith("CL") else name[0].upper()
    if solute_elements != sdf_elements:
        raise RuntimeError(f"PDB/SDF element ordering differs: {source.molecule_id}")
    if not set(solute_elements).issubset(ALLOWED_ELEMENTS):
        raise RuntimeError(f"unsupported element in {source.molecule_id}")

    raw_solute = positions[solute_topology_indices]
    solute = unwrap_solute(raw_solute, bonds, box)
    ranked = []
    for source_index, indices in water_groups(topology):
        water = positions[indices]
        displacements = minimum_image(water[0] - raw_solute, box)
        nearest = int(np.argmin(np.linalg.norm(displacements, axis=1)))
        oxygen = solute[nearest] + displacements[nearest]
        hydrogen_1 = oxygen + minimum_image(water[1] - water[0], box)
        hydrogen_2 = oxygen + minimum_image(water[2] - water[0], box)
        distance = float(np.linalg.norm(displacements[nearest]))
        ranked.append((distance, source_index, np.vstack((oxygen, hydrogen_1, hydrogen_2))))
    ranked.sort(key=lambda item: (item[0], item[1]))
    if len(ranked) < 12:
        raise RuntimeError(f"fewer than 12 waters: {source.molecule_id}")
    return {
        "topology_bytes": topology_bytes,
        "state_bytes": state_bytes,
        "sdf_bytes": sdf_bytes,
        "solute_elements": solute_elements,
        "solute_positions": solute,
        "ranked": ranked,
        "box": box,
        "time_ps": time_ps,
        "openmm_version": openmm_version,
        "n_parent_atoms": len(topology),
        "n_parent_waters": len(ranked),
        "n_parent_ions": sum(atom["resname"] in {"POT", "CLA"} for atom in topology),
        "n_virtual_sites_excluded": sum(
            str(atom["element"]).upper() == "DU" or str(atom["atom_name"]).startswith("LP")
            for atom in topology
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results/posthoc_downstream_response_confirmatory/configurations",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    bases = []
    w4_probes = []
    outer_probes = []
    base_plus_w4 = []
    case_rows = []
    probe_rows = []
    source_records = []
    for source in SOURCES:
        values = materialize(source)
        solute = Atoms(values["solute_elements"], positions=values["solute_positions"])
        first_twelve = values["ranked"][:12]
        base = solute.copy()
        for _, _, coordinates in first_twelve[:3]:
            base += Atoms(("O", "H", "H"), positions=coordinates)
        config_id = f"confirmatory__{source.molecule_id}__base3"
        base.info = {
            "config_id": config_id,
            "case_id": config_id,
            "molecule_id": source.molecule_id,
            "family": source.family,
            "regime": "equilibrated_explicit_water",
            "n_solute_atoms": len(solute),
            "n_waters": 3,
            "source_frame_index": 0,
        }
        bases.append(base)
        for rank, (distance, source_index, coordinates) in enumerate(first_twelve, 1):
            if rank <= 3:
                continue
            probe_id = f"{config_id}__W{rank}"
            probe = Atoms(("O", "H", "H"), positions=coordinates)
            probe.info = {
                "probe_id": probe_id,
                "case_id": config_id,
                "base_config_id": config_id,
                "molecule_id": source.molecule_id,
                "water_rank": rank,
                "source_water_index": source_index,
                "oxygen_distance_A": distance,
            }
            outer_probes.append(probe)
            if rank == 4:
                w4_probes.append(probe.copy())
                joined = base + probe
                joined.info = dict(base.info)
                joined.info.update({"probe_id": probe_id, "n_base_atoms": len(base)})
                base_plus_w4.append(joined)
            probe_rows.append(
                {
                    "probe_id": probe_id,
                    "case_id": config_id,
                    "molecule_id": source.molecule_id,
                    "water_rank": rank,
                    "source_water_index": source_index,
                    "oxygen_distance_A": distance,
                }
            )

        case_row = {
            "case_id": config_id,
            "molecule_id": source.molecule_id,
            "family": source.family,
            "canonical_smiles": source.canonical_smiles,
            "n_solute_atoms": len(solute),
            "n_base_atoms": len(base),
            "source_repository": str(SOURCE_REPOSITORY),
            "source_commit": source.commit,
            "source_topology_path": source.topology_path,
            "source_topology_sha256": sha256_bytes(values["topology_bytes"]),
            "source_state_path": source.state_path,
            "source_state_sha256": sha256_bytes(values["state_bytes"]),
            "source_solute_sdf_path": source.solute_sdf_path,
            "source_solute_sdf_sha256": sha256_bytes(values["sdf_bytes"]),
            "equilibration_state_time_ps": values["time_ps"],
            "openmm_version": values["openmm_version"],
            "n_parent_atoms": values["n_parent_atoms"],
            "n_parent_waters": values["n_parent_waters"],
            "n_parent_ions_excluded": values["n_parent_ions"],
            "n_virtual_sites_excluded": values["n_virtual_sites_excluded"],
            "response_supervision_identity_match": False,
            "selection_used_model_or_qm_performance": False,
        }
        for rank, (distance, source_index, _) in enumerate(first_twelve, 1):
            case_row[f"w{rank}_source_index"] = source_index
            case_row[f"w{rank}_distance_A"] = distance
        case_rows.append(case_row)
        source_records.append({**asdict(source), **{k: case_row[k] for k in case_row if k.endswith("sha256")}})

    paths = {
        "base_3water.extxyz": bases,
        "heldout_w4.extxyz": w4_probes,
        "outer_w4_w12.extxyz": outer_probes,
        "base_plus_w4.extxyz": base_plus_w4,
    }
    for name, frames in paths.items():
        write(args.output / name, frames, format="extxyz")
    pd.DataFrame(case_rows).sort_values("molecule_id").to_csv(
        args.output / "configurations.csv", index=False
    )
    pd.DataFrame(probe_rows).sort_values(["molecule_id", "water_rank"]).to_csv(
        args.output / "outer_probes.csv", index=False
    )
    with (args.output / "source_records.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_records[0].keys())
        writer.writeheader()
        writer.writerows(source_records)

    manifest = {
        "experiment": "post hoc held-out-water confirmatory breadth and spatial range",
        "source_repository": str(SOURCE_REPOSITORY),
        "n_solutes": len(SOURCES),
        "n_families": len({source.family for source in SOURCES}),
        "snapshots_per_solute": 1,
        "n_outer_probes": len(outer_probes),
        "water_ranks": list(range(4, 13)),
        "base_definition": "solute plus W1-W3",
        "ordering_rule": "minimum-image oxygen distance to nearest solute atom, then source water index",
        "source_rule": "all unique eligible completed explicit-water equilibration states found in the workspace/Git audit; one 2OJ9 identity retained after alias/tautomer collapse",
        "outputs": {
            path.name: sha256(path)
            for path in sorted(args.output.iterdir())
            if path.is_file() and path.name != "extraction_manifest.json"
        },
        "performance_based_selection": False,
        "reference_or_prediction_values_accessed": False,
        "prospective_claim": False,
    }
    (args.output / "extraction_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
