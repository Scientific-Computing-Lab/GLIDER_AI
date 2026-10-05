"""Construction and serialization of the frozen H0 Hamiltonian."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import openmm
from openff.toolkit import Molecule
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator


TEMPERATURE = 298.15 * unit.kelvin
PRESSURE = 1.0 * unit.bar
CUTOFF = 1.0 * unit.nanometer
SWITCH = 0.9 * unit.nanometer
PADDING = 1.8 * unit.nanometer
TIMESTEP = 1.0 * unit.femtosecond
FORCEFIELD = "openff-2.2.1.offxml"


@dataclass
class H0System:
    molecule: Molecule
    topology: app.Topology
    positions: unit.Quantity
    aqueous_system: openmm.System
    gas_system: openmm.System
    n_solute_atoms: int


def ensure_ambertools_on_path() -> None:
    """OpenFF discovers AmberTools from PATH at import/registry construction."""
    prefix_bin = str(Path(__file__).resolve().parents[2] / ".conda_hlift" / "bin")
    if Path(prefix_bin, "antechamber").exists() and prefix_bin not in os.environ.get("PATH", "").split(":"):
        os.environ["PATH"] = prefix_bin + ":" + os.environ.get("PATH", "")


def _set_nonbonded_options(system: openmm.System) -> None:
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            if force.getNonbondedMethod() != openmm.NonbondedForce.NoCutoff:
                force.setUseSwitchingFunction(True)
                force.setSwitchingDistance(SWITCH)
                force.setUseDispersionCorrection(True)
                force.setEwaldErrorTolerance(1.0e-4)


def _parameterized_forcefield(molecule: Molecule, include_water: bool) -> app.ForceField:
    forcefield = app.ForceField("tip3p.xml") if include_water else app.ForceField()
    generator = SMIRNOFFTemplateGenerator(molecules=[molecule], forcefield=FORCEFIELD)
    forcefield.registerTemplateGenerator(generator.generator)
    return forcefield


def build_h0(smiles: str, seed: int) -> H0System:
    ensure_ambertools_on_path()
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=False)
    molecule.generate_conformers(n_conformers=1)
    molecule.assign_partial_charges("am1bcc")
    if not np.all(np.isfinite(molecule.partial_charges.m)):
        raise RuntimeError("Non-finite AM1-BCC charge")
    if abs(float(molecule.partial_charges.sum().m) - float(molecule.total_charge.m)) > 1.0e-6:
        raise RuntimeError("AM1-BCC charge sum does not match formal charge")

    solute_topology = molecule.to_topology().to_openmm()
    solute_positions = molecule.conformers[0].to_openmm()
    n_solute_atoms = molecule.n_atoms

    aqueous_ff = _parameterized_forcefield(molecule, include_water=True)
    modeller = app.Modeller(solute_topology, solute_positions)
    modeller.addSolvent(
        aqueous_ff,
        model="tip3p",
        padding=PADDING,
        boxShape="dodecahedron",
    )
    aqueous_system = aqueous_ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=CUTOFF,
        constraints=app.HBonds,
        rigidWater=True,
        removeCMMotion=True,
    )
    _set_nonbonded_options(aqueous_system)
    aqueous_system.addForce(openmm.MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))

    gas_ff = _parameterized_forcefield(molecule, include_water=False)
    gas_system = gas_ff.createSystem(
        solute_topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        removeCMMotion=True,
    )
    _set_nonbonded_options(gas_system)

    return H0System(
        molecule=molecule,
        topology=modeller.topology,
        positions=modeller.positions,
        aqueous_system=aqueous_system,
        gas_system=gas_system,
        n_solute_atoms=n_solute_atoms,
    )


def serialize_h0(bundle: H0System, output_dir: Path, metadata: dict) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    aqueous_xml = output_dir / "aqueous_system.xml"
    gas_xml = output_dir / "gas_system.xml"
    topology_pdb = output_dir / "initial_solvated.pdb"
    solute_sdf = output_dir / "solute.sdf"
    charges_csv = output_dir / "am1bcc_charges.csv"

    aqueous_xml.write_text(openmm.XmlSerializer.serialize(bundle.aqueous_system))
    gas_xml.write_text(openmm.XmlSerializer.serialize(bundle.gas_system))
    with topology_pdb.open("w") as handle:
        app.PDBFile.writeFile(bundle.topology, bundle.positions, handle, keepIds=True)
    bundle.molecule.to_file(solute_sdf, file_format="sdf")
    with charges_csv.open("w") as handle:
        handle.write("atom_index,element,partial_charge_e\n")
        for atom, charge in zip(bundle.molecule.atoms, bundle.molecule.partial_charges.m):
            handle.write(f"{atom.molecule_atom_index},{atom.symbol},{float(charge):.12f}\n")

    paths = [aqueous_xml, gas_xml, topology_pdb, solute_sdf, charges_csv]
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    record = {
        **metadata,
        "n_solute_atoms": bundle.n_solute_atoms,
        "n_aqueous_atoms": bundle.aqueous_system.getNumParticles(),
        "formal_charge_e": float(bundle.molecule.total_charge.m),
        "am1bcc_charge_sum_e": float(bundle.molecule.partial_charges.sum().m),
        "files_sha256": hashes,
    }
    (output_dir / "system_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def cuda_platform() -> tuple[openmm.Platform, dict[str, str]]:
    platform = openmm.Platform.getPlatformByName("CUDA")
    return platform, {"Precision": "mixed", "DeterministicForces": "true"}
