"""Force-field parameterization with OpenFF + OpenMM SystemGenerator."""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openff.toolkit.topology import Topology
from openmmforcefields.generators import SystemGenerator

from .config import Config
from .molecules import atom_count, create_off_molecule

log = logging.getLogger(__name__)

WATER_NUM_ATOMS = 3


# ── Internal helpers ──────────────────────────────────────────────────


def _build_water_topology(num_water: int) -> app.Topology:
    topo = app.Topology()
    chain = topo.addChain()
    oxygen = app.Element.getBySymbol("O")
    hydrogen = app.Element.getBySymbol("H")
    for i in range(num_water):
        res = topo.addResidue("HOH", chain, str(i + 1))
        o = topo.addAtom("O", oxygen, res)
        h1 = topo.addAtom("H1", hydrogen, res)
        h2 = topo.addAtom("H2", hydrogen, res)
        topo.addBond(o, h1)
        topo.addBond(o, h2)
    return topo


def _rename_residues(topology: app.Topology, labels: list[str]) -> None:
    """Set residue names from an ordered label list."""
    for idx, (_, res) in enumerate(
        (ch, r) for ch in topology.chains() for r in ch.residues()
    ):
        if idx < len(labels):
            res.name = labels[idx]


# ── Public API ────────────────────────────────────────────────────────


def parameterize_system(packed_pdb: Path, config: Config):
    """Read PACKMOL output, assign FF parameters, return OpenMM objects.

    Returns ``(system, topology, positions)``.
    """
    sys_cfg = config.system
    ff_cfg = config.forcefield

    # ── 1. Load raw positions ─────────────────────────────────────────
    log.info("Loading packed coordinates from %s", packed_pdb)
    pdb = app.PDBFile(str(packed_pdb))
    all_pos_nm = np.array(pdb.getPositions().value_in_unit(unit.nanometers))

    # ── 2. Slice positions by PACKMOL ordering ────────────────────────
    #   solute1×n1 | solute2×n2 | … | counterion_type1 | counterion_type2 | … | water
    offset = 0
    for spec in sys_cfg.solutes:
        offset += spec.count * atom_count(spec.smiles)
    for ci, n in sys_cfg.counterion_counts:
        offset += n * atom_count(ci.smiles)

    n_solute_total = offset
    n_water_atoms = sys_cfg.num_water * WATER_NUM_ATOMS

    solute_pos = all_pos_nm[:n_solute_total] * unit.nanometers
    water_pos = all_pos_nm[n_solute_total:n_solute_total + n_water_atoms] * unit.nanometers

    log.info("Solute atoms: %d  |  Water atoms: %d  |  Total: %d",
             n_solute_total, n_water_atoms, n_solute_total + n_water_atoms)

    # ── 3. Build solute topology via OpenFF ───────────────────────────
    log.info("Building solute topology with OpenFF ...")
    off_molecules = []   # unique molecules for SystemGenerator
    solute_mol_list = []  # ordered list matching PACKMOL layout
    residue_labels = []  # for renaming

    seen_smiles: dict[str, object] = {}
    for spec in sys_cfg.solutes:
        if spec.smiles not in seen_smiles:
            mol = create_off_molecule(spec.smiles)
            mol.name = spec.resname
            seen_smiles[spec.smiles] = mol
            off_molecules.append(mol)
        solute_mol_list.extend([seen_smiles[spec.smiles]] * spec.count)
        residue_labels.extend([spec.resname] * spec.count)

    # Counterions (one or more types)
    for ci, n in sys_cfg.counterion_counts:
        if ci.smiles not in seen_smiles:
            mol = create_off_molecule(ci.smiles)
            mol.name = ci.resname
            seen_smiles[ci.smiles] = mol
            off_molecules.append(mol)
        solute_mol_list.extend([seen_smiles[ci.smiles]] * n)
        residue_labels.extend([ci.resname] * n)
        log.info("  counterion %s: %d ions", ci.resname, n)

    off_topo = Topology.from_molecules(solute_mol_list)
    omm_solute_topo = off_topo.to_openmm()

    # ── 4. Combine solute + water via Modeller ────────────────────────
    modeller = app.Modeller(omm_solute_topo, solute_pos)
    water_topo = _build_water_topology(sys_cfg.num_water)
    modeller.add(water_topo, water_pos)

    residue_labels.extend(["HOH"] * sys_cfg.num_water)
    _rename_residues(modeller.topology, residue_labels)

    # ── 5. Set box vectors (must precede create_system for PME) ───────
    bx, by, bz = [v * 0.1 for v in sys_cfg.box_angstrom]
    modeller.topology.setPeriodicBoxVectors((
        mm.Vec3(bx, 0, 0), mm.Vec3(0, by, 0), mm.Vec3(0, 0, bz),
    ))

    log.info("Topology: %d atoms, %d residues, box %.2f x %.2f x %.2f nm",
             modeller.topology.getNumAtoms(),
             sum(1 for _ in modeller.topology.residues()), bx, by, bz)

    # ── 6. SystemGenerator ────────────────────────────────────────────
    log.info("Creating SystemGenerator  (%d unique molecules, water: %s)",
             len(off_molecules), ff_cfg.water_model)

    generator = SystemGenerator(
        forcefields=[ff_cfg.water_model],
        small_molecule_forcefield=ff_cfg.small_molecule,
        molecules=off_molecules,
        forcefield_kwargs={"constraints": app.HBonds},
        periodic_forcefield_kwargs={
            "nonbondedMethod": app.PME,
            "nonbondedCutoff": 1.0 * unit.nanometers,
        },
    )

    log.info("Generating OpenMM System ...")
    system = generator.create_system(modeller.topology)

    # ── 7. Sanity checks ──────────────────────────────────────────────
    if not system.usesPeriodicBoundaryConditions():
        raise RuntimeError("System is not periodic -- PME was not applied.")

    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            q = sum(force.getParticleParameters(i)[0].value_in_unit(
                    unit.elementary_charge) for i in range(force.getNumParticles()))
            log.info("Net charge: %.4f e  |  Constraints: %d",
                     q, system.getNumConstraints())
            break

    return system, modeller.topology, modeller.positions


# ── Serialization ─────────────────────────────────────────────────────


def save_system(system, topology, positions, output_dir: Path) -> None:
    param_dir = output_dir / "parameterize"
    param_dir.mkdir(parents=True, exist_ok=True)
    (param_dir / "system.xml").write_text(mm.XmlSerializer.serialize(system))
    with (param_dir / "topology.pdb").open("w") as fh:
        app.PDBFile.writeFile(topology, positions, fh)
    log.info("Saved system.xml + topology.pdb -> %s", param_dir)


def load_system(output_dir: Path):
    param_dir = output_dir / "parameterize"
    system = mm.XmlSerializer.deserialize(
        (param_dir / "system.xml").read_text())
    pdb = app.PDBFile(str(param_dir / "topology.pdb"))
    return system, pdb.topology, pdb.positions
