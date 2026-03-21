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
from .molecules import (
    STEARATE_NUM_ATOMS,
    SODIUM_NUM_ATOMS,
    WATER_NUM_ATOMS,
    create_stearate_off,
    create_sodium_off,
)

log = logging.getLogger(__name__)


# ── Internal helpers ──────────────────────────────────────────────────


def _build_water_topology(num_water: int) -> tuple[app.Topology, np.ndarray]:
    """Return a clean water-only OpenMM Topology (O, H1, H2 per residue).

    Also returns a boolean mask (not used here but kept for symmetry).
    """
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


def _rename_residues(topology: app.Topology,
                     n_stearate: int, n_sodium: int, n_water: int) -> None:
    """Ensure every residue has the name the SystemGenerator expects."""
    idx = 0
    for chain in topology.chains():
        for res in chain.residues():
            if idx < n_stearate:
                res.name = "STL"
            elif idx < n_stearate + n_sodium:
                res.name = "NA"
            else:
                res.name = "HOH"
            idx += 1


# ── Public API ────────────────────────────────────────────────────────


def parameterize_system(packed_pdb: Path, config: Config):
    """Read PACKMOL output, assign force-field parameters, return OpenMM objects.

    Returns
    -------
    system : openmm.System
    topology : openmm.app.Topology
    positions : openmm.unit.Quantity
    """
    sys_cfg = config.system
    ff_cfg = config.forcefield

    # ── 1. Load raw positions from PACKMOL PDB ───────────────────────
    log.info("Loading packed coordinates from %s", packed_pdb)
    pdb = app.PDBFile(str(packed_pdb))
    all_pos_nm = np.array(pdb.getPositions().value_in_unit(unit.nanometers))

    # ── 2. Slice positions by known PACKMOL ordering ─────────────────
    #   stearate(1..N) | sodium(1..N) | water(1..N)
    ns = sys_cfg.num_stearate
    nw = sys_cfg.num_water
    n_stearate_atoms = ns * STEARATE_NUM_ATOMS
    n_sodium_atoms = ns * SODIUM_NUM_ATOMS
    n_solute_atoms = n_stearate_atoms + n_sodium_atoms
    n_water_atoms = nw * WATER_NUM_ATOMS

    solute_pos = all_pos_nm[:n_solute_atoms] * unit.nanometers
    water_pos = all_pos_nm[n_solute_atoms:n_solute_atoms + n_water_atoms] * unit.nanometers

    log.info("Atoms: %d stearate + %d sodium + %d water = %d total",
             n_stearate_atoms, n_sodium_atoms, n_water_atoms,
             n_stearate_atoms + n_sodium_atoms + n_water_atoms)

    # ── 3. Build solute topology via OpenFF ───────────────────────────
    log.info("Building solute topology with OpenFF ...")
    stearate_mol = create_stearate_off()
    sodium_mol = create_sodium_off()
    stearate_mol.name = "STL"
    sodium_mol.name = "NA"

    solute_mols = [stearate_mol] * ns + [sodium_mol] * ns
    off_topo = Topology.from_molecules(solute_mols)
    omm_solute_topo = off_topo.to_openmm()

    # ── 4. Combine solute + water via Modeller ────────────────────────
    modeller = app.Modeller(omm_solute_topo, solute_pos)
    water_topo = _build_water_topology(nw)
    modeller.add(water_topo, water_pos)

    _rename_residues(modeller.topology, ns, ns, nw)
    log.info("Modeller topology: %d atoms, %d residues",
             modeller.topology.getNumAtoms(),
             sum(1 for _ in modeller.topology.residues()))

    # ── 5. SystemGenerator ────────────────────────────────────────────
    log.info("Creating SystemGenerator  (solute: %s, water: %s)",
             ff_cfg.small_molecule, ff_cfg.water_model)

    generator = SystemGenerator(
        forcefields=[ff_cfg.water_model],
        small_molecule_forcefield=ff_cfg.small_molecule,
        molecules=[stearate_mol, sodium_mol],
        forcefield_kwargs={"constraints": app.HBonds},
        periodic_forcefield_kwargs={
            "nonbondedMethod": app.PME,
            "nonbondedCutoff": 1.0 * unit.nanometers,
        },
    )

    log.info("Generating OpenMM System (this may take a moment) ...")
    system = generator.create_system(modeller.topology)

    # ── 6. Periodic box vectors ───────────────────────────────────────
    box_nm = sys_cfg.box_size_angstrom * 0.1  # A -> nm
    system.setDefaultPeriodicBoxVectors(
        mm.Vec3(box_nm, 0, 0) * unit.nanometers,
        mm.Vec3(0, box_nm, 0) * unit.nanometers,
        mm.Vec3(0, 0, box_nm) * unit.nanometers,
    )

    # ── 7. Quick sanity checks ────────────────────────────────────────
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            q_total = sum(
                force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                for i in range(force.getNumParticles())
            )
            log.info("Net charge: %.4f e  |  Constraints: %d",
                     q_total, system.getNumConstraints())
            break

    return system, modeller.topology, modeller.positions


# ── Serialization ─────────────────────────────────────────────────────


def save_system(system: mm.System, topology: app.Topology,
                positions, output_dir: Path) -> None:
    """Persist the parameterized system to disk."""
    param_dir = output_dir / "parameterize"
    param_dir.mkdir(parents=True, exist_ok=True)

    xml_path = param_dir / "system.xml"
    xml_path.write_text(mm.XmlSerializer.serialize(system))

    pdb_path = param_dir / "topology.pdb"
    with pdb_path.open("w") as fh:
        app.PDBFile.writeFile(topology, positions, fh)

    log.info("Saved system.xml + topology.pdb -> %s", param_dir)


def load_system(output_dir: Path):
    """Load a previously-saved parameterized system.

    Returns (system, topology, positions).
    """
    param_dir = output_dir / "parameterize"
    system = mm.XmlSerializer.deserialize(
        (param_dir / "system.xml").read_text()
    )
    pdb = app.PDBFile(str(param_dir / "topology.pdb"))
    return system, pdb.topology, pdb.positions
