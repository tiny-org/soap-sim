"""Molecule generation with RDKit and OpenFF interop.

All functions accept arbitrary SMILES so the system is not limited
to stearate -- any combination of surfactants can be simulated.
"""
from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit.topology import Molecule

log = logging.getLogger(__name__)

_RDKIT_SEED = 42


# ── Generic helpers ───────────────────────────────────────────────────


@lru_cache(maxsize=32)
def atom_count(smiles: str) -> int:
    """Number of atoms (including explicit H) for a SMILES string."""
    return Chem.AddHs(Chem.MolFromSmiles(smiles)).GetNumAtoms()


@lru_cache(maxsize=32)
def _embed(smiles: str) -> Chem.Mol:
    """RDKit Mol with 3-D coordinates for any SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    if mol.GetNumAtoms() == 1:
        # Single atom (ion) -- just place at origin
        conf = Chem.Conformer(1)
        mol.AddConformer(conf, assignId=True)
        return mol
    params = AllChem.ETKDGv3()
    params.randomSeed = _RDKIT_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"RDKit embedding failed for {smiles}")
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def create_off_molecule(smiles: str) -> Molecule:
    """OpenFF Molecule with atom ordering matching :func:`generate_pdb`."""
    return Molecule.from_rdkit(_embed(smiles), allow_undefined_stereo=True)


# ── PDB writers ───────────────────────────────────────────────────────


def generate_solute_pdb(smiles: str, resname: str, path: Path) -> None:
    """Write a 3-D monomer PDB for any SMILES."""
    mol = _embed(smiles)
    Chem.MolToPDBFile(mol, str(path))
    text = path.read_text()
    path.write_text(text.replace("UNL", resname[:3].upper()))
    log.info("Wrote %s -> %s  (%d atoms)", resname, path, mol.GetNumAtoms())


def generate_counterion_pdb(smiles: str, path: Path) -> None:
    """Write a PDB for a single-atom ion (Na+, K+, Cl-, ...)."""
    mol = Chem.MolFromSmiles(smiles)
    elem = mol.GetAtomWithIdx(0).GetSymbol()
    resname = elem[:3].upper()
    # Fixed-width PDB ATOM record
    path.write_text(
        f"ATOM      1 {elem:>3s}  {resname:<3s} A   1"
        f"      0.000   0.000   0.000  1.00  0.00"
        f"          {elem:>2s}\n"
        f"END\n"
    )
    log.info("Wrote counterion (%s) -> %s", elem, path)


def generate_water_pdb(path: Path) -> None:
    """TIP3P water geometry: O-H = 0.9572 A, H-O-H = 104.52 deg."""
    path.write_text(
        "ATOM      1  O   HOH A   1"
        "      0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      2  H1  HOH A   1"
        "      0.757   0.000   0.586  1.00  0.00           H\n"
        "ATOM      3  H2  HOH A   1"
        "     -0.757   0.000   0.586  1.00  0.00           H\n"
        "END\n"
    )
    log.info("Wrote TIP3P water -> %s", path)


def generate_all_monomers(config, output_dir: Path) -> dict[str, Path]:
    """Write all monomer PDBs for the configured system."""
    output_dir.mkdir(parents=True, exist_ok=True)
    sys_cfg = config.system
    paths: dict[str, Path] = {}

    for spec in sys_cfg.solutes:
        p = output_dir / f"{spec.name}.pdb"
        generate_solute_pdb(spec.smiles, spec.resname, p)
        paths[spec.name] = p

    for i, (ci, count) in enumerate(sys_cfg.counterion_counts):
        p = output_dir / f"counterion_{ci.resname.lower()}.pdb"
        generate_counterion_pdb(ci.smiles, p)
        paths[f"counterion_{ci.resname}"] = p

    p = output_dir / "water.pdb"
    generate_water_pdb(p)
    paths["water"] = p

    return paths
