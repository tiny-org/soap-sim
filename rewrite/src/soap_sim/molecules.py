"""Molecule generation with RDKit and OpenFF interop."""
from __future__ import annotations

import logging
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit.topology import Molecule

log = logging.getLogger(__name__)

# ── SMILES & constants ────────────────────────────────────────────────

STEARATE_SMILES = "CCCCCCCCCCCCCCCCCC(=O)[O-]"  # C18H35O2-
SODIUM_SMILES = "[Na+]"

STEARATE_NUM_ATOMS = 55   # 18C + 2O + 35H
SODIUM_NUM_ATOMS = 1
WATER_NUM_ATOMS = 3

_RDKIT_SEED = 42  # deterministic conformer generation


# ── RDKit helpers ─────────────────────────────────────────────────────


def _embed_stearate() -> Chem.Mol:
    """Return an RDKit Mol with 3-D coordinates for the stearate anion."""
    mol = Chem.MolFromSmiles(STEARATE_SMILES)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {STEARATE_SMILES}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = _RDKIT_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("RDKit embedding failed for stearate")
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


# ── OpenFF molecules ──────────────────────────────────────────────────
# Built from the same RDKit objects used for PDB generation so that the
# atom ordering is guaranteed to match.


def create_stearate_off() -> Molecule:
    """OpenFF Molecule with atom ordering matching :func:`generate_stearate_pdb`."""
    return Molecule.from_rdkit(_embed_stearate(), allow_undefined_stereo=True)


def create_sodium_off() -> Molecule:
    return Molecule.from_smiles(SODIUM_SMILES)


# ── PDB writers ───────────────────────────────────────────────────────


def generate_stearate_pdb(path: Path) -> None:
    mol = _embed_stearate()
    Chem.MolToPDBFile(mol, str(path))
    # Normalise residue name for downstream identification
    text = path.read_text()
    path.write_text(text.replace("UNL", "STL"))
    log.info("Wrote stearate monomer -> %s  (%d atoms)", path, mol.GetNumAtoms())


def generate_sodium_pdb(path: Path) -> None:
    path.write_text(
        "ATOM      1  NA  NA  A   1"
        "      0.000   0.000   0.000  1.00  0.00          NA\n"
        "END\n"
    )
    log.info("Wrote sodium ion -> %s", path)


def generate_water_pdb(path: Path) -> None:
    # TIP3P geometry: O-H = 0.9572 A, H-O-H = 104.52 deg
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


def generate_all_monomers(output_dir: Path) -> dict[str, Path]:
    """Write monomer PDBs and return ``{name: path}`` mapping."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "stearate": output_dir / "stearate.pdb",
        "sodium":   output_dir / "sodium.pdb",
        "water":    output_dir / "water.pdb",
    }
    generate_stearate_pdb(paths["stearate"])
    generate_sodium_pdb(paths["sodium"])
    generate_water_pdb(paths["water"])
    return paths
