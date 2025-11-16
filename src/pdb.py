"""Generate a sodium stearate structure and export it as a PDB file."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
from openmm import unit
from openmm.app import PDBFile
from openff.toolkit.topology import Molecule, Topology
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

STEARATE_SMILES = "CCCCCCCCCCCCCCCCC(=O)[O-]"
SODIUM_SMILES = "[Na+]"
DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "c18-na.pdb"


def _embed_with_rdkit(smiles: str, add_hs: bool = True) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES string: {smiles}")

    if add_hs:
        mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 2024
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise RuntimeError(f"RDKit embedding failed for {smiles} (status {status})")

    # MMFF optimization provides a reasonable initial geometry
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _rdkit_positions(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array(conf.GetPositions(), dtype=np.float64)


def _carboxylate_indices(stearate: Chem.Mol) -> Tuple[int, Tuple[int, int]]:
    carbon_index = None
    oxygen_indices: Tuple[int, int] | None = None

    for atom in stearate.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        oxygen_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"]
        if len(oxygen_neighbors) == 2:
            carbon_index = atom.GetIdx()
            oxygen_indices = (oxygen_neighbors[0], oxygen_neighbors[1])
            break

    if carbon_index is None or oxygen_indices is None:
        raise ValueError("Unable to locate the carboxylate carbon and oxygens in stearate")

    return carbon_index, oxygen_indices


def _place_sodium_near_headgroup(stearate: Chem.Mol) -> np.ndarray:
    positions = _rdkit_positions(stearate)
    carbon_idx, oxygen_indices = _carboxylate_indices(stearate)

    carbon_coord = positions[carbon_idx]
    oxygen_coords = positions[list(oxygen_indices)]
    oxygen_centroid = oxygen_coords.mean(axis=0)

    direction = oxygen_centroid - carbon_coord
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        direction = np.array([1.0, 0.0, 0.0])
    else:
        direction /= norm

    return oxygen_centroid + direction * 2.5  # 2.5 Å separation


def build_sodium_stearate() -> Tuple[Topology, unit.Quantity]:
    stearate_rdkit = _embed_with_rdkit(STEARATE_SMILES)
    sodium_coord = _place_sodium_near_headgroup(stearate_rdkit)

    sodium_rdkit = Chem.MolFromSmiles(SODIUM_SMILES)
    if sodium_rdkit is None:
        raise ValueError("Failed to construct sodium molecule from SMILES")

    sodium_conformer = Chem.Conformer(1)
    sodium_conformer.SetAtomPosition(0, Point3D(*sodium_coord))
    sodium_rdkit.AddConformer(sodium_conformer, assignId=True)

    stearate_positions = _rdkit_positions(stearate_rdkit)
    sodium_positions = np.array([sodium_coord])
    combined_positions = np.vstack([stearate_positions, sodium_positions]) * unit.angstroms

    off_molecules = [
        Molecule.from_rdkit(stearate_rdkit, allow_undefined_stereo=True),
        Molecule.from_rdkit(sodium_rdkit, allow_undefined_stereo=True),
    ]
    off_topology = Topology.from_molecules(off_molecules)

    return off_topology, combined_positions


def write_pdb(topology: Topology, positions: unit.Quantity, output_path: Path) -> None:
    openmm_topology = topology.to_openmm()

    # Rename residues for readability (OpenMM allows mutation here)
    residues = list(openmm_topology.residues())
    if residues:
        residues[0].name = "STR"
    if len(residues) > 1:
        residues[1].name = "NA"

    with output_path.open("w") as handle:
        PDBFile.writeFile(openmm_topology, positions, handle)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a sodium stearate PDB structure")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Path to the PDB file to be written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    topology, positions = build_sodium_stearate()
    write_pdb(topology, positions, args.output)
    print(f"Wrote sodium stearate structure to {args.output}")


if __name__ == "__main__":
    main()
