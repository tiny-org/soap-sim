import numpy as np
from openmm import unit
from openmm.app import PDBFile, Modeller, ForceField

# OpenFF imports
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField as OpenFFForceField

# RDKit imports for structure generation
from rdkit import Chem
from rdkit.Chem import AllChem

# --- 1. Define Molecules and SMILES ---
# Stearate anion (C18H35O2-, i.e. CH3(CH2)16COO-) and Sodium cation (Na+)
STEARATE_SMILES = 'CCCCCCCCCCCCCCCCCC(=O)[O-]'
SODIUM_SMILES = '[Na+]'
ION_COORDS = [] # List to hold RDKit-generated coordinates

# --- 2. Generate RDKit Molecules with 3D Coordinates ---
# The goal is to create one clean set of coordinates for all ions.

# Stearate Generation
stearate_mol = Chem.MolFromSmiles(STEARATE_SMILES)
stearate_mol = Chem.AddHs(stearate_mol) # Add explicit hydrogens (CRITICAL for MD)
AllChem.EmbedMolecule(stearate_mol, AllChem.ETKDGv2())
AllChem.MMFFOptimizeMolecule(stearate_mol) # Perform simple optimization
ION_COORDS.append(stearate_mol.GetConformer().GetPositions())

# Sodium Ion Generation (single atom)
sodium_mol = Chem.MolFromSmiles(SODIUM_SMILES)

# FIX: Explicitly add a single conformer to the sodium molecule object
sodium_mol.AddConformer(Chem.Conformer(1)) 

# Position the sodium ion near the carboxylate head
carboxylate_C_coord = ION_COORDS[0][-1] 
sodium_coord = carboxylate_C_coord + np.array([0.3, 0.0, 0.0]) 

# Set the position of the first (and only) atom in the conformer
# The index 0 refers to the first conformer, which is what SetAtomPosition uses by default.
sodium_mol.GetConformer().SetAtomPosition(0, sodium_coord)
ION_COORDS.append(np.array([sodium_coord]))

# Combine coordinates into a single array
all_coords = np.concatenate(ION_COORDS)
# Convert coordinates to OpenMM compatible units (nanometers)
positions = all_coords * unit.angstroms

# --- 3. Create OpenMM Topology and OpenFF Topology ---
# OpenFF needs the Molecule objects (for chemical identity) and the OpenMM Topology (for atom connectivity).

# Create the combined OpenMM Topology manually (since RDKit PDB output can be messy)
# We use the OpenFF toolkit's internal capabilities for this.
off_mols = [Molecule(stearate_mol), Molecule(sodium_mol)]
off_topology = Topology.from_molecules(off_mols)

# Fixed code:
# --- 4. Parameterization with OpenFF (Sage 2.1.0) ---
print("Applying Open Force Field (Sage 2.1.0) parameters...")

# Load the desired SMIRNOFF force field
openff_forcefield = OpenFFForceField('openff-2.1.0.offxml')

# Apply parameters to create the OpenMM System without extra arguments
system = openff_forcefield.create_openmm_system(off_topology)

# --- 5. Solvation and Final Setup (OpenMM Modeller) ---
print("Creating solvated system...")

# Initialize Modeller using the OpenFF Topology and the RDKit coordinates
modeller = Modeller(off_topology.to_openmm(), positions)

# Add Solvent (Water and Ions)
# Note: For water, we load the standard TIP3P parameters often used alongside OpenFF.
# You might need to use a separate XML if the default Sage doesn't cover TIP3P entirely,
# but often it's included or handled by OpenFF's internal logic.
# For simplicity, we assume TIP3P is covered by a standard XML available to OpenMM.
# We'll load a standard OpenMM force field for the solvent.
openmm_forcefield = ForceField('amber14-all.xml', 'tip3p.xml')

modeller.addSolvent(
    openmm_forcefield,
    model='tip3p',
    # 1.0 nm padding to create a cubic box around the solute
    padding=1.0 * unit.nanometer,
    # System is neutral, so we don't need to add extra salt or neutralize
    neutralize=False
)

# --- 6. Final Outputs ---
print(f"\nSolvated system size: {modeller.topology.getNumAtoms()} atoms")
print("System setup complete. Ready for Minimization and Equilibration.")

# Save the final prepared system (optional, but good for inspection)
with open('solvated_sodium_stearate.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

# system, modeller.topology, and modeller.positions are now ready for the OpenMM simulation context.
