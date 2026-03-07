import os
import subprocess
import numpy as np
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# --- Simulation Parameters ---
# Reduced system size for "small" simulation
NUM_STEARATE = 50
NUM_SODIUM = 50
NUM_WATER = 850
# Estimated box size for this smaller system:
# Volume approx: (50 * ~300 A^3) + (850 * ~30 A^3) = 15000 + 25500 = 40500 A^3
# Cube root(40500) ~= 34.3 A. We use 40.0 A (4.0 nm) to be safe and allow relaxation.
BOX_SIZE = 40.0 # Angstroms (4.0 nm)

# --- Filenames ---
STEARATE_PDB = 'stearate_monomer.pdb'
SODIUM_PDB = 'sodium_ion.pdb'
WATER_PDB = 'water_monomer.pdb'
PACKMOL_INPUT = 'packmol_input.inp'
OUTPUT_PDB = 'system_coordinates.pdb'

# --- 1. Functions to Generate Monomer PDBs ---

def generate_stearate_pdb():
    """Generates a 3D structure for the Stearate ion using OpenFF."""
    # SMILES for the Stearate ion (C18H35O2-) - FIXED SMILES (18 Carbons total)
    stearate_smiles = "CCCCCCCCCCCCCCCCCC(=O)[O-]"
    print(f"Generating 3D coordinates for {stearate_smiles}...")

    # Create RDKit Molecule object and generate an initial conformer
    try:
        mol = Chem.MolFromSmiles(stearate_smiles)
        mol = Chem.AddHs(mol)
        # Generate a reasonable 3D conformer
        AllChem.EmbedMolecule(mol)
        
        # Write to PDB format
        Chem.MolToPDBFile(mol, STEARATE_PDB)

        # Post-process the PDB to ensure a clear residue name (STL) that PACKMOL
        # and the downstream parameterization script can reliably detect.
        try:
            with open(STEARATE_PDB, 'r') as f:
                lines = f.readlines()

            fixed = []
            for line in lines:
                # Replace common unnamed residue labels (UNL) with STL.
                # This is a simple textual substitution that targets the residue name
                # field that many OpenFF PDB writers use; it's sufficient for our
                # generated monomer PDBs.
                fixed.append(line.replace('UNL', 'STL'))

            with open(STEARATE_PDB, 'w') as f:
                f.writelines(fixed)
        except Exception:
            # If post-processing fails, continue but warn the user.
            print('Warning: could not normalize residue name in stearate PDB.')

        print(f"Successfully generated and saved {STEARATE_PDB}")
        return True
    except Exception as e:
        print(f"Error generating stearate PDB: {e}")
        print("Please ensure rdkit is installed (conda install -c conda-forge rdkit)")
        return False

def generate_simple_pdbs():
    """Writes simple PDB files for the Sodium ion (Na+) and Water (TIP3P-like)."""
    
    # --- Sodium Ion (Na+) PDB ---
    # Atom name 'NA' and residue name 'NA' so downstream code can detect Na residues
    # reliably. Use a standard PDB ATOM record.
    na_pdb_content = (
        "ATOM      1  NA  NA  A   1      0.000   0.000   0.000  1.00  0.00          NA\n"
        "END\n"
    )
    with open(SODIUM_PDB, 'w') as f:
        f.write(na_pdb_content)
    print(f"Generated {SODIUM_PDB}")

    # --- Water (H2O) PDB (TIP3P-like geometry for compatibility) ---
    # Use residue name 'HOH' and atom names O, H1, H2 to match common OpenMM/Gromacs
    # conventions and the downstream parameterizer.
    water_pdb_content = (
        "ATOM      1  O   HOH A   1      0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      2  H1  HOH A   1      0.000   0.800   0.600  1.00  0.00           H\n"
        "ATOM      3  H2  HOH A   1      0.000  -0.800   0.600  1.00  0.00           H\n"
        "END\n"
    )
    with open(WATER_PDB, 'w') as f:
        f.write(water_pdb_content)
    print(f"Generated {WATER_PDB}")

# --- 2. Function to Write PACKMOL Input File ---

def write_packmol_input():
    """Writes the PACKMOL input control file (.inp)."""
    box_min = 0.0
    box_max = BOX_SIZE
    
    # PACKMOL coordinates are usually in Angstroms if using PDB/XYZ inputs
    
    packmol_input = f"""
# PACKMOL Input File: Sodium Stearate Solution
# Total Atoms: ~111,000
# Box Size: {BOX_SIZE / 10:.2f} nm per side

# Global tolerance for atom-atom separation (in Angstroms)
tolerance 2.0
output {OUTPUT_PDB}
filetype pdb

# The box limits for the 'inside box' command are:
# X_min={box_min}, Y_min={box_min}, Z_min={box_min}
# X_max={box_max}, Y_max={box_max}, Z_max={box_max}

# --- Stearate Ions (C18H35O2-) ---
structure {STEARATE_PDB}
  number {NUM_STEARATE}
  # Syntax: inside box x_min y_min z_min x_max y_max z_max
  inside box {box_min} {box_min} {box_min} {box_max} {box_max} {box_max}
end structure

# --- Sodium Ions (Na+) ---
structure {SODIUM_PDB}
  number {NUM_SODIUM}
  inside box {box_min} {box_min} {box_min} {box_max} {box_max} {box_max}
end structure

# --- Water Molecules (H2O) ---
structure {WATER_PDB}
  number {NUM_WATER}
  inside box {box_min} {box_min} {box_min} {box_max} {box_max} {box_max}
end structure
"""
    with open(PACKMOL_INPUT, 'w') as f:
        f.write(packmol_input.strip())
    print(f"\nGenerated PACKMOL input file: {PACKMOL_INPUT}")
    print(f"Target system size: {NUM_STEARATE} stearates, {NUM_SODIUM} sodium ions, {NUM_WATER} water molecules.")

# --- 3. Function to Run PACKMOL ---

def run_packmol():
    """Executes the PACKMOL program."""
    print("\n--- Running PACKMOL ---")
    print(f"This may take several minutes to place {NUM_STEARATE + NUM_SODIUM + NUM_WATER} molecules.")
    try:
        # Assuming 'packmol' executable is in your system's PATH
        result = subprocess.run(['packmol'], input=PACKMOL_INPUT, text=True, capture_output=True, check=True)
        print("\nPACKMOL completed successfully!")
        print(f"Output saved to: {OUTPUT_PDB}")
        
        # Optional: Print a snippet of the PACKMOL output
        # print("\nPACKMOL Output Snippet:")
        # print(result.stdout[-1000:]) 
        
    except subprocess.CalledProcessError as e:
        print("\n!!! PACKMOL ERROR !!!")
        print("Execution failed. Ensure 'packmol' is installed and in your system PATH.")
        print(f"Error details: {e}")
        print(f"PACKMOL Output: {e.stdout}")
        print(f"PACKMOL Error: {e.stderr}")
    except FileNotFoundError:
        print("\n!!! PACKMOL EXECUTABLE NOT FOUND !!!")
        print("Please ensure 'packmol' is installed and available in your system's PATH.")


# --- Main Execution Block ---
if __name__ == "__main__":
    
    # 1. Generate all required monomer PDB files
    if not generate_stearate_pdb():
        print("Cannot proceed without the stearate monomer structure.")
    else:
        generate_simple_pdbs()
        
        # 2. Write the PACKMOL control file
        write_packmol_input()
        
        # 3. Run PACKMOL to generate the final system PDB
        # This step requires the 'packmol' executable to be installed.
        # run_packmol()
        
        print("\n--- Next Step ---")
        print(f"A new stearate monomer PDB has been generated.")
        print(f"You MUST manually execute PACKMOL again to generate a new system_coordinates.pdb with the C18 molecule:")
        print(f"    packmol < {PACKMOL_INPUT}")
        print(f"Once done, the {OUTPUT_PDB} file will contain the coordinates for the entire system.")
        print(f"You can then proceed to the parameterization script using the {OUTPUT_PDB} file.")