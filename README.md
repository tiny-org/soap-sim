# soap-sim

Soap molecular simulation with OpenMM

## Installation

| Step | Explanation               | Command                                              |
| ---- | ------------------------- | ---------------------------------------------------- |
| 1    | Install micromamba        | `"${SHELL}" <(curl -L micro.mamba.pm/install.sh)`    |
| 2    | Alias conda -> micromamba | `alias conda='micromamba'`                           |
| 3    | Create python env         | `conda create -n soap_sim_env python=3.12`           |
| 4    | Activate env              | `conda activate soap_sim_env`                        |
| 3    | Install openmm            | `conda install -c conda-forge openmm openmmforcefields openff-toolkit packmol` |
| 5    | Test installation         | `python -m openmm.testInstallation`                  |
| 6    | Set CPU-Mode (on M1)      | `export OPENMM_DEFAULT_PLATFORM=CPU`                 |
| 7    | Go to example dir         | `cd example`                                         |
| 8    | Run example               | `python simulatePdb.py`                              |
| 9    | Exit python env           | `conda deactivate`                                   |

* cd src/system
* python builder.py
* packmol < packmol_input.inp

------

The simulation runs in three stages: environment setup, system building, and running the simulation. Here are the steps:

1. Environment Setup

# Install micromamba (if not already installed)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
alias conda='micromamba'

# Create and activate the environment
conda create -n soap_sim_env python=3.12
conda activate soap_sim_env

# Install all dependencies
conda install -c conda-forge openmm openmmforcefields openff-toolkit rdkit parmed numpy packmol

# On Apple Silicon, force CPU mode
export OPENMM_DEFAULT_PLATFORM=CPU

# Verify OpenMM works
python -m openmm.testInstallation

2. Build the System

This generates monomer PDB files, writes the PACKMOL input, and runs PACKMOL to pack 50 stearate ions + 50 sodium ions + 850 water molecules into a 4 nm box.

cd src/system
python builder.py

Output: system_coordinates.pdb

3. Parameterize the System

This assigns OpenFF (Sage 2.2.0) force field parameters to the stearate/sodium, combines with TIP3P water, and sets up PME electrostatics.

# Still in src/system/
python forcefields.py

Output: openmm_system.pdb (the parameterized coordinates ready for simulation)

4. Run the Simulation

A quick 500-step validation run with energy minimization + Langevin dynamics at 300 K:

# Still in src/system/
python run_quick_sim.py

This will print energy and temperature every 100 steps so you can verify the system is stable.

## Links

Smiles to PDB Converter: https://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html


* Generate the force field parameters for the stearate ion
* Use PACKMOL or a similar tool to arrange the stearate ions
* Use ParmEd to combine the PDB/coordinates with the force fields (Stearate + Water) and generate the final Amber .prmtop file