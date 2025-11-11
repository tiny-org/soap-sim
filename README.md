# soap-sim

Soap molecular simulation with OpenMM

## Installation

| Step | Explanation               | Command                                                         |
| ---- | ------------------------- | --------------------------------------------------------------- |
| 1    | Install micromamba        | `"${SHELL}" <(curl -L micro.mamba.pm/install.sh)`               |
| 2    | Alias conda -> micromamba | `alias conda micromamba`                                        |
| 3    | Create python env         | `conda create -n soap_sim_env python=3.12`                      |
| 4    | Activate env              | `conda activate soap_sim_env`                                   |
| 3    | Install openmm            | `conda install -c conda-forge openmm openff-toolkit-base rdkit` |
| 5    | Test installation         | `python -m openmm.testInstallation`                             |
| 6    | Set CPU-Mode (on M1)      | `export OPENMM_DEFAULT_PLATFORM=CPU`                            |
| 7    | Go to example dir         | `cd example`                                                    |
| 8    | Run example               | `python simulatePdb.py`                                         |
| 8    | Exit python env           | `conda deactivate`                                              |

## Links

Smiles to PDB Converter: https://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html
