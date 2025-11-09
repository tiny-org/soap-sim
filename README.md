# soap-sim

Soap molecular simulation with OpenMM

## Installation

| Step | Explanation               | Command                                           |
| ---- | ------------------------- | ------------------------------------------------- |
| 1    | Install micromamba        | `"${SHELL}" <(curl -L micro.mamba.pm/install.sh)` |
| 2    | Alias conda -> micromamba | `alias conda micromamba`                          |
| 2    | Install openmm            | `conda install -c conda-forge openmm`             |
| 2    | Activate base env         | `conda activate base`                             |
| 5    | Test installation         | `python -m openmm.testInstallation`               |
