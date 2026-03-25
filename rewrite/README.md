# soap-sim

Molecular dynamics simulation of sodium stearate / water mixtures using OpenMM.

## What it does

Simulates a 50 % by weight sodium stearate (soap) / water system at 80 °C
and 1 bar using the OpenFF Sage 2.2 force field and TIP3P water.

The simulation protocol:

1. **Build** -- generate monomer structures (RDKit), pack into a periodic box (PACKMOL)
2. **Parameterize** -- assign force-field parameters (OpenFF + SystemGenerator)
3. **Simulate** -- energy minimization, NVT equilibration, NPT production (OpenMM)
4. **Analyze** -- summarize temperature, energy, and density from the production log

## Installation

All heavy dependencies (OpenMM, RDKit, OpenFF) are conda-only:

```bash
# 1. Create environment
micromamba create -n soap python=3.12
micromamba activate soap

# 2. Install dependencies
micromamba install -c conda-forge openmm openmmforcefields openff-toolkit rdkit parmed numpy packmol pymol-open-source

# 3. On Apple Silicon, force CPU mode
export OPENMM_DEFAULT_PLATFORM=CPU

# 4. Install this package (editable)
cd rewrite
pip install -e .
```

Verify OpenMM works: `python -m openmm.testInstallation`

## Usage

All commands read from `config.toml` in the current directory:

```bash
# Remove output if major config changes
rm -rf output/

# Run everything (build + parameterize + simulate + analyze)
python -m soap_sim run

# Or run steps individually:
python -m soap_sim build           # -> output/build/packed.pdb
python -m soap_sim parameterize    # -> output/parameterize/system.xml
python -m soap_sim simulate        # -> output/simulate/trajectory.dcd
python -m soap_sim analyze         # print summary statistics
```

Show 3D model of simulation
```bash
pymol output/parameterize/topology.pdb output/simulate/trajectory.dcd
```

Use a different config file:

```bash
python -m soap_sim run -c my_config.toml
```

Enable verbose logging:

```bash
python -m soap_sim run -v
```

## Configuration

Edit `config.toml` to change any parameter:

```toml
[system]
num_stearate    = 50      # ion pairs
weight_fraction = 0.50    # sodium stearate mass fraction
# water count is computed automatically (851 for 50 NaSt at 50%)

[simulation]
temperature_celsius = 80.0
pressure_bar        = 1.0

[simulation.production]
steps        = 500_000    # 1 ns
dcd_interval = 1000
```

## Output structure

```
output/
├── build/
│   ├── stearate.pdb, sodium.pdb, water.pdb   (monomers)
│   ├── packmol.inp
│   └── packed.pdb                             (packed system)
├── parameterize/
│   ├── system.xml                             (OpenMM System)
│   └── topology.pdb                           (positions + topology)
└── simulate/
    ├── minimized.pdb
    ├── equilibrated.pdb
    ├── trajectory.dcd                         (production trajectory)
    ├── production.csv                         (energy / T / density log)
    ├── checkpoint.chk
    └── final.pdb
```

## Architecture

| Module          | Responsibility                                      |
| --------------- | --------------------------------------------------- |
| `config.py`     | TOML loading, dataclass validation, derived values  |
| `molecules.py`  | RDKit monomer generation, OpenFF molecule creation   |
| `packing.py`    | PACKMOL input generation and execution              |
| `parameterize.py` | OpenFF/SystemGenerator parameterization, serialization |
| `simulate.py`   | Minimization, NVT equilibration, NPT production     |
| `analysis.py`   | Log parsing and summary statistics                  |
| `__main__.py`   | CLI with subcommands                                |

Each module exposes a clean Python API so you can import and customise
individual steps in notebooks or scripts.
