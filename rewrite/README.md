# soap-sim

Molecular dynamics simulation of soap / water mixtures using OpenMM.

Supports arbitrary surfactant mixtures defined by SMILES -- configure
any combination of fatty-acid soaps, counterions, and water via a single
TOML file.

## Installation

```bash
micromamba create -n soap python=3.12
micromamba activate soap
micromamba install -c conda-forge \
    openmm openmmforcefields openff-toolkit rdkit parmed numpy packmol \
    pymol-open-source
export OPENMM_DEFAULT_PLATFORM=CPU   # Apple Silicon only
cd rewrite && pip install -e .
python -m openmm.testInstallation
```

## Usage

```bash
# Full pipeline: build -> parameterize -> simulate -> analyze
python -m soap_sim -c config_soap.toml run

# Individual steps
python -m soap_sim -c config_soap.toml build
python -m soap_sim -c config_soap.toml parameterize
python -m soap_sim -c config_soap.toml simulate
python -m soap_sim -c config_soap.toml analyze

# Verbose logging
python -m soap_sim -v -c config_soap.toml run
```

## Configuration

All parameters live in a TOML file. Define solutes by SMILES:

```toml
[[system.solutes]]
name   = "laurate"
smiles = "CCCCCCCCCCCC(=O)[O-]"
count  = 69

[[system.solutes]]
name   = "stearate"
smiles = "CCCCCCCCCCCCCCCCCC(=O)[O-]"
count  = 50

[system]
num_water      = 851          # explicit count
# water_weight_fraction = 0.20  # alternative: auto-compute from mass fraction
target_density = 0.85         # g/cm^3, for initial box size estimate
# box_dimensions = [50, 50, 200]  # override auto box (Angstroms)

[simulation]
temperature_celsius = 80.0
pressure_bar        = 1.0
barostat            = "isotropic"   # or "anisotropic" for lamellar phases

[simulation.production]
steps        = 5_000_000
dcd_interval = 5000

[output]
directory = "output_soap"
```

See `config_soap.toml` (4-soap mixture) and `config_lamellar.toml` (lamellar phase)
for complete examples.

## Visualization

Use PyMOL to inspect structures and play back trajectories:

```bash
# View the final frame
pymol output/simulate/final.pdb

# Play back the production trajectory
pymol output/parameterize/topology.pdb output/simulate/trajectory.dcd
```

Useful PyMOL commands once inside the viewer:

```
hide everything              # clear the view
show sticks, resn LAU+MYR+PAL+STE   # show soap chains as sticks
show spheres, resn NA        # show counterions as spheres
show lines, resn HOH         # show water as thin lines
set sphere_scale, 0.4        # shrink spheres
color green, resn LAU+MYR+PAL+STE
color purple, resn NA
preset.ball_and_stick()      # alternative: ball-and-stick for everything
```

## Output

```
output/
  build/           monomers + packmol.inp + packed.pdb
  parameterize/    system.xml + topology.pdb
  simulate/        minimized.pdb, equilibrated.pdb, trajectory.dcd,
                   production.csv, checkpoint.chk, final.pdb
```

## Architecture

| Module           | Responsibility                                       |
| ---------------- | ---------------------------------------------------- |
| `config.py`      | TOML loading, dataclass validation, derived values   |
| `molecules.py`   | RDKit monomer generation from arbitrary SMILES       |
| `packing.py`     | PACKMOL input generation and execution               |
| `parameterize.py`| OpenFF/SystemGenerator parameterization              |
| `simulate.py`    | Minimize -> NVT equilibration -> NPT production      |
| `analysis.py`    | Log parsing and summary statistics                   |
| `__main__.py`    | CLI with subcommands                                 |
