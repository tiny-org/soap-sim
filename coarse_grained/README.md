# Coarse-grained soap simulation (MARTINI 2 + GROMACS)

Microsecond-scale simulation of soap / water / polyol mixtures using
the MARTINI 2 coarse-grained force field and GROMACS.

This complements the atomistic OpenMM simulation in `../rewrite/` by
trading molecular detail for speed -- enabling observation of slow
processes like crystallization, phase separation, and mesophase ordering.

## When to use this vs the atomistic simulation

| | Atomistic (OpenMM) | Coarse-grained (GROMACS) |
|---|---|---|
| **Resolution** | All atoms | ~4 atoms per bead |
| **Timescale** | 1-50 ns | 1-100 us |
| **Length scale** | 5-10 nm | 15-50 nm |
| **Use for** | Local structure, H-bonds, charge details | Phase behavior, self-assembly, crystallization |

## Installation

```bash
micromamba activate soap   # reuse the existing environment
micromamba install -c conda-forge gromacs
bash setup.sh             # download MARTINI 2.2 force field files
```

## Usage

```bash
python build.py           # generate system.gro + topol.top
bash run.sh               # minimize -> NVT (1 ns) -> NPT (1 us)
```

## Configuration

Edit `config.toml` to change the system composition:

```toml
[[molecules]]
name     = "STR"          # matches [moleculetype] in molecules/stearate.itp
count    = 200
n_beads  = 5              # CG beads per molecule
is_chain = true           # linear chain (soaps) vs single bead

[system]
box = [15.0, 15.0, 15.0]  # nm
```

Edit `mdp/npt.mdp` to change production run length, temperature, or pressure.

## MARTINI bead mapping

### Fatty acid soaps (anionic)

| Molecule | Config name | Beads | Mapping |
|---|---|---|---|
| Laurate (C12) | `LAU` | 4 | Qa(COO⁻) + 3×C1 |
| Myristate (C14) | `MYR` | 4 | Qa(COO⁻) + 3×C1 |
| Palmitate (C16) | `PAL` | 5 | Qa(COO⁻) + 4×C1 |
| Stearate (C18) | `STR` | 5 | Qa(COO⁻) + 4×C1 |
| Oleate (C18:1) | `OLE` | 5 | Qa + C1 + C3(cis kink) + 2×C1 |
| Ricinoleate (C18:1, 12-OH) | `RIC` | 5 | Qa + C1 + C3(cis kink) + P1(OH) + C1 |

### Counterion

| Molecule | Config name | Beads | Mapping |
|---|---|---|---|
| Na⁺ | `NA+` | 1 | Qd |

### Solvents and polyols

| Molecule | Config name | Beads | Mapping |
|---|---|---|---|
| Water (4 H₂O) | `W` | 1 | P4 |
| Ethanol | `ETH` | 1 | P1 |
| Propylene glycol | `PGO` | 1 | P1 |
| Glycerol | `GLY` | 1 | P4 |
| Sorbitol | `SOR` | 2 | P1 + P1 |
| Sucrose | `SUC` | 3 | P1 + P1 + P1 |

Each C1 bead represents ~4 CH₂ groups. C3 represents an unsaturated
segment (cis double bond → 120° kink). Each W bead represents 4 water molecules.

## Output

```
output/
  system.gro     initial CG coordinates
  topol.top      system topology
  em.gro         minimized structure
  nvt.gro        post-equilibration
  npt.xtc        production trajectory
  npt.edr        energy data
  npt.gro        final frame
```

## Visualization

```bash
# Convert trajectory to PDB for PyMOL (fix PBC wrapping)
echo 0 | gmx trjconv -f output/npt.xtc -s output/npt.tpr -o output/traj.pdb -pbc mol

pymol output/traj.pdb
```

In PyMOL:
```
hide everything
show spheres
color green, resn LAU+MYR+PAL+STR
color orange, resn PGO
color purple, resn NA+
color grey80, resn W
set sphere_scale, 0.5
```

## Comparing crystallization with/without polyol

Run two simulations:

1. **With PG**: use `config.toml` as-is (50% propylene glycol)
2. **Without PG**: edit `config.toml`, remove the PGO entry, increase water
   to compensate

Cool both systems from 80°C to 25°C by editing `ref_t` in `mdp/npt.mdp`
and running 1 μs at each temperature. Compare the final structures --
the system with PG should remain disordered while the pure soap/water
system shows lamellar crystal ordering.
