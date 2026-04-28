# AGENTS.md

Agent-focused notes for working in this repository. For human-oriented intro,
see `README.md`.

## Project overview

`soap-sim` is a molecular-dynamics playground for soap (fatty-acid salt) /
water / polyol mixtures. The repository contains three tracks of code, plus a
toy smoke test:

| Path              | Status   | Stack                          | Purpose |
|-------------------|----------|--------------------------------|---------|
| `rewrite/`        | **Active** | OpenMM + OpenFF + PACKMOL    | Atomistic MD of arbitrary surfactant mixtures, TOML-driven CLI |
| `coarse_grained/` | **Active** | GROMACS + MARTINI 2.2        | Coarse-grained MD for microsecond / mesophase behaviour |
| `src/`            | Legacy   | OpenMM + OpenFF + PACKMOL      | First prototype hardcoded to sodium stearate. Kept for reference; do not extend. |
| `example/`        | Demo     | OpenMM only                    | 14-line standalone OpenMM smoke test on a prebuilt PDB |

When asked to "add a feature" or "fix a bug", default to working in
`rewrite/` (atomistic) or `coarse_grained/` (CG) -- whichever fits the task.
Do **not** modify `src/` unless the user explicitly asks; treat it as
historical.

## Environment setup

### Heavy dependencies are conda-only

`openmm`, `openmmforcefields`, `openff-toolkit`, `rdkit`, `packmol`, `gromacs`,
and `pymol-open-source` are **not on PyPI** for our target platform. Always
install them via micromamba/conda. `pip install` of these will fail or pull in
a broken stack.

```bash
# One-time: install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
alias conda='micromamba'   # the rest of this doc uses `conda`

# Create env shared by atomistic and CG tracks
conda create -n soap python=3.12
conda activate soap
conda install -c conda-forge \
    openmm openmmforcefields openff-toolkit rdkit numpy packmol \
    pymol-open-source

# Atomistic track: install the local package in editable mode
cd rewrite && pip install -e . && cd ..

# CG track: also install GROMACS and download MARTINI 2.2
conda install -c conda-forge gromacs
bash coarse_grained/setup.sh

# Verify OpenMM
python -m openmm.testInstallation
```

### Apple Silicon (M-series)

OpenMM's GPU backends are unreliable on macOS ARM. Force CPU before running
any OpenMM script:

```bash
export OPENMM_DEFAULT_PLATFORM=CPU
```

The `simulate.py` platform selector falls back to CPU automatically when
`platform = "auto"` in the config, but the env var is the surest fix.

## Atomistic track (`rewrite/`)

### Layout

```
rewrite/
  pyproject.toml                  # editable install; only `numpy` listed
  config_soap.toml                # 4-soap reference mix
  config_soap_pg.toml             # 4 soaps + propylene glycol
  config_lamellar.toml            # 50% NaStearate, anisotropic barostat
  src/soap_sim/
    __main__.py                   # CLI: build / parameterize / simulate / analyze / run
    config.py                     # TOML loader + frozen dataclasses
    molecules.py                  # RDKit -> OpenFF Molecule generation from SMILES
    packing.py                    # PACKMOL input + execution
    parameterize.py               # OpenFF SystemGenerator -> serialised system.xml
    simulate.py                   # minimize -> NVT -> NPT
    analysis.py                   # parse production.csv, summarise
```

### Running

The CLI is `python -m soap_sim`. Always pass `-c <config.toml>`; the default
`config.toml` does not exist.

```bash
cd rewrite

# Full pipeline: build -> parameterize -> simulate -> analyze
python -m soap_sim -c config_soap.toml run

# Individual stages (each reads outputs of the previous one)
python -m soap_sim -c config_soap.toml build
python -m soap_sim -c config_soap.toml parameterize
python -m soap_sim -c config_soap.toml simulate
python -m soap_sim -c config_soap.toml analyze

# Verbose logging
python -m soap_sim -v -c config_soap.toml run
```

The console-script alias `soap-sim` (declared in `pyproject.toml`) is also
installed, but `python -m soap_sim` is preferred in scripts because it is
unambiguous about which interpreter runs.

### Output layout

Stages communicate by writing into `<output.directory>` from the config:

```
<output_dir>/
  build/         monomers/*.pdb, packmol.inp, packed.pdb
  parameterize/  system.xml, topology.pdb
  simulate/      minimized.pdb, equilibrated.pdb, equilibration.csv,
                 trajectory.dcd, production.csv, checkpoint.chk, final.pdb
```

The top-level `output*/` directories are gitignored. Never commit them.

### Config conventions

- TOML; loaded by `soap_sim.config.load_config`.
- Solutes are a list of SMILES + count; counterions are auto-computed from
  net charge. Multiple counterion species supported via fractions that must
  sum to `1.0`.
- Either `num_water` (explicit) **or** `water_weight_fraction` (auto-compute)
  must be set. If both appear, the loader rejects the unused key as an
  unknown setting -- pick one.
- `box_dimensions` overrides the auto-sized cubic box (use Angstroms; CG track
  uses nm -- do not mix).
- `barostat = "anisotropic"` is required for lamellar / mesophase systems.
  Default `isotropic` is wrong for layered structures.
- `temperature_celsius`, not Kelvin (converted internally).
- All dataclasses are `frozen=True`; never mutate config at runtime.

### Adding a new molecule

The atomistic track is fully SMILES-driven -- there is no separate parameter
file to maintain. To add e.g. coconut soap, just add another
`[[system.solutes]]` block to a TOML. The first 3 chars of `name` become the
PDB residue name, or set `residue_name = "FOO"` explicitly. Keep residue
names unique within a config; PyMOL colouring keys off them.

## Coarse-grained track (`coarse_grained/`)

### Layout

```
coarse_grained/
  config.toml                     # composition + box
  build.py                        # writes output/system.gro + output/topol.top
  setup.sh                        # downloads MARTINI 2.2 + martini_v2.0_ions ITPs
  run.sh                          # gmx grompp/mdrun for em -> nvt -> npt
  mdp/{em,nvt,npt}.mdp            # GROMACS run-control
  molecules/*.itp                 # MARTINI topology fragments per species
  forcefield/                     # downloaded by setup.sh, gitignored
  output/                         # all artefacts, gitignored
```

### Running

```bash
cd coarse_grained
bash setup.sh                     # one-time: GROMACS check + MARTINI download
python build.py                   # generates output/system.gro + topol.top
bash run.sh                       # em -> nvt (1 ns) -> npt (1 us)
```

`build.py` does *not* call PACKMOL; it places chains as random straight
segments inside the box, then GROMACS minimisation cleans up the geometry.
Do not import the atomistic `soap_sim` package here -- the two tracks are
independent.

### Adding a new molecule (CG)

1. Add an `[moleculetype]` ITP in `coarse_grained/molecules/<name>.itp`.
2. Register it in `build.py`:
   - bead names in `BEAD_NAMES`
   - ITP file mapping in `ITP_MAP`
   - charge in `CHARGE` (only if it is a counterion or charged solute)
3. Reference it from `config.toml` under `[[solutes]]` / `[[solvents]]` /
   `[[counterions]]`. The `name` must match the `[moleculetype]` molname in
   the ITP (e.g. `STR`, not `stearate`).

## Code style

The codebase is small and intentionally idiomatic Python. Match the existing
patterns rather than introducing new frameworks.

- Python 3.12+. `from __future__ import annotations` at the top of every
  module under `rewrite/src/soap_sim/`.
- Frozen `@dataclass` for config; `field(default_factory=...)` for nested
  defaults.
- TOML via stdlib `tomllib`; do not add `tomli`/`pydantic`/`click`.
- CLI is `argparse` only -- no Click, Typer, etc.
- `logging` for diagnostics (`log = logging.getLogger(__name__)`); print only
  at top-level CLI boundaries (`__main__.py`).
- `@lru_cache` for repeated RDKit calls keyed on SMILES.
- Type hints on public functions; PEP 8 spacing; single-line docstrings on
  helpers, multi-line on public APIs.
- Section comments `# -- Section ---` separate logical groups within a module.
- No test framework set up; do not introduce one without asking.

## Verification (no formal test suite)

There are no unit tests, lint configs, or CI. To validate a change:

1. Smoke-test the atomistic CLI end-to-end on a tiny config:
   ```bash
   cd rewrite
   # Create a throwaway config with small counts and few steps.
   # (See "Quick smoke config" below.)
   python -m soap_sim -v -c smoke.toml run
   ```
   A successful run prints `Production run summary` with non-NaN temperature
   / density / energies. Anything else (PACKMOL nonzero exit, "PME was not
   applied", net charge != 0) is a bug.

2. Smoke-test the CG track:
   ```bash
   cd coarse_grained
   python build.py        # must end with "System: N CG beads in ... nm box"
   bash run.sh            # em + nvt + npt should finish; npt is slow, ok to
                          # truncate `nsteps` in mdp/npt.mdp for a quick check
   ```

3. Check that PACKMOL/OpenFF parameterisation succeed with `-v`. Verbose logs
   print the net charge after parameterisation -- it must round to `0.0` for
   neutralised systems.

### Quick smoke config

When iterating on rewrite/ code, copy this into `rewrite/smoke.toml` (it is
gitignored under `output*/` rules; do not commit):

```toml
[[system.solutes]]
name = "stearate"
smiles = "CCCCCCCCCCCCCCCCCC(=O)[O-]"
count = 10

[system]
num_water = 200
target_density = 0.9

[simulation]
temperature_celsius = 80.0
[simulation.minimize]
max_iterations = 100
[simulation.equilibrate]
steps = 500
[simulation.production]
steps = 1000
dcd_interval = 100
log_interval = 100
checkpoint_interval = 500

[output]
directory = "output_smoke"
```

This finishes in well under a minute on CPU and exercises every code path.

## Visualization

PyMOL is the standard viewer; install via conda (`pymol-open-source`).

```bash
# Atomistic
pymol rewrite/output_soap/parameterize/topology.pdb \
      rewrite/output_soap/simulate/trajectory.dcd

# CG: convert XTC -> PDB first to fix periodic-boundary wrapping
echo 0 | gmx trjconv -f coarse_grained/output/npt.xtc \
                     -s coarse_grained/output/npt.tpr \
                     -o coarse_grained/output/traj.pdb -pbc mol
pymol coarse_grained/output/traj.pdb
```

Standard colour scheme (paste into PyMOL command line; adjust resnames to
match the active config):

```
hide everything
show sticks, resn LAU+MYR+PAL+STE
show sticks, resn PGO
show spheres, resn NA+
show lines, resn HOH
color forest, resn STE
color chartreuse, resn PAL
color limon, resn MYR
color tv_green, resn LAU
color orange, resn PGO
color purple, resn NA+
color grey80, resn HOH
set sphere_scale, 0.4
```

## Common gotchas

- **PACKMOL hangs / fails silently**: the Fortran runtime needs a seekable
  stdin; `packing.py` shells out via `packmol < input.inp`. Never call
  `subprocess.run([...], stdin=PIPE)` here -- it will deadlock.
- **"System is not periodic -- PME was not applied"**: box vectors must be
  set on `modeller.topology` *before* `SystemGenerator.create_system`. See
  `parameterize.py` step 5.
- **Lamellar / membrane systems collapse**: switch to
  `barostat = "anisotropic"` so the z-axis can relax independently.
- **Net charge not zero after parameterisation**: counterion auto-count
  divides by `abs(counterion charge)`. Mixed +1/+2 counterions are not
  supported; use only +1 ions (Na+, K+).
- **CG bead/atom name mismatch**: residue names in MARTINI ITPs must match
  the `name` in `coarse_grained/config.toml` and the keys in `BEAD_NAMES` /
  `ITP_MAP` in `build.py`. Keep all three in sync.
- **Apple Silicon "Operation not permitted"**: see Apple Silicon section --
  set `OPENMM_DEFAULT_PLATFORM=CPU`.

## Pull requests / commits

- Commit messages are short, lowercase, imperative ("add potassium counter
  ions for coarse sim", "fix counterions"). Match this style; no Conventional
  Commits prefixes.
- Do not commit anything under `output/`, `output_*/`, `movies/`,
  `coarse_grained/forcefield/`, or `*.dcd` / `*.xtc` / `*.edr` -- they are
  gitignored for a reason (large binary trajectories).
- Do not commit `*.egg-info/` (it appears under `rewrite/src/` after
  `pip install -e .`; gitignored).
- No CI is wired up. There is nothing to "make green"; manual smoke tests are
  the bar.

## What not to do

- Do not refactor `src/` (the legacy prototype). Improvements belong in
  `rewrite/`.
- Do not add `pip` dependencies for `openmm`, `rdkit`, `openff-toolkit`,
  `packmol`, `gromacs`, or `pymol`. They are conda-only.
- Do not introduce new top-level packages or rename the `soap_sim`
  module without first checking the user; the CLI entry point and console
  script depend on it.
- Do not switch atomistic and CG units silently -- atomistic uses Angstroms
  in TOML and nm in OpenMM; CG uses nm everywhere.
