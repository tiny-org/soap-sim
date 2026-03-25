"""PACKMOL integration: input generation, execution, box-size estimation."""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from .config import Config
from .molecules import generate_all_monomers

log = logging.getLogger(__name__)


def _write_packmol_input(config: Config, monomer_dir: Path,
                         output_pdb: Path, inp_path: Path) -> None:
    """Write the PACKMOL control file."""
    sys = config.system
    bx, by, bz = sys.box_angstrom

    inp_path.write_text(
        f"# PACKMOL input -- sodium stearate / water mixture\n"
        f"tolerance 2.0\n"
        f"discale 1.5\n"
        f"movebadrandom\n"
        f"nloop 200\n"
        f"output {output_pdb}\n"
        f"filetype pdb\n"
        f"\n"
        f"# Stearate ions (C18H35O2-)\n"
        f"structure {monomer_dir / 'stearate.pdb'}\n"
        f"  number {sys.num_stearate}\n"
        f"  inside box 0.0 0.0 0.0 {bx} {by} {bz}\n"
        f"end structure\n"
        f"\n"
        f"# Sodium ions (Na+)\n"
        f"structure {monomer_dir / 'sodium.pdb'}\n"
        f"  number {sys.num_sodium}\n"
        f"  inside box 0.0 0.0 0.0 {bx} {by} {bz}\n"
        f"end structure\n"
        f"\n"
        f"# Water (TIP3P)\n"
        f"structure {monomer_dir / 'water.pdb'}\n"
        f"  number {sys.num_water}\n"
        f"  inside box 0.0 0.0 0.0 {bx} {by} {bz}\n"
        f"end structure\n"
    )
    log.info("Wrote PACKMOL input -> %s  (box %.1f x %.1f x %.1f A)",
             inp_path, bx, by, bz)


def _run_packmol(inp_path: Path) -> None:
    """Execute PACKMOL with shell redirection (Fortran needs seekable stdin)."""
    log.info("Running PACKMOL (this may take a minute) ...")
    result = subprocess.run(
        f"packmol < {inp_path}",
        shell=True, text=True, capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"PACKMOL failed (exit {result.returncode}).\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    log.info("PACKMOL finished successfully.")


# ── Public API ────────────────────────────────────────────────────────


def build_system(config: Config) -> Path:
    """Generate monomers, pack them, return path to the packed PDB.

    Directory layout under *config.output.directory*::

        build/
            stearate.pdb
            sodium.pdb
            water.pdb
            packmol.inp
            packed.pdb
    """
    build_dir = Path(config.output.directory) / "build"
    build_dir.mkdir(parents=True, exist_ok=True)

    monomer_dir = build_dir
    generate_all_monomers(monomer_dir)

    packed_pdb = build_dir / "packed.pdb"
    inp_path = build_dir / "packmol.inp"

    _write_packmol_input(config, monomer_dir, packed_pdb, inp_path)
    _run_packmol(inp_path)

    n_total = (config.system.num_stearate * 55
               + config.system.num_sodium
               + config.system.num_water * 3)
    log.info("Packed system: %d stearate, %d Na+, %d water  (%d atoms)",
             config.system.num_stearate, config.system.num_sodium,
             config.system.num_water, n_total)
    return packed_pdb
