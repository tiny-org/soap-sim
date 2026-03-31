"""PACKMOL integration: input generation, execution, box-size estimation."""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from .config import Config
from .molecules import generate_all_monomers, atom_count

log = logging.getLogger(__name__)


def _write_packmol_input(config: Config, monomer_dir: Path,
                         output_pdb: Path, inp_path: Path) -> None:
    """Write the PACKMOL control file."""
    sys = config.system
    bx, by, bz = sys.box_angstrom
    box_str = f"{bx} {by} {bz}"

    lines = [
        "# PACKMOL input -- soap / water mixture",
        "tolerance 2.0",
        "discale 1.5",
        "movebadrandom",
        "nloop 200",
        f"output {output_pdb}",
        "filetype pdb",
        "",
    ]

    # One block per solute type
    for spec in sys.solutes:
        lines += [
            f"# {spec.name} ({spec.smiles})",
            f"structure {monomer_dir / f'{spec.name}.pdb'}",
            f"  number {spec.count}",
            f"  inside box 0.0 0.0 0.0 {box_str}",
            "end structure",
            "",
        ]

    # One block per counterion type
    for ci, count in sys.counterion_counts:
        lines += [
            f"# Counterion ({ci.smiles}, {ci.resname})",
            f"structure {monomer_dir / f'counterion_{ci.resname.lower()}.pdb'}",
            f"  number {count}",
            f"  inside box 0.0 0.0 0.0 {box_str}",
            "end structure",
            "",
        ]

    # Water
    lines += [
        "# Water (TIP3P)",
        f"structure {monomer_dir / 'water.pdb'}",
        f"  number {sys.num_water}",
        f"  inside box 0.0 0.0 0.0 {box_str}",
        "end structure",
    ]

    inp_path.write_text("\n".join(lines) + "\n")
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
    """Generate monomers, pack them, return path to the packed PDB."""
    build_dir = Path(config.output.directory) / "build"
    build_dir.mkdir(parents=True, exist_ok=True)

    generate_all_monomers(config, build_dir)

    packed_pdb = build_dir / "packed.pdb"
    inp_path = build_dir / "packmol.inp"

    _write_packmol_input(config, build_dir, packed_pdb, inp_path)
    _run_packmol(inp_path)

    sys = config.system
    n_total = (
        sum(s.count * atom_count(s.smiles) for s in sys.solutes)
        + sum(n * atom_count(ci.smiles) for ci, n in sys.counterion_counts)
        + sys.num_water * 3
    )
    for s in sys.solutes:
        log.info("  %s: %d molecules", s.name, s.count)
    for ci, n in sys.counterion_counts:
        log.info("  %s: %d ions", ci.resname, n)
    log.info("  water: %d  |  total atoms: %d", sys.num_water, n_total)
    return packed_pdb
