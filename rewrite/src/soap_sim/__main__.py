"""CLI entry point:  python -m soap_sim <command> [-c config.toml]"""
from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path

from .config import Config, load_config

# Silence FutureWarning from torch via openff-interchange
warnings.filterwarnings("ignore", category=FutureWarning,
                        module=r"openff\.interchange")
# Silence pint redefining units
warnings.filterwarnings("ignore", message=r"Redefining",
                        module=r"pint")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(levelname)-8s %(name)s  %(message)s",
        level=level,
        stream=sys.stdout,
    )
    # Suppress noisy third-party loggers
    if not verbose:
        for name in (
            "openff.interchange",
            "openmmforcefields",
            "pint",
        ):
            logging.getLogger(name).setLevel(logging.WARNING)


def _load(args: argparse.Namespace) -> Config:
    if args.config.exists():
        return load_config(args.config)
    logging.warning("Config %s not found -- using defaults", args.config)
    return Config()


# ── Subcommands ───────────────────────────────────────────────────────


def cmd_build(args: argparse.Namespace) -> None:
    from .packing import build_system

    config = _load(args)
    sys_cfg = config.system
    bx, by, bz = sys_cfg.box_angstrom
    print(f"System: {sys_cfg.num_stearate} NaStearate + {sys_cfg.num_water} "
          f"H2O  ({sys_cfg.weight_fraction:.0%} by weight)")
    print(f"Box: {bx:.1f} x {by:.1f} x {bz:.1f} A "
          f"({bx/10:.2f} x {by/10:.2f} x {bz/10:.2f} nm)")

    packed = build_system(config)
    print(f"\nPacked system written to {packed}")


def cmd_parameterize(args: argparse.Namespace) -> None:
    from .parameterize import parameterize_system, save_system

    config = _load(args)
    packed_pdb = Path(config.output.directory) / "build" / "packed.pdb"
    if not packed_pdb.exists():
        sys.exit(f"ERROR: {packed_pdb} not found.  Run 'build' first.")

    system, topology, positions = parameterize_system(packed_pdb, config)
    save_system(system, topology, positions, Path(config.output.directory))
    print("\nParameterization complete.")


def cmd_simulate(args: argparse.Namespace) -> None:
    from .parameterize import load_system
    from .simulate import run_simulation

    config = _load(args)
    out = Path(config.output.directory)
    if not (out / "parameterize" / "system.xml").exists():
        sys.exit("ERROR: system.xml not found.  Run 'parameterize' first.")

    system, topology, positions = load_system(out)
    run_simulation(system, topology, positions, config)


def cmd_analyze(args: argparse.Namespace) -> None:
    from .analysis import summarize

    config = _load(args)
    summarize(config)


def cmd_run(args: argparse.Namespace) -> None:
    """Run the full pipeline: build -> parameterize -> simulate -> analyze."""
    cmd_build(args)
    cmd_parameterize(args)
    cmd_simulate(args)
    cmd_analyze(args)


# ── Argument parser ───────────────────────────────────────────────────


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="soap-sim",
        description="Molecular dynamics of sodium stearate / water mixtures",
    )
    parser.add_argument(
        "-c", "--config", type=Path, default=Path("config.toml"),
        help="Path to TOML configuration file (default: config.toml)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable debug-level logging",
    )

    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("build",        help="Generate monomers and pack with PACKMOL")
    sub.add_parser("parameterize", help="Assign force-field parameters")
    sub.add_parser("simulate",     help="Run MD simulation")
    sub.add_parser("analyze",      help="Summarize production log")
    sub.add_parser("run",          help="Full pipeline (build+param+sim+analyze)")

    args = parser.parse_args(argv)
    _setup_logging(args.verbose)

    dispatch = {
        "build":        cmd_build,
        "parameterize": cmd_parameterize,
        "simulate":     cmd_simulate,
        "analyze":      cmd_analyze,
        "run":          cmd_run,
    }
    dispatch[args.command](args)


if __name__ == "__main__":
    main()
