"""OpenMM simulation: minimization, NVT equilibration, NPT production."""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from .config import Config

log = logging.getLogger(__name__)


# ── Platform selection ────────────────────────────────────────────────


def _select_platform(name: str) -> mm.Platform:
    """Return the requested (or fastest available) OpenMM platform."""
    if name.lower() != "auto":
        return mm.Platform.getPlatformByName(name)
    for candidate in ("CUDA", "OpenCL", "CPU", "Reference"):
        try:
            plat = mm.Platform.getPlatformByName(candidate)
            log.info("Selected platform: %s", plat.getName())
            return plat
        except Exception:
            continue
    raise RuntimeError("No OpenMM platform available")


# ── Reporter helpers ──────────────────────────────────────────────────


def _add_equilibration_reporters(sim: app.Simulation, config: Config,
                                 log_path: Path) -> None:
    interval = config.simulation.equilibrate.log_interval
    sim.reporters.append(app.StateDataReporter(
        str(log_path), interval,
        step=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, speed=True,
    ))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, interval,
        step=True, temperature=True, potentialEnergy=True, speed=True,
    ))


def _add_production_reporters(sim: app.Simulation, config: Config,
                              sim_dir: Path) -> None:
    prod = config.simulation.production
    sim.reporters.append(app.DCDReporter(
        str(sim_dir / "trajectory.dcd"), prod.dcd_interval,
    ))
    sim.reporters.append(app.StateDataReporter(
        str(sim_dir / "production.csv"), prod.log_interval,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, density=True, speed=True,
    ))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, prod.log_interval,
        step=True, temperature=True, density=True, speed=True,
    ))
    sim.reporters.append(app.CheckpointReporter(
        str(sim_dir / "checkpoint.chk"), prod.checkpoint_interval,
    ))


# ── Simulation phases ────────────────────────────────────────────────


def _minimize(simulation: app.Simulation, config: Config,
              sim_dir: Path) -> None:
    log.info("Energy minimization (max %d iterations) ...",
             config.simulation.minimize.max_iterations)
    simulation.minimizeEnergy(
        maxIterations=config.simulation.minimize.max_iterations,
    )
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    log.info("Minimized potential energy: %s",
             state.getPotentialEnergy())

    with (sim_dir / "minimized.pdb").open("w") as fh:
        app.PDBFile.writeFile(
            simulation.topology, state.getPositions(), fh,
        )


def _equilibrate_nvt(simulation: app.Simulation, config: Config,
                     sim_dir: Path) -> None:
    """Short NVT run to thermalise the system."""
    steps = config.simulation.equilibrate.steps
    log.info("NVT equilibration: %d steps (%.1f ps) at %.1f K ...",
             steps,
             steps * config.simulation.timestep_fs / 1000,
             config.simulation.temperature_kelvin)

    _add_equilibration_reporters(simulation, config,
                                 sim_dir / "equilibration.csv")
    simulation.step(steps)

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    log.info("Post-equilibration potential: %s", state.getPotentialEnergy())

    with (sim_dir / "equilibrated.pdb").open("w") as fh:
        app.PDBFile.writeFile(
            simulation.topology, state.getPositions(), fh,
        )

    # Clear reporters before production
    simulation.reporters.clear()


def _production_npt(simulation: app.Simulation, system: mm.System,
                    config: Config, sim_dir: Path) -> None:
    """NPT production run with Monte Carlo barostat."""
    sim_cfg = config.simulation
    steps = sim_cfg.production.steps

    # Add barostat for NPT
    barostat = mm.MonteCarloBarostat(
        sim_cfg.pressure_bar * unit.bar,
        sim_cfg.temperature_kelvin * unit.kelvin,
        25,  # attempt frequency (steps)
    )
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)

    log.info("NPT production: %d steps (%.1f ps) at %.1f K / %.1f bar ...",
             steps,
             steps * sim_cfg.timestep_fs / 1000,
             sim_cfg.temperature_kelvin,
             sim_cfg.pressure_bar)

    _add_production_reporters(simulation, config, sim_dir)
    simulation.step(steps)

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    log.info("Final potential energy: %s", state.getPotentialEnergy())

    with (sim_dir / "final.pdb").open("w") as fh:
        app.PDBFile.writeFile(
            simulation.topology, state.getPositions(), fh,
        )


# ── Public API ────────────────────────────────────────────────────────


def run_simulation(system: mm.System, topology: app.Topology,
                   positions, config: Config) -> None:
    """Run the full simulation protocol: minimize -> NVT -> NPT."""
    sim_cfg = config.simulation
    sim_dir = Path(config.output.directory) / "simulate"
    sim_dir.mkdir(parents=True, exist_ok=True)

    # ── Integrator ────────────────────────────────────────────────────
    integrator = mm.LangevinMiddleIntegrator(
        sim_cfg.temperature_kelvin * unit.kelvin,
        sim_cfg.friction_per_ps / unit.picoseconds,
        sim_cfg.timestep_fs * unit.femtoseconds,
    )

    platform = _select_platform(sim_cfg.platform)
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)

    log.info("System: %d atoms | Platform: %s",
             topology.getNumAtoms(), platform.getName())

    # ── Run protocol ──────────────────────────────────────────────────
    _minimize(simulation, config, sim_dir)
    _equilibrate_nvt(simulation, config, sim_dir)
    _production_npt(simulation, system, config, sim_dir)

    log.info("Simulation complete.  Output in %s/", sim_dir)
