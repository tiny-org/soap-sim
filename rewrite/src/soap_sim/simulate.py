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


# ── Public API ────────────────────────────────────────────────────────


def run_simulation(system: mm.System, topology: app.Topology,
                   positions, config: Config) -> None:
    """Run the full simulation protocol: minimize -> NVT -> NPT."""
    sim_cfg = config.simulation
    sim_dir = Path(config.output.directory) / "simulate"
    sim_dir.mkdir(parents=True, exist_ok=True)

    platform = _select_platform(sim_cfg.platform)

    def _make_integrator():
        return mm.LangevinMiddleIntegrator(
            sim_cfg.temperature_kelvin * unit.kelvin,
            sim_cfg.friction_per_ps / unit.picoseconds,
            sim_cfg.timestep_fs * unit.femtoseconds,
        )

    # ── Phase 1: Minimize + NVT equilibration (no barostat) ──────────
    nvt_sim = app.Simulation(topology, system, _make_integrator(), platform)
    nvt_sim.context.setPositions(positions)

    log.info("System: %d atoms | Platform: %s",
             topology.getNumAtoms(), platform.getName())

    _minimize(nvt_sim, config, sim_dir)
    _equilibrate_nvt(nvt_sim, config, sim_dir)

    # Save full state (positions, velocities, box vectors)
    nvt_state = nvt_sim.context.getState(
        getPositions=True, getVelocities=True,
    )
    del nvt_sim  # free GPU/OpenCL resources

    # ── Phase 2: NPT production (with barostat) ─────────────────────
    pressure = sim_cfg.pressure_bar * unit.bar
    temperature = sim_cfg.temperature_kelvin * unit.kelvin

    if sim_cfg.barostat == "anisotropic":
        barostat = mm.MonteCarloAnisotropicBarostat(
            mm.Vec3(pressure, pressure, pressure),
            temperature,
            True, True, True,  # scale x, y, z independently
            25,
        )
        log.info("Using anisotropic barostat (each axis scales independently)")
    else:
        barostat = mm.MonteCarloBarostat(pressure, temperature, 25)

    system.addForce(barostat)

    npt_sim = app.Simulation(topology, system, _make_integrator(), platform)
    npt_sim.context.setPositions(nvt_state.getPositions())
    npt_sim.context.setVelocities(nvt_state.getVelocities())
    npt_sim.context.setPeriodicBoxVectors(
        *nvt_state.getPeriodicBoxVectors()
    )

    steps = sim_cfg.production.steps
    log.info("NPT production: %d steps (%.1f ps) at %.1f K / %.1f bar ...",
             steps,
             steps * sim_cfg.timestep_fs / 1000,
             sim_cfg.temperature_kelvin,
             sim_cfg.pressure_bar)

    _add_production_reporters(npt_sim, config, sim_dir)
    npt_sim.step(steps)

    final_state = npt_sim.context.getState(
        getEnergy=True, getPositions=True,
    )
    log.info("Final potential energy: %s", final_state.getPotentialEnergy())

    with (sim_dir / "final.pdb").open("w") as fh:
        app.PDBFile.writeFile(
            topology, final_state.getPositions(), fh,
        )

    log.info("Simulation complete.  Output in %s/", sim_dir)
