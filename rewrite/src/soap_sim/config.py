"""Configuration dataclasses and TOML loading."""
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# Molecular weights (g/mol)
MW_SODIUM_STEARATE = 306.45  # Na+ + C18H35O2-
MW_WATER = 18.015


def _num_water(num_stearate: int, weight_fraction: float) -> int:
    """Compute water molecules needed for a target weight fraction."""
    mass_stearate = num_stearate * MW_SODIUM_STEARATE
    mass_water = mass_stearate * (1.0 - weight_fraction) / weight_fraction
    return round(mass_water / MW_WATER)


def _box_size_angstrom(num_stearate: int, num_water: int,
                       target_density: float) -> float:
    """Estimate cubic box side length (A) from density."""
    total_mass_amu = num_stearate * MW_SODIUM_STEARATE + num_water * MW_WATER
    total_mass_g = total_mass_amu * 1.6605e-24
    volume_ang3 = (total_mass_g / target_density) * 1e24
    # 5 % padding so PACKMOL can place molecules comfortably
    return round(volume_ang3 ** (1 / 3) * 1.05, 1)


# ── Dataclasses ───────────────────────────────────────────────────────


@dataclass(frozen=True)
class SystemConfig:
    num_stearate: int = 50
    weight_fraction: float = 0.50
    target_density: float = 1.0
    box_dimensions: Optional[tuple[float, float, float]] = None  # A, override

    @property
    def num_sodium(self) -> int:
        return self.num_stearate

    @property
    def num_water(self) -> int:
        return _num_water(self.num_stearate, self.weight_fraction)

    @property
    def box_size_angstrom(self) -> float:
        """Cubic box side length (used when box_dimensions is None)."""
        return _box_size_angstrom(self.num_stearate, self.num_water,
                                  self.target_density)

    @property
    def box_angstrom(self) -> tuple[float, float, float]:
        """(x, y, z) box dimensions in Angstroms."""
        if self.box_dimensions is not None:
            return self.box_dimensions
        s = self.box_size_angstrom
        return (s, s, s)


@dataclass(frozen=True)
class ForceFieldConfig:
    small_molecule: str = "openff-2.2.0.offxml"
    water_model: str = "tip3p.xml"


@dataclass(frozen=True)
class MinimizeConfig:
    max_iterations: int = 1000


@dataclass(frozen=True)
class EquilibrateConfig:
    steps: int = 50_000
    log_interval: int = 500


@dataclass(frozen=True)
class ProductionConfig:
    steps: int = 500_000
    dcd_interval: int = 1000
    log_interval: int = 1000
    checkpoint_interval: int = 10_000


@dataclass(frozen=True)
class SimulationConfig:
    temperature_celsius: float = 80.0
    pressure_bar: float = 1.0
    timestep_fs: float = 2.0
    friction_per_ps: float = 1.0
    platform: str = "auto"
    barostat: str = "isotropic"  # isotropic | anisotropic
    minimize: MinimizeConfig = field(default_factory=MinimizeConfig)
    equilibrate: EquilibrateConfig = field(default_factory=EquilibrateConfig)
    production: ProductionConfig = field(default_factory=ProductionConfig)

    @property
    def temperature_kelvin(self) -> float:
        return self.temperature_celsius + 273.15


@dataclass(frozen=True)
class OutputConfig:
    directory: str = "output"


@dataclass(frozen=True)
class Config:
    system: SystemConfig = field(default_factory=SystemConfig)
    forcefield: ForceFieldConfig = field(default_factory=ForceFieldConfig)
    simulation: SimulationConfig = field(default_factory=SimulationConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


# ── Loader ────────────────────────────────────────────────────────────


def load_config(path: Path) -> Config:
    """Read *path* and return a fully-validated :class:`Config`."""
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)

    sys_raw = dict(raw.get("system", {}))
    box_dim = sys_raw.pop("box_dimensions", None)
    if box_dim is not None:
        box_dim = tuple(box_dim)
    system = SystemConfig(**sys_raw, box_dimensions=box_dim)
    forcefield = ForceFieldConfig(**raw.get("forcefield", {}))

    sim_raw = dict(raw.get("simulation", {}))
    minimize = MinimizeConfig(**sim_raw.pop("minimize", {}))
    equilibrate = EquilibrateConfig(**sim_raw.pop("equilibrate", {}))
    production = ProductionConfig(**sim_raw.pop("production", {}))
    simulation = SimulationConfig(
        **sim_raw,
        minimize=minimize,
        equilibrate=equilibrate,
        production=production,
    )

    output = OutputConfig(**raw.get("output", {}))
    return Config(system=system, forcefield=forcefield,
                  simulation=simulation, output=output)
