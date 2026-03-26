"""Configuration dataclasses and TOML loading."""
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

MW_WATER = 18.015


@lru_cache(maxsize=32)
def _mw_with_h(smiles: str) -> float:
    """Molecular weight including implicit hydrogens (uses RDKit)."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    return Descriptors.MolWt(Chem.AddHs(Chem.MolFromSmiles(smiles)))


# ── Dataclasses ───────────────────────────────────────────────────────


@dataclass(frozen=True)
class SoluteSpec:
    """One solute species (e.g. one fatty-acid soap)."""
    name: str
    smiles: str
    count: int

    @property
    def residue_name(self) -> str:
        return self.name[:3].upper()


@dataclass(frozen=True)
class SystemConfig:
    solutes: tuple[SoluteSpec, ...] = ()
    counterion_smiles: str = "[Na+]"
    num_water: int = 0
    target_density: float = 1.0
    box_dimensions: tuple[float, float, float] | None = None

    @property
    def num_counterions(self) -> int:
        return sum(s.count for s in self.solutes)

    @property
    def box_angstrom(self) -> tuple[float, float, float]:
        if self.box_dimensions is not None:
            return self.box_dimensions
        total_mass_g = self._total_mass_amu() * 1.6605e-24
        volume_ang3 = (total_mass_g / self.target_density) * 1e24
        side = round(volume_ang3 ** (1 / 3) * 1.05, 1)
        return (side, side, side)

    def _total_mass_amu(self) -> float:
        mass = self.num_water * MW_WATER
        mass += self.num_counterions * _mw_with_h(self.counterion_smiles)
        for s in self.solutes:
            mass += s.count * _mw_with_h(s.smiles)
        return mass


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


def _num_water_from_fraction(solutes, counterion_smiles, weight_fraction):
    """Compute water count so that water is *weight_fraction* of total mass."""
    soap_mass = sum(
        s.count * (_mw_with_h(s.smiles) + _mw_with_h(counterion_smiles))
        for s in solutes
    )
    mass_water = soap_mass * weight_fraction / (1.0 - weight_fraction)
    return round(mass_water / MW_WATER)


def load_config(path: Path) -> Config:
    """Read a TOML config file and return a validated :class:`Config`."""
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)

    # ── System ────────────────────────────────────────────────────────
    sys_raw = dict(raw.get("system", {}))

    box_dim = sys_raw.pop("box_dimensions", None)
    if box_dim is not None:
        box_dim = tuple(box_dim)

    solutes = tuple(
        SoluteSpec(**s) for s in sys_raw.pop("solutes", [])
    )
    if not solutes:
        raise ValueError("Config must define at least one [[system.solutes]] entry")

    counterion = sys_raw.pop("counterion_smiles", "[Na+]")

    if "num_water" in sys_raw:
        num_water = sys_raw.pop("num_water")
    elif "water_weight_fraction" in sys_raw:
        wf = sys_raw.pop("water_weight_fraction")
        num_water = _num_water_from_fraction(solutes, counterion, wf)
    else:
        raise ValueError("Config must set either num_water or water_weight_fraction")

    system = SystemConfig(
        solutes=solutes,
        counterion_smiles=counterion,
        num_water=num_water,
        target_density=sys_raw.pop("target_density", 1.0),
        box_dimensions=box_dim,
    )

    # ── Force field ───────────────────────────────────────────────────
    forcefield = ForceFieldConfig(**raw.get("forcefield", {}))

    # ── Simulation ────────────────────────────────────────────────────
    sim_raw = dict(raw.get("simulation", {}))
    simulation = SimulationConfig(
        **{k: v for k, v in sim_raw.items()
           if k not in ("minimize", "equilibrate", "production")},
        minimize=MinimizeConfig(**sim_raw.get("minimize", {})),
        equilibrate=EquilibrateConfig(**sim_raw.get("equilibrate", {})),
        production=ProductionConfig(**sim_raw.get("production", {})),
    )

    output = OutputConfig(**raw.get("output", {}))
    return Config(system=system, forcefield=forcefield,
                  simulation=simulation, output=output)
