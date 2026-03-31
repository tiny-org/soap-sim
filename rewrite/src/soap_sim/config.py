"""Configuration dataclasses and TOML loading."""
from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

MW_WATER = 18.015


@lru_cache(maxsize=32)
def _formal_charge(smiles: str) -> int:
    """Net formal charge of a molecule from its SMILES."""
    from rdkit import Chem
    return Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))


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
    residue_name: str = ""

    @property
    def resname(self) -> str:
        if self.residue_name:
            return self.residue_name.upper()
        return self.name[:3].upper()


@dataclass(frozen=True)
class CounterionSpec:
    """One counterion species with its fraction of the total."""
    smiles: str
    fraction: float = 1.0  # fraction of total counterions (0-1)

    @property
    def resname(self) -> str:
        from rdkit import Chem
        elem = Chem.MolFromSmiles(self.smiles).GetAtomWithIdx(0).GetSymbol()
        return elem[:3].upper()


@dataclass(frozen=True)
class SystemConfig:
    solutes: tuple[SoluteSpec, ...] = ()
    counterions: tuple[CounterionSpec, ...] = (CounterionSpec("[Na+]"),)
    num_water: int = 0
    target_density: float = 1.0
    box_dimensions: tuple[float, float, float] | None = None

    @property
    def total_counterions(self) -> int:
        """Total counterions needed to neutralize the solute charge."""
        total_solute_charge = sum(
            s.count * _formal_charge(s.smiles) for s in self.solutes
        )
        if not self.counterions or total_solute_charge == 0:
            return 0
        ci_charge = _formal_charge(self.counterions[0].smiles)
        if ci_charge == 0:
            return 0
        return abs(total_solute_charge) // abs(ci_charge)

    @property
    def counterion_counts(self) -> list[tuple[CounterionSpec, int]]:
        """List of (spec, count) for each counterion type."""
        total = self.total_counterions
        if total == 0:
            return []
        result = []
        assigned = 0
        for i, ci in enumerate(self.counterions):
            if i == len(self.counterions) - 1:
                n = total - assigned  # last one gets the remainder
            else:
                n = round(total * ci.fraction)
            result.append((ci, n))
            assigned += n
        return result

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
        for ci, n in self.counterion_counts:
            mass += n * _mw_with_h(ci.smiles)
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
    barostat: str = "isotropic"
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


def _num_water_from_fraction(solutes, counterions, weight_fraction):
    """Compute water count so that water is *weight_fraction* of total mass."""
    soap_mass = 0.0
    for s in solutes:
        soap_mass += s.count * _mw_with_h(s.smiles)
    # Approximate: use average counterion MW weighted by fraction
    avg_ci_mw = sum(ci.fraction * _mw_with_h(ci.smiles) for ci in counterions)
    total_anions = sum(s.count for s in solutes
                       if _formal_charge(s.smiles) < 0)
    soap_mass += total_anions * avg_ci_mw
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

    # Counterions: list or legacy single string
    ci_raw = sys_raw.pop("counterions", None)
    if ci_raw is not None:
        counterions = tuple(CounterionSpec(**c) for c in ci_raw)
        frac_sum = sum(c.fraction for c in counterions)
        if abs(frac_sum - 1.0) > 0.01:
            raise ValueError(
                f"Counterion fractions must sum to 1.0, got {frac_sum:.3f}")
    else:
        ci_smiles = sys_raw.pop("counterion_smiles", "[Na+]")
        counterions = (CounterionSpec(ci_smiles),)

    if "num_water" in sys_raw:
        num_water = sys_raw.pop("num_water")
    elif "water_weight_fraction" in sys_raw:
        wf = sys_raw.pop("water_weight_fraction")
        num_water = _num_water_from_fraction(solutes, counterions, wf)
    else:
        raise ValueError("Config must set either num_water or water_weight_fraction")

    system = SystemConfig(
        solutes=solutes,
        counterions=counterions,
        num_water=num_water,
        target_density=sys_raw.pop("target_density", 1.0),
        box_dimensions=box_dim,
    )

    if sys_raw:
        raise ValueError(f"Unknown [system] config keys: {list(sys_raw.keys())}")

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
