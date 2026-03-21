"""Post-simulation analysis: parse logs, compute summary statistics."""
from __future__ import annotations

import csv
import logging
import statistics
from pathlib import Path

from .config import Config

log = logging.getLogger(__name__)


def _parse_csv(path: Path) -> dict[str, list[float]]:
    """Read a StateDataReporter CSV into ``{column: [values]}``."""
    columns: dict[str, list[float]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            for key, val in row.items():
                key = key.strip().strip('"')
                columns.setdefault(key, []).append(float(val))
    return columns


def _stats(values: list[float]) -> str:
    if not values:
        return "n/a"
    mn = min(values)
    mx = max(values)
    avg = statistics.mean(values)
    std = statistics.stdev(values) if len(values) > 1 else 0.0
    return f"mean={avg:.2f}  std={std:.2f}  min={mn:.2f}  max={mx:.2f}"


def summarize(config: Config) -> None:
    """Print summary statistics from the production log."""
    sim_dir = Path(config.output.directory) / "simulate"
    prod_csv = sim_dir / "production.csv"

    if not prod_csv.exists():
        log.warning("Production log not found at %s", prod_csv)
        return

    data = _parse_csv(prod_csv)

    print("\n=== Production run summary ===\n")
    for label, key in [
        ("Temperature (K)",        "Temperature (K)"),
        ("Potential Energy (kJ/mol)", 'Potential Energy (kJ/mole)'),
        ("Total Energy (kJ/mol)",  'Total Energy (kJ/mole)'),
        ("Density (g/mL)",         "Density (g/mL)"),
    ]:
        vals = data.get(key, [])
        if vals:
            print(f"  {label:32s}  {_stats(vals)}")
        else:
            # Try without the exact key (header names vary slightly)
            for k, v in data.items():
                if label.split("(")[0].strip().lower() in k.lower():
                    print(f"  {label:32s}  {_stats(v)}")
                    break
            else:
                print(f"  {label:32s}  (not found in log)")

    if "Speed (ns/day)" in data:
        speeds = data["Speed (ns/day)"]
        print(f"\n  {'Performance':32s}  {_stats(speeds)} ns/day")

    print()
