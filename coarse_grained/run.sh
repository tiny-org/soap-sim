#!/usr/bin/env bash
# Run the MARTINI CG simulation: minimize -> NVT -> NPT
set -euo pipefail
cd "$(dirname "$0")"

OUT=output
TOP="$OUT/topol.top"
GRO="$OUT/system.gro"

if [ ! -f "$GRO" ] || [ ! -f "$TOP" ]; then
    echo "ERROR: $GRO or $TOP not found.  Run 'python build.py' first."
    exit 1
fi

echo "=== Energy minimization ==="
gmx grompp -f mdp/em.mdp -c "$GRO" -p "$TOP" -o "$OUT/em.tpr" -maxwarn 2
gmx mdrun -deffnm "$OUT/em" -v

echo ""
echo "=== NVT equilibration (1 ns) ==="
gmx grompp -f mdp/nvt.mdp -c "$OUT/em.gro" -p "$TOP" -o "$OUT/nvt.tpr" -maxwarn 2
gmx mdrun -deffnm "$OUT/nvt" -v

echo ""
echo "=== NPT production (1 us) ==="
gmx grompp -f mdp/npt.mdp -c "$OUT/nvt.gro" -p "$TOP" -t "$OUT/nvt.cpt" -o "$OUT/npt.tpr" -maxwarn 2
gmx mdrun -deffnm "$OUT/npt" -v

echo ""
echo "=== Done ==="
echo "Trajectory: $OUT/npt.xtc"
echo "Energy:     $OUT/npt.edr"
echo ""
echo "View in PyMOL:"
echo "  gmx trjconv -f $OUT/npt.xtc -s $OUT/npt.tpr -o $OUT/traj.pdb -pbc mol <<< 0"
echo "  pymol $OUT/traj.pdb"
