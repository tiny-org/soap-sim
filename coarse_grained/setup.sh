#!/usr/bin/env bash
# Download MARTINI 2.2 force field files and verify GROMACS is available.
set -euo pipefail
cd "$(dirname "$0")"

echo "=== Checking GROMACS ==="
if ! command -v gmx &>/dev/null; then
    echo "GROMACS not found.  Install with:"
    echo "  micromamba install -c conda-forge gromacs"
    exit 1
fi
gmx --version | head -3

echo ""
echo "=== Downloading MARTINI 2.2 force field ==="
mkdir -p forcefield
FF=forcefield/martini_v2.2.itp
IONS=forcefield/martini_v2.0_ions.itp

if [ ! -f "$FF" ]; then
    curl -fSL -o "$FF" \
        "http://cgmartini.nl/images/parameters/ITP/martini_v2.2.itp"
    echo "Downloaded $FF"
else
    echo "$FF already exists, skipping."
fi

if [ ! -f "$IONS" ]; then
    curl -fSL -o "$IONS" \
        "http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp"
    echo "Downloaded $IONS"
else
    echo "$IONS already exists, skipping."
fi

echo ""
echo "=== Setup complete ==="
echo "Run:  python build.py"
echo "Then: bash run.sh"
