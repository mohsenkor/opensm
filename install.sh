#!/usr/bin/env bash
# install.sh — Full installation of opensm and its NAMD dependencies.
#
# Usage:
#   bash install.sh        # pip install + compile fssh library
#
# Requirements:
#   - Python 3.8+ with pip
#   - git
#   - A C/Fortran compiler (gfortran) for the fssh library

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=============================================="
echo " OpenSM installation"
echo "=============================================="
echo " Repo: $REPO_DIR"
echo "=============================================="

# ── 0. Upgrade pip ───────────────────────────────
echo ""
echo "[0/2] Upgrading pip..."
pip install --upgrade pip
echo "  Done."

# ── 1. Install opensm + all dependencies ─────────
echo ""
echo "[1/2] Installing opensm and dependencies..."
pip install "$REPO_DIR"
echo "  Done."

# ── 2. Compile fssh library ───────────────────────
echo ""
echo "[2/2] Compiling PyRAI2MD fssh library (pyrai2md update)..."
pyrai2md update
echo "  Done."

echo ""
echo "=============================================="
echo " Installation complete."
echo ""
echo " Verify with:"
echo "   python -c \"from opensm import SimulationManager\""
echo "   python -c \"from opensm import SimulationManager; SimulationManager.check_dependencies()\""
echo "=============================================="
