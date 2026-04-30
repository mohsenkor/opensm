"""
File parsers for Molden frequency files, multi-frame XYZ, velocity files,
and OpenQP transition tables.
"""

import re
import numpy as np

from .constants import Z_TO_ELEMENT


def parse_molden(filepath):
    """
    Parse a Molden frequency file.

    Returns
    -------
    atoms : list of str
        Element symbols.
    coords : ndarray, shape (natoms, 3)
        Equilibrium coordinates in Bohr.
    freqs : ndarray, shape (nmodes,)
        Vibrational frequencies in cm⁻¹.
    modes : ndarray, shape (nmodes, natoms, 3)
        Cartesian normal-mode displacement vectors (Bohr).
    """
    with open(filepath, "r") as f:
        lines = f.read().splitlines()

    atoms = []
    coords = []
    freqs = []
    modes = []

    section = None
    current_mode = []

    for line in lines:
        stripped = line.strip()

        # Detect section headers
        if stripped.upper() == "[FREQ]":
            section = "freq"
            continue
        elif stripped.upper() == "[FR-COORD]":
            section = "fr-coord"
            continue
        elif stripped.upper() == "[FR-NORM-COORD]":
            section = "fr-norm-coord"
            continue
        elif stripped.startswith("["):
            if section == "fr-norm-coord" and current_mode:
                modes.append(current_mode)
                current_mode = []
            section = None
            continue

        if not stripped:
            continue

        if section == "freq":
            try:
                freqs.append(float(stripped))
            except ValueError:
                pass

        elif section == "fr-coord":
            parts = stripped.split()
            if len(parts) >= 4:
                sym = parts[0]
                if sym.isdigit():
                    sym = Z_TO_ELEMENT.get(int(sym), sym)
                else:
                    sym = sym.capitalize()
                atoms.append(sym)
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        elif section == "fr-norm-coord":
            if stripped.lower().startswith("vibration"):
                if current_mode:
                    modes.append(current_mode)
                current_mode = []
            else:
                parts = stripped.split()
                if len(parts) >= 3:
                    current_mode.append(
                        [float(parts[0]), float(parts[1]), float(parts[2])]
                    )

    # Flush last mode
    if current_mode:
        modes.append(current_mode)

    natoms = len(atoms)
    coords = np.array(coords, dtype=float)
    freqs = np.array(freqs, dtype=float)
    modes = np.array(modes, dtype=float)

    if modes.shape[0] != len(freqs):
        raise ValueError(
            f"Molden parse error: {len(freqs)} frequencies but {modes.shape[0]} modes"
        )
    if modes.shape[1] != natoms:
        raise ValueError(
            f"Molden parse error: {natoms} atoms but mode vectors have {modes.shape[1]} entries"
        )

    return atoms, coords, freqs, modes


def parse_refxyz(filepath):
    """Parse a reference XYZ file and return element symbols."""
    with open(filepath, "r") as f:
        lines = f.read().splitlines()

    atoms = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            try:
                float(parts[1]); float(parts[2]); float(parts[3])
                sym = parts[0]
                if sym.isdigit():
                    sym = Z_TO_ELEMENT.get(int(sym), sym)
                atoms.append(sym)
            except ValueError:
                continue
    return atoms


def read_xyz_frames(filepath):
    """Read a multi-frame XYZ file.
    Returns list of (atoms, coords_array) tuples."""
    with open(filepath, "r") as f:
        lines = f.read().splitlines()

    frames = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            natoms = int(line)
        except ValueError:
            i += 1
            continue

        i += 1   # skip comment line
        if i < len(lines):
            i += 1  # move past comment

        atoms = []
        coords = []
        for _ in range(natoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            if len(parts) >= 4:
                atoms.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            i += 1

        if len(atoms) == natoms:
            frames.append((atoms, np.array(coords, dtype=float)))

    return frames


def read_velo_frames(filepath):
    """Read a Wigner velocity file.
    Returns list of (natoms, 3) numpy arrays."""
    with open(filepath, "r") as f:
        lines = f.read().splitlines()

    frames = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if "[Bohr" in line or "[bohr" in line:
            i += 1
            velo = []
            while i < len(lines):
                vline = lines[i].strip()
                if not vline or "[Bohr" in vline or "[bohr" in vline:
                    break
                parts = vline.split()
                if len(parts) >= 3:
                    try:
                        velo.append([float(parts[0]), float(parts[1]), float(parts[2])])
                    except ValueError:
                        break
                i += 1
            if velo:
                frames.append(np.array(velo, dtype=float))
        else:
            i += 1

    return frames


def read_openqp_transitions(filename, initial_state=1):
    """
    Parse an OpenQP log/output file for excitation energies and
    oscillator strengths.

    Returns list of (energy_ev, osc_strength, final_state) tuples.
    """
    transitions = []
    pattern = re.compile(
        r"^\s*(\d+)\s*->\s*(\d+)\s+"
        r"([-+]?\d+\.\d+)\s+"
        r"[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+"
        r"[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+"
        r"([-+]?\d+\.\d+)"
    )

    with open(filename, "r", errors="ignore") as f:
        for line in f:
            match = pattern.match(line)
            if match:
                istate = int(match.group(1))
                fstate = int(match.group(2))
                energy_ev = float(match.group(3))
                osc = float(match.group(4))
                if istate == initial_state:
                    transitions.append((energy_ev, osc, fstate))

    return transitions
