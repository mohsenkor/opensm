"""
Utility functions: file resolution, dependency checking, XYZ helpers.
"""

import os
from pathlib import Path

from .constants import ELEMENT_TO_Z


def resolve_content(value):
    """If *value* is a path to an existing file, read its text.
    Otherwise treat as a raw string."""
    if isinstance(value, (str, Path)):
        s = str(value)
        if "\n" not in s and len(s) < 260:
            p = Path(s)
            try:
                if p.is_file():
                    return p.read_text(encoding="utf-8"), True
            except OSError:
                pass
    return str(value), False


def check_package(name):
    """Check if a Python package is importable.
    Returns (available, version_or_error)."""
    try:
        mod = __import__(name)
        version = getattr(mod, "__version__", "installed")
        return True, version
    except ImportError:
        return False, "not installed"


# Which packages each mode requires
MODE_DEPENDENCIES = {
    "qp":        ["oqp"],
    "qmmm":      ["oqp", "openmm"],
    "namd":      ["oqp", "PyRAI2MD"],
    "namd_qmmm": ["oqp", "openmm", "PyRAI2MD"],
}


def check_all_dependencies():
    """Check all packages and return a dict of results."""
    packages = {
        "oqp":        "OpenQP      (QM engine)",
        "openmm":     "OpenMM      (MM engine)",
        "PyRAI2MD":   "PyRAI2MD    (NAMD driver)",
    }
    results = {}
    for pkg, label in packages.items():
        available, info = check_package(pkg)
        results[pkg] = {"label": label, "available": available, "info": info}
    return results


def check_mode_dependencies(mode):
    """Check that required packages for a given mode are installed.
    Raises ImportError if any are missing."""
    required = MODE_DEPENDENCIES.get(mode, [])
    missing = []
    for pkg in required:
        available, _ = check_package(pkg)
        if not available:
            missing.append(pkg)
    if missing:
        pkg_names = {"oqp": "OpenQP", "openmm": "OpenMM", "PyRAI2MD": "PyRAI2MD"}
        names = [f"{pkg_names.get(p, p)} ({p})" for p in missing]
        raise ImportError(
            f"Mode '{mode}' requires the following packages that are not installed:\n"
            + "\n".join(f"  - {n}" for n in names)
            + "\n\nUse SimulationManager.check_dependencies() to see all package status."
        )


def xyz_to_system_block(xyz):
    """Convert XYZ string/file → OpenQP system= lines with atomic numbers."""
    text, _ = resolve_content(xyz)
    system_lines = []
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) >= 4:
            try:
                float(parts[1]); float(parts[2]); float(parts[3])
            except ValueError:
                continue
            z = ELEMENT_TO_Z.get(parts[0], parts[0])
            system_lines.append(f"   {z}   {parts[1]}   {parts[2]}   {parts[3]}")
    return "\n".join(system_lines)
