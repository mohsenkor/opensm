"""
Input file builders for OpenQP and PyRAI2MD.
"""

import os
import importlib.util
from pathlib import Path

from .utils import xyz_to_system_block


def _find_oqp_path():
    """Return the oqp package directory, or '' if not found."""
    spec = importlib.util.find_spec("oqp")
    if spec is not None and spec.origin is not None:
        return str(Path(spec.origin).parent)
    return ""


def sections_to_inp(sections, xyz=None):
    """Convert a dict of sections → OpenQP .inp / .openqp file content."""
    lines = []
    for section_name, params in sections.items():
        lines.append(f"[{section_name}]")
        if section_name == "input" and xyz is not None:
            system_block = xyz_to_system_block(xyz)
            lines.append(f"system=\n{system_block}")
        for key, value in params.items():
            if isinstance(value, bool):
                value = str(value).lower()
            lines.append(f"{key}={value}")
        lines.append("")
    return "\n".join(lines)


def namd_dict_to_control(namd_params, title):
    """
    Convert the namd dict → PyRAI2MD control file content.

    The 'openqp' key in &openqp is resolved in this order:
      1. Explicit value already set in namd_params
      2. oqp package directory (auto-detected via importlib)
      3. $OPENQP_ROOT environment variable
    Raises EnvironmentError if none of the above succeeds.
    """
    section_order = ["control", "molecule", "openqp", "md"]
    section_labels = {
        "control": "&CONTROL",
        "molecule": "&MOLECULE",
        "openqp": "&openqp",
        "md": "&MD",
    }

    lines = []
    for sec in section_order:
        if sec not in namd_params:
            continue

        lines.append(section_labels[sec])

        if sec == "control":
            lines.append(f"title    {title}")

        params = namd_params[sec]
        for key, value in params.items():
            if isinstance(value, bool):
                value = str(value).lower()

            if sec == "openqp" and key == "openqp" and (value == "" or value is None):
                value = (
                    _find_oqp_path()
                    or os.environ.get("OPENQP_ROOT", "")
                )
                if not value:
                    raise EnvironmentError(
                        "Cannot locate the oqp package. Either:\n"
                        "  1. Install OpenQP via pip (pip install pyopenqp), or\n"
                        "  2. Set the $OPENQP_ROOT environment variable, or\n"
                        "  3. Pass the path explicitly:\n"
                        '     manager.set_namd_config("openqp", "openqp", "/path/to/oqp")'
                    )

            if value == "" or value is None:
                lines.append(f"{key}")
            else:
                lines.append(f"{key} {value}")

        lines.append("")

    return "\n".join(lines)
