"""
Input file builders for OpenQP and PyRAI2MD.
"""

import os

from .utils import xyz_to_system_block


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

    The 'openqp' key in &openqp section is auto-filled from
    the $OPENQP_ROOT environment variable if not set explicitly.
    """
    section_order = ["control", "molecule", "openqp", "md"]
    section_labels = {
        "control": "&CONTROL",
        "molecule": "&MOLECULE",
        "openqp": "&openqp",
        "md": "&MD",
    }

    openqp_root = os.environ.get("OPENQP_ROOT", "")

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
                if openqp_root:
                    value = openqp_root
                else:
                    raise EnvironmentError(
                        "OpenQP path not set. Either:\n"
                        "  1. Set the $OPENQP_ROOT environment variable, or\n"
                        "  2. Pass it explicitly:\n"
                        '     manager.set_namd_config("openqp", "openqp", "/path/to/openqp")'
                    )

            if value == "" or value is None:
                lines.append(f"{key}")
            else:
                lines.append(f"{key} {value}")

        lines.append("")

    return "\n".join(lines)
