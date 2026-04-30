"""
Built-in simulation templates (loaded from the templates/ folder) and the
internal Hessian template used by Wigner sampling.
"""

import json
from pathlib import Path

# ====================================================================
#  Hessian template for Wigner sampling (internal — not user-facing)
# ====================================================================
WIGNER_HESS_TEMPLATE = {
    "input": {
        "runtype": "hess",
        "charge": 0,
        "functional": "bhhlyp",
        "basis": "6-31g*",
        "method": "tdhf",
    },
    "scf": {
        "type": "rohf",
        "multiplicity": 3,
    },
    "tdhf": {
        "type": "mrsf",
        "nstate": 2,
    },
    "hess": {
        "state": 1,
    },
}


# ====================================================================
#  Built-in simulation templates (one JSON file per template)
# ====================================================================
TEMPLATES_DIR = Path(__file__).parent / "templates"


def _load_builtin_templates():
    templates = {}
    if TEMPLATES_DIR.is_dir():
        for f in sorted(TEMPLATES_DIR.glob("*.json")):
            try:
                data = json.loads(f.read_text(encoding="utf-8"))
                templates[f.stem] = data
            except (json.JSONDecodeError, OSError):
                pass
    return templates


BUILTIN_TEMPLATES = _load_builtin_templates()
