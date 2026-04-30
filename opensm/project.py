"""
Project persistence: save / load manager state across Python sessions.
"""

import json
from copy import deepcopy
from pathlib import Path


def save_project(manager, filepath=None):
    """
    Save the full manager state to a JSON file.

    Parameters
    ----------
    manager : SimulationManager
        The manager instance to save.
    filepath : str or Path, optional
        Where to save (default: <folder>/<title>.project.json).

    Returns
    -------
    str
        Path to the saved project file.
    """
    proj_path = Path(filepath) if filepath else (
        Path(manager.folder) / f"{manager.title}.project.json"
    )
    proj_path.parent.mkdir(parents=True, exist_ok=True)

    state = {
        "_version": 1,
        "mode": manager.mode,
        "folder": manager.folder,
        "title": manager.title,
        "template_source": manager._template_source,
        "sections": manager._sections,
        "namd_params": manager._namd_params,
        "xyz": manager._xyz,
        "pdb": manager._pdb,
        "qm_atoms": manager._qm_atoms,
        "generated": manager._generated,
        "configs": manager.configs,
        "workflow": getattr(manager, "_workflow", {}),
    }

    with open(proj_path, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2, ensure_ascii=False, default=str)

    print(f"Project saved: {proj_path}")
    return str(proj_path)


def load_project_state(filepath):
    """
    Load raw project state from a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to the .project.json file.

    Returns
    -------
    dict
        Raw state dictionary.
    """
    proj_path = Path(filepath)
    if not proj_path.is_file():
        raise FileNotFoundError(f"Project file not found: {filepath}")

    with open(proj_path, "r", encoding="utf-8") as f:
        return json.load(f)


def show_workflow(workflow, title, mode, folder):
    """
    Print a summary of completed workflow steps.

    Parameters
    ----------
    workflow : dict
        The _workflow dict from the manager.
    title, mode, folder : str
        Manager metadata for display.
    """
    print(f"Workflow status for '{title}':")
    print(f"  Mode: {mode} | Folder: {folder}")
    print("-" * 60)

    if "wigner" in workflow:
        w = workflow["wigner"]
        print(f"  \u2713 wigner_sampling     "
              f"{w.get('nsamples', '?')} samples, "
              f"{w.get('nmodes', '?')} modes")
        print(f"      geom \u2192 {w.get('geom_file', '?')}")
        print(f"      velo \u2192 {w.get('velo_file', '?')}")
    else:
        print("  \u2717 wigner_sampling     (not yet run)")

    if "absorption_batch" in workflow:
        ab = workflow["absorption_batch"]
        print(f"  \u2713 absorption_batch    "
              f"{ab.get('n_samples', '?')} samples")
        print(f"      folder \u2192 {ab.get('base_folder', '?')}")
    else:
        print("  \u2717 absorption_batch    (not yet run)")

    if "absorption_spectrum" in workflow:
        sp = workflow["absorption_spectrum"]
        print(f"  \u2713 absorption_spectrum "
              f"{sp.get('n_converged', '?')} converged, "
              f"{sp.get('n_failed', '?')} failed")
        print(f"      data \u2192 {sp.get('spectrum_file', '?')}")
        if sp.get("plot_file"):
            print(f"      plot \u2192 {sp['plot_file']}")
    else:
        print("  \u2717 absorption_spectrum (not yet run \u2014 "
              "run after batch jobs finish)")

    if "selected" in workflow:
        sel = workflow["selected"]
        print(f"  \u2713 select_wigner_window "
              f"{sel.get('n_selected', '?')} / {sel.get('n_total', '?')} "
              f"in [{sel.get('e_min', '?')}, {sel.get('e_max', '?')}] eV")
        print(f"      geom \u2192 {sel.get('geom_file', '?')}")
        print(f"      velo \u2192 {sel.get('velo_file', '?')}")
    else:
        print("  \u2717 select_wigner_window (not yet run)")


def get_workflow(workflow, step=None):
    """
    Retrieve saved workflow results.

    Parameters
    ----------
    workflow : dict
        The _workflow dict.
    step : str, optional
        One of 'wigner', 'absorption_batch', 'absorption_spectrum',
        'selected'.  If None, returns the entire dict.

    Returns
    -------
    dict or None
    """
    if step is None:
        return workflow

    data = workflow.get(step)
    if data is None:
        return None

    data = deepcopy(data)

    # Restore integer keys for all_transitions
    if step == "absorption_spectrum" and "all_transitions" in data:
        data["all_transitions"] = {
            int(k): v for k, v in data["all_transitions"].items()
        }

    return data
