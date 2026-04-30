"""
OpenSM — Unified simulation manager for multi-scale simulations in OpenQP.

Usage
-----
    from opensm import SimulationManager

    manager = SimulationManager(
        mode="namd", folder="run", title="mol", xyz="mol.xyz",
    )

    # Full workflow:
    wigner = manager.wigner_sampling(nsamples=1000)
    abs_info = manager.absorption_batch(wigner['geom_file'])
    manager.save_project()

    # --- after batch jobs finish, in a new session ---
    manager = SimulationManager.load_project("run/mol.project.json")
    spec = manager.absorption_spectrum()
    selected = manager.select_wigner_window(...)
    manager.batch_from_wigner(selected['geom_file'], selected['velo_file'])
"""

from .manager import SimulationManager
from .engines import OpenQPEngine, PyRAI2MDEngine
from .analysis import TrajectoryAnalysis

__all__ = [
    "SimulationManager",
    "OpenQPEngine",
    "PyRAI2MDEngine",
    "TrajectoryAnalysis",
]
