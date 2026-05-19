# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
OpenSM is a Python package (`opensm/`) that provides a unified interface for multi-scale quantum chemistry simulations using OpenQP. It orchestrates QM, QM/MM, NAMD, and NAMD+QM/MM workflows through a single `SimulationManager` class.

**Author**: Mohsen Mazaherifar (computational chemistry researcher)
**Backends**: OpenQP (QM engine), OpenMM (MM engine), PyRAI2MD (NAMD driver)
**Method**: MRSF-TDDFT (multi-reference spin-flip TDDFT) is the primary electronic structure method

## Development Setup

```bash
# Import check
python -c "from opensm import SimulationManager"

# Check which backends are available for each mode
python -c "from opensm import SimulationManager; SimulationManager.check_dependencies()"
```

No environment variable is needed. OpenQP path is auto-detected from the installed `oqp` package.

There is no formal test suite. Development testing is done inline:
```python
from opensm import SimulationManager
manager = SimulationManager(mode="namd", folder="test_run", title="mol", xyz="mol.xyz")
manager.show_config()
# verify output files in test_run/
```

## Package Structure
```
opensm/
├── __init__.py      # Exports: SimulationManager, OpenQPEngine, PyRAI2MDEngine, TrajectoryAnalysis
├── constants.py     # ELEMENT_TO_Z, ELEMENT_TO_MASS, Z_TO_ELEMENT, physical constants (a.u.)
├── utils.py         # resolve_content(), dependency checking, xyz_to_system_block()
├── parsers.py       # parse_molden(), parse_refxyz(), read_xyz_frames(), read_velo_frames(), read_openqp_transitions()
├── builders.py      # sections_to_inp(), namd_dict_to_control(), _find_oqp_path()
├── templates.py     # BUILTIN_TEMPLATES (loaded from templates/), WIGNER_HESS_TEMPLATE
├── templates/       # Built-in JSON templates (qp_default, qmmm_default, namd_default, namd_qmmm_default)
├── engines.py       # OpenQPEngine, PyRAI2MDEngine classes
├── wigner.py        # build_hess_input(), run_hessian(), wigner_sampling()
├── spectrum.py      # build_spectrum_input(), absorption_batch(), absorption_spectrum(), select_wigner_window()
├── batch.py         # generate_run_script(), generate_job_script(), generate_master_script(), submit_jobs()
├── project.py       # save_project(), load_project_state(), show_workflow(), get_workflow()
├── analysis.py      # TrajectoryAnalysis class — parses PyRAI2MD output, computes populations/energies/hops
└── manager.py       # SimulationManager class (thin orchestrator delegating to modules)
```

## Architecture Decisions
- **manager.py is a thin orchestrator** — each method delegates to standalone functions in focused modules so they can be debugged/tested independently
- **Templates are the single source of truth** — all default configs live in `templates/*.json`; `_default_sections()` and `_default_namd()` in `manager.py` load from `BUILTIN_TEMPLATES`, there are no hardcoded dicts
- **Workflow-aware defaults** — methods auto-fill file arguments from previous workflow steps (wigner → absorption → select → batch → analyze), so after `load_project()` users don't need to manually pass paths
- **Project persistence** — `save_project()` / `load_project()` serialize full state to JSON for cross-session workflows
- **Auto-save** — `_save_workflow()` auto-updates the project file if it exists on disk after each workflow step
- **All SLURM/PBS paths must be absolute** — `batch.py` uses `Path.resolve()` for `--output`, `--error`, and `cd` paths
- **State-resolved absorption spectrum** — decomposed by final state (S₁, S₂, …) with per-state Gaussian-broadened curves
- **`TrajectoryAnalysis` is standalone** — can be constructed directly without `SimulationManager`; `manager.analyze_trajectories()` is a convenience wrapper

## Non-Obvious Cross-File Behaviors

**`resolve_content()` heuristic** (`utils.py`): Every argument that accepts a file path or inline string (xyz, pdb, openqp_content, etc.) is passed through `resolve_content()`. If the value has no newlines and is under 260 chars, it tries to open it as a file; otherwise it treats the value as literal content. This means passing a short string without newlines could accidentally be interpreted as a file path.

**`PyRAI2MDEngine._set_runtype()`** (`engines.py`): When writing NAMD input files, the engine silently overrides `runtype` in the generated `.openqp` file — setting it to `prop` for `namd` mode and `qmmm` for `namd_qmmm`. This happens regardless of what `runtype` is set to in `sections`. Do not set `runtype` in sections for NAMD modes.

**NAMD file naming** (`engines.py`, `batch.py`): For NAMD modes, the control file is written as `title` (no extension), the QM input as `title.openqp`, geometry as `title.xyz`, and velocities as `title.velo`. OpenQP modes use `title.inp` instead.

**`namd_dict_to_control()` OpenQP path resolution** (`builders.py`): If the `openqp` key inside the `&openqp` section is empty or `None`, the path is resolved in this order: (1) auto-detect via `importlib.util.find_spec("oqp")` — returns the installed `oqp` package directory; (2) `$OPENQP_ROOT` environment variable; (3) raises `EnvironmentError`. No environment variable is needed when OpenQP is installed via pip.

**`xyz_to_system_block()`** (`utils.py`): Converts element symbols in XYZ to atomic numbers for the OpenQP `system=` block. Only numeric atomic numbers appear in the final `.inp` file, not element symbols.

**`batch_from_wigner()` mutates `self.folder`** (`manager.py`): During trajectory creation the method temporarily sets `self.folder` to each trajectory subdirectory, then resets it to `base` at the end. If an exception occurs mid-loop, `self.folder` will be left pointing to a trajectory subdir.

**Absorption spectrum sample-index parsing** (`spectrum.py`): When reading log files, `absorption_spectrum()` infers sample index from the last `_`-separated token in the subfolder name (e.g. `mol_abs_0042` → index 41). Subfolders not matching this pattern fall back to sequential numbering.

**`TrajectoryAnalysis` state parsing** (`analysis.py`): Current state at each step is read from the comment line of `title.md.xyz` using the pattern `state (\d+)`. Surface-hop events are parsed separately from `title.sh.xyz` using `coord N state F to T CI`. Trajectories that have no xyz data (e.g. crashed on step 1) are silently skipped and counted as failed.

**Built-in templates loaded from disk** (`templates.py`): `BUILTIN_TEMPLATES` is populated at import time by `_load_builtin_templates()`, which reads all `*.json` files from the `templates/` subdirectory. `WIGNER_HESS_TEMPLATE` remains a Python dict (not a JSON file) because it is internal-only.

## Simulation Modes
- `qp`: Pure QM single-point (OpenQP only) — needs `oqp`
- `qmmm`: QM/MM (OpenQP + OpenMM via QMMM_MD) — needs `oqp`, `openmm`
- `namd`: Non-adiabatic MD (PyRAI2MD + OpenQP) — needs `oqp`, `PyRAI2MD`
- `namd_qmmm`: NAMD with QM/MM potential — needs `oqp`, `openmm`, `PyRAI2MD`

## Template System

Built-in templates (`qp_default`, `qmmm_default`, `namd_default`, `namd_qmmm_default`) live as JSON files in `templates/` and are loaded dynamically at import time. **All default configuration lives in these files — there are no hardcoded defaults in Python code.**

```python
# Instantiate from a built-in or custom JSON template
manager = SimulationManager.from_template("namd_default", folder="run", title="mol", xyz="mol.xyz")

# Inspect available templates
SimulationManager.list_templates()               # built-ins + any *.json in current dir
SimulationManager.list_templates("/path/dir")    # scan a specific directory

# Save current settings as a reusable template
manager.save_template("my_template.json", description="Project-specific settings")

# Create a template from scratch (no manager needed)
SimulationManager.create_template("foo.json", mode="namd", sections={...}, namd={...})
```

## Runtime Config Modification

After construction, sections can be modified without reconstructing the manager:

```python
manager.set_config("tdhf", "nstate", 8)          # modify OpenQP section, regenerates .inp
manager.set_namd_config("md", "step", 200)        # modify NAMD control file
manager.set_xyz("new.xyz")                        # swap geometry (rewrites .xyz and .inp for qp)
manager.set_pdb("new.pdb")                        # swap PDB for QM/MM modes
manager.show_config()                             # print all generated file contents
manager.show_config("sections")                   # return _sections dict
manager.show_config("openqp")                     # return raw .inp content string
```

## Complete Workflow (NAMD with absorption spectrum)

```python
from opensm import SimulationManager

# 1. Create manager
manager = SimulationManager(mode="namd", folder="run", title="azomethane", xyz="azo.xyz")

# 2. Wigner sampling (auto-runs OpenQP Hessian if no molden_file given)
manager.wigner_sampling(nsamples=1000, temp=300)

# 3. Absorption batch (geom_file auto-filled from wigner step)
manager.absorption_batch(partition="compute", cpus_per_task=16, modules=["anaconda3"])
manager.submit()         # submits absorption jobs
manager.save_project()   # save before closing Python

# --- new session after jobs finish ---
manager = SimulationManager.load_project("run/azomethane.project.json")
manager.absorption_spectrum(sigma=0.15)           # state-resolved plot (x_unit="eV" or "nm")

# 4. Select samples from energy window (all args auto-filled)
manager.select_wigner_window(e_min=3.5, e_max=4.5, target_state=2)
# or by wavelength:
# manager.select_wigner_window(w_min=275, w_max=354, target_state=2)

# 5. Create NAMD trajectories (auto uses selected geom/velo)
manager.batch_from_wigner(base_folder="namd_trajectories", partition="compute")
manager.submit()         # submits NAMD jobs

# --- new session after NAMD jobs finish ---
# 6. Analyze trajectory ensemble
ana = manager.analyze_trajectories()   # auto-detects base_folder from workflow
ana.summary()
ana.save()                             # writes population.dat, energy.dat, hop_statistics.dat, survival.dat
ana.plot_population("population.png")
ana.plot_energy("energy.png")
ana.plot_survival("survival.png")
```

## Key Files and Formats

### OpenQP Input (.inp)
INI-style sections: `[input]`, `[scf]`, `[tdhf]`, `[hess]`, `[properties]`, `[nac]`, `[qmmm]`, etc.
Generated by `sections_to_inp()` in `builders.py`. The `system=` block with atomic coordinates is injected from XYZ.

### PyRAI2MD Control File
Sections: `&CONTROL`, `&MOLECULE`, `&openqp`, `&MD`
Generated by `namd_dict_to_control()` in `builders.py`. OpenQP path auto-detected from installed `oqp` package via `importlib.util.find_spec("oqp")`.

### PyRAI2MD Trajectory Output (parsed by `analysis.py`)
Each trajectory folder contains:
- `title.md.energies` — columns: `time(au)  Epot(au)  Ekin(au)  Etot(au)  E_S1(au)  E_S2(au) ...`; header line contains "time"
- `title.md.xyz` — multi-frame XYZ; comment format: `title coord N state S`
- `title.sh.xyz` — hop-event frames only; comment format: `title coord N state F to T CI`
- `title.sh.energies` — same column format as `md.energies` but only at hop frames
- Time unit: atomic units (1 a.u. = 0.02418884 fs); energy unit: Hartree (1 Ha = 27.211386 eV)

### Wigner Output
- `wigner_geom.xyz`: multi-frame XYZ in Angstrom (natoms / [Angstrom] / element x y z)
- `wigner_velo.xyz`: velocities in Bohr/time_au (index [Bohr / time_au] / vx vy vz)

### Absorption Output
- `absorption_spectrum.dat`: columns = Energy(eV), Total, State_1->2, State_1->3, …
- `absorption_transitions.dat`: sample_idx, initial_state, final_state, energy_eV, osc_strength
- `absorption_spectrum.png`: state-resolved plot (eV or nm x-axis via `x_unit` parameter)

### Project File (.project.json)
Contains: mode, folder, title, sections, namd_params, xyz/pdb paths, generated inputs, workflow state

## OpenQP-Specific Details
- Hessian template (`WIGNER_HESS_TEMPLATE`): runtype=hess with MRSF-TDDFT, merged with manager's electronic-structure sections
- Absorption template: runtype=energy, strips dynamics-only sections (properties, nac, qmmm, hess)
- Transition table regex pattern: `istate -> fstate  E(eV)  col  col  col  col  osc_strength`
- Molden parser: reads `[FREQ]`, `[FR-COORD]`, `[FR-NORM-COORD]` sections
- Coordinates in Molden are in Bohr; output geometries in Angstrom

## Wigner Sampling Algorithm
For each normal mode i with frequency ω_i:
- σ_Q = sqrt(1/(2ω) · coth(ω/(2·k_B·T)))  (position width)
- σ_P = sqrt(ω/2 · coth(ω/(2·k_B·T)))     (momentum width)
- Sample Q_i ~ N(0, σ_Q) and P_i ~ N(0, σ_P)
- Transform back to Cartesian via mass-weighted normal modes
- `sigma` parameter in `absorption_spectrum()` is Gaussian broadening width in eV (typical: 0.10–0.30)

## Coding Conventions
- Python 3.8+ compatible
- numpy for numerical work, matplotlib for plotting
- pathlib.Path for all file operations
- deepcopy for template/section merging
- All batch scripts use absolute paths (Path.resolve())
- OpenQP imports are lazy (only when actually running calculations)
- Module functions are standalone-testable without SimulationManager
