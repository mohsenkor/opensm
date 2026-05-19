# OpenSM

Unified interface for multi-scale quantum chemistry simulations using OpenQP.
Orchestrates QM, QM/MM, NAMD, and NAMD+QM/MM workflows through a single
`SimulationManager` class.

**Method**: MRSF-TDDFT (multi-reference spin-flip TDDFT)  
**Backends**: OpenQP (QM), OpenMM (MM), PyRAI2MD (NAMD driver)

---

## Installation

### Requirements

- Python ≥ 3.8
- git
- A C/Fortran compiler (gfortran) for the fssh library

### Install

```bash
git clone git@github.com:mohsenkor/opensm.git
cd opensm
pip install .
pyrai2md update        # compile the fssh surface-hopping library (one-time step)
```

Or use the provided script which does both steps:

```bash
bash install.sh
```

### Check dependencies

```python
from opensm import SimulationManager
SimulationManager.check_dependencies()
```

---

## Simulation Modes

| Mode | Description | Requires |
|------|-------------|---------|
| `qp` | Pure QM single-point | OpenQP |
| `qmmm` | QM/MM | OpenQP + OpenMM |
| `namd` | Non-adiabatic MD | OpenQP + PyRAI2MD |
| `namd_qmmm` | NAMD with QM/MM | OpenQP + OpenMM + PyRAI2MD |

---

## Full NAMD Workflow

The typical workflow for a NAMD simulation with absorption spectrum
selection consists of six steps across two sessions.

### Session 1 — Setup, sampling, and job submission

```python
from opensm import SimulationManager

# 1. Create the manager
manager = SimulationManager(
    mode="namd",
    folder="formaldehyde_gas",
    title="test",
    xyz="ch2o.xyz",
)

# 2. Wigner sampling (runs an OpenQP Hessian automatically if no molden file)
wigner = manager.wigner_sampling(nsamples=600)

# 3. Create absorption single-point batch
abs_info = manager.absorption_batch(wigner['geom_file'], partition="ryzn")

# 4. Save project and submit
manager.save_project()
manager.submit()
```

### Session 2 — After absorption jobs finish

```python
from opensm import SimulationManager

# 5. Load project
manager = SimulationManager.load_project("formaldehyde_gas/test.project.json")

# 6. Tune electronic-structure and MD settings
manager.set_config("properties", "grad", 2)
manager.set_namd_config("md", "step", 10000)
manager.set_namd_config("md", "direct", 10000)
manager.set_namd_config("md", "root", 2)
manager.set_namd_config("md", "deco", 0.1)
manager.set_namd_config("molecule", "ci", 2)
manager.set_namd_config("molecule", "coupling", "1 2")

# 7. Inspect generated input files
manager.show_config()

# 8. Build the absorption spectrum
manager.absorption_spectrum(sigma=0.2)

# 9. Select Wigner samples from an energy window
manager.select_wigner_window(e_min=4.0, e_max=4.28, target_state=2)

# 10. Create NAMD trajectory folders and submit
manager.batch_from_wigner(
    base_folder="namd_trajectories",
    partition="ryzn,trd",
    cpus_per_task=16,
    time="24:00:00",
)
manager.submit("namd")
manager.save_project()
```

### Session 3 — After NAMD jobs finish

```python
from opensm import TrajectoryAnalysis

ana = TrajectoryAnalysis("namd_trajectories", title="test")
ana.load()
ana.summary()
ana.save()                              # writes population.dat, energy.dat, etc.
ana.plot_population("population.png")
ana.plot_energy("energy.png")
ana.plot_survival("survival.png")
```

Or via the manager (auto-detects the trajectory folder from the workflow):

```python
ana = manager.analyze_trajectories()
ana.summary()
ana.save()
```

---

## Configuration Reference

### `set_config(section, key, value)` — OpenQP input settings

```python
manager.set_config("tdhf", "nstate", 5)
manager.set_config("properties", "grad", 2)
manager.set_config("scf", "conv", 1e-8)
```

### `set_namd_config(section, key, value)` — PyRAI2MD control settings

```python
manager.set_namd_config("md", "step", 10000)   # number of MD steps
manager.set_namd_config("md", "size", 20.67)   # time step in a.u.
manager.set_namd_config("md", "root", 2)        # initial state
manager.set_namd_config("md", "deco", 0.1)      # decoherence correction
manager.set_namd_config("molecule", "ci", 2)    # number of states
manager.set_namd_config("molecule", "coupling", "1 2")  # NAC pairs
```

### Absorption spectrum options

```python
manager.absorption_spectrum(
    sigma=0.2,       # Gaussian broadening width in eV
    x_unit="eV",     # or "nm" for wavelength axis
)
```

### Wigner window selection — eV or nm

```python
# by energy (eV)
manager.select_wigner_window(e_min=4.0, e_max=4.28, target_state=2)

# by wavelength (nm)
manager.select_wigner_window(w_min=290, w_max=310, target_state=2)
```

---

## Template System

All default settings live in JSON files inside `opensm/templates/`.
Built-in templates are loaded automatically at import time — there are no
hardcoded defaults in Python code.

| Template | Mode | Description |
|----------|------|-------------|
| `qp_default` | `qp` | Pure QM single-point |
| `qmmm_default` | `qmmm` | QM/MM with OpenMM |
| `namd_default` | `namd` | Non-adiabatic MD |
| `namd_qmmm_default` | `namd_qmmm` | NAMD with QM/MM |

### Use a built-in template

```python
manager = SimulationManager.from_template(
    "namd_default", folder="run", title="mol", xyz="mol.xyz"
)
```

`SimulationManager(mode="namd", ...)` without a `template=` argument also
loads `namd_default` automatically.

### List available templates

```python
SimulationManager.list_templates()           # built-ins + any *.json in cwd
SimulationManager.list_templates("/my/dir")  # scan a specific directory
```

### Save current settings as a new template

```python
manager.save_template("my_mol_template.json", description="Formaldehyde NAMD")
```

This writes a JSON file you can reuse in any future project.

### Create a template from scratch (no manager needed)

```python
SimulationManager.create_template(
    "custom_namd.json",
    mode="namd",
    sections={
        "input":  {"functional": "bhhlyp", "basis": "6-31G*",
                   "method": "tdhf", "charge": 0},
        "scf":    {"type": "rohf", "multiplicity": 3,
                   "converger_type": "soscf", "conv": "1e-8"},
        "tdhf":   {"type": "mrsf", "nstate": 3},
        "dftgrid":{"rad_npts": 96, "ang_npts": 302, "pruned": ""},
        "properties": {"export": True, "nac": "nacme",
                       "back_door": True, "grad": 3},
        "nac": {},
    },
    namd={
        "control":  {"qc_ncpu": 16, "jobtype": "md",
                     "qm": "openqp", "abinit": "openqp"},
        "molecule": {"ci": 3, "spin": 0, "coupling": "1 2, 1 3, 2 3"},
        "openqp":   {"threads": 16, "use_hpc": -1, "align_mo": 1},
        "md": {"step": 5000, "size": 20.67, "root": 2,
               "sfhp": "fssh", "nactype": "dcm",
               "deco": "OFF", "substep": 20},
    },
    description="Custom 3-state NAMD for small molecules",
)
```

### Use a custom template

```python
# From a JSON file in the current directory or any path
manager = SimulationManager.from_template(
    "custom_namd.json", folder="run", title="mol", xyz="mol.xyz"
)

# Or pass the path explicitly
manager = SimulationManager.from_template(
    "/path/to/custom_namd.json", folder="run", title="mol", xyz="mol.xyz"
)
```

### Template file format

Templates are plain JSON files with this structure:

```json
{
  "description": "Human-readable description",
  "mode": "namd",
  "sections": {
    "input":   { "functional": "bhhlyp", "basis": "6-31G*", "method": "tdhf", "charge": 0 },
    "scf":     { "type": "rohf", "multiplicity": 3, "converger_type": "soscf", "conv": "1e-8" },
    "tdhf":    { "type": "mrsf", "nstate": 5 },
    "dftgrid": { "rad_npts": 96, "ang_npts": 302, "pruned": "" },
    "properties": { "export": true, "nac": "nacme", "back_door": true, "grad": 5 },
    "nac": {}
  },
  "namd": {
    "control":  { "qc_ncpu": 16, "jobtype": "md", "qm": "openqp", "abinit": "openqp" },
    "molecule": { "ci": 5, "spin": 0, "coupling": "1 2, 1 3, 2 3" },
    "openqp":   { "threads": 16, "use_hpc": -1, "align_mo": 1 },
    "md": {
      "step": 5000, "size": 20.67, "root": 2,
      "sfhp": "fssh", "nactype": "dcm", "deco": "OFF",
      "substep": 20, "initcond": 0, "randvelo": 0
    }
  }
}
```

The `namd` block is only required for `namd` and `namd_qmmm` modes.
You can omit any key — missing keys fall back to the built-in template defaults.

---

## Output Files

| File | Description |
|------|-------------|
| `wigner_geom.xyz` | Wigner-sampled geometries (Angstrom) |
| `wigner_velo.xyz` | Wigner-sampled velocities (Bohr/time_au) |
| `absorption_spectrum.dat` | Energy(eV), Total, State_1→2, … |
| `absorption_transitions.dat` | Per-sample transition energies and oscillator strengths |
| `absorption_spectrum.png` | State-resolved absorption spectrum plot |
| `*.project.json` | Serialized project state for cross-session workflows |
| `population.dat` | Ensemble state populations vs time |
| `energy.dat` | Ensemble average energies vs time |
| `hop_statistics.dat` | Surface-hopping event summary |
| `survival.dat` | Survival probability of the initial state |

---

## License

MIT
