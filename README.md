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

### Set environment variable (required for NAMD modes)

```bash
export OPENQP_ROOT=/path/to/openqp
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

Built-in templates (`qp_default`, `qmmm_default`, `namd_default`,
`namd_qmmm_default`) provide reasonable defaults for each mode.

```python
# Instantiate from a built-in template
manager = SimulationManager.from_template(
    "namd_default", folder="run", title="mol", xyz="mol.xyz"
)

# List available templates
SimulationManager.list_templates()

# Save current settings as a reusable template
manager.save_template("my_template.json")
```

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
