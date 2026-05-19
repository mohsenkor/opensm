"""
SimulationManager — Unified Python interface for multi-scale simulations
in OpenQP.

Supports four modes:
    qp        : Pure QM single-point (OpenQP)
    qmmm      : QM/MM (OpenQP + OpenMM via QMMM_MD)
    namd      : Non-adiabatic MD (PyRAI2MD + OpenQP)
    namd_qmmm : NAMD with QM/MM potential (PyRAI2MD + OpenQP + OpenMM)

Minimal usage
-------------
    from opensm import SimulationManager

    manager = SimulationManager.from_template(
        "namd_default",
        folder="Nac_run",
        title="test1-1",
        xyz="test1-1.xyz",
    )
    manager.show_config()
    manager.run()
"""

import json
import os
from copy import deepcopy
from pathlib import Path

from .constants import ELEMENT_TO_Z
from .utils import (
    resolve_content, check_all_dependencies, check_package,
    check_mode_dependencies, MODE_DEPENDENCIES,
)
from .builders import sections_to_inp, namd_dict_to_control
from .templates import BUILTIN_TEMPLATES
from .engines import OpenQPEngine, PyRAI2MDEngine
from . import wigner as _wigner
from . import spectrum as _spectrum
from . import batch as _batch
from . import project as _project
from . import analysis as _analysis


class SimulationManager:
    """
    Unified interface for QM, QM/MM, NAMD, and NAMD-QM/MM simulations.

    Parameters
    ----------
    mode : str
        One of 'qp', 'qmmm', 'namd', 'namd_qmmm'.
    folder : str
        Working directory.
    title : str
        Project name (used for all filenames).
    xyz : str or path
        XYZ coordinates.
    pdb : str or path, optional
        PDB file for QM/MM modes.
    qm_atoms : str, optional
        QM atom selection for QM/MM (e.g. '0-2').
    sections : dict, optional
        OpenQP section overrides.
    namd : dict, optional
        PyRAI2MD section overrides (control, molecule, openqp, md).
    template : str, optional
        JSON template path or built-in name.
    openqp_content, xyz_content, input_content : str, optional
        Pre-built inputs (bypass auto-generation).
    configs : dict, optional
        Engine-level config overrides.
    """

    VALID_MODES = {"qp", "qmmm", "namd", "namd_qmmm"}

    def __init__(
        self,
        mode=None,
        folder=".",
        title="test1-1",
        xyz=None,
        pdb=None,
        qm_atoms=None,
        sections=None,
        namd=None,
        template=None,
        openqp_content=None,
        xyz_content=None,
        input_content=None,
        configs=None,
    ):
        # ── Load template ──
        tpl = {}
        if template is not None:
            tpl = self._load_template(template)

        self.mode = mode or tpl.get("mode")
        if self.mode not in self.VALID_MODES:
            raise ValueError(
                f"Unsupported mode '{self.mode}'. Valid: {sorted(self.VALID_MODES)}"
            )

        self.folder = folder
        self.title = title
        self.configs = configs or {}
        self._template_source = template
        self._workflow = {}

        self._xyz = xyz
        self._pdb = pdb
        self._qm_atoms = qm_atoms

        # ── Build OpenQP sections ──
        self._sections = deepcopy(tpl.get("sections", {}))
        if sections:
            for sec_name, sec_params in sections.items():
                if sec_name in self._sections:
                    self._sections[sec_name].update(sec_params)
                else:
                    self._sections[sec_name] = deepcopy(sec_params)
        if not self._sections:
            self._sections = self._default_sections()

        # inject molecule-specific values
        if self._pdb is not None and "qmmm" in self._sections:
            pdb_path = Path(self._pdb)
            self._sections["qmmm"]["pdb_file"] = (
                pdb_path.name if pdb_path.is_file() else f"{self.title}.pdb"
            )
        if self._qm_atoms is not None and "qmmm" in self._sections:
            self._sections["qmmm"]["qm_atoms"] = self._qm_atoms

        # ── Build NAMD params ──
        if self.mode in {"namd", "namd_qmmm"}:
            self._namd_params = self._default_namd()
            for sec_name, sec_params in tpl.get("namd", {}).items():
                if sec_name in self._namd_params and isinstance(sec_params, dict):
                    self._namd_params[sec_name].update(deepcopy(sec_params))
                else:
                    self._namd_params[sec_name] = deepcopy(sec_params)
            if namd:
                for sec_name, sec_params in namd.items():
                    if sec_name in self._namd_params and isinstance(sec_params, dict):
                        self._namd_params[sec_name].update(sec_params)
                    else:
                        self._namd_params[sec_name] = (
                            deepcopy(sec_params) if isinstance(sec_params, dict) else sec_params
                        )
        else:
            self._namd_params = {}

        # ── Generate inputs ──
        self._generated = {}
        self._generate_inputs(openqp_content, xyz_content, input_content)

        # ── Engines ──
        self.openqp = None
        self.pyrai2md = None
        self._setup_engines()
        self._setup_files()

    # ================================================================
    #  Default sections / NAMD params
    # ================================================================
    def _default_sections(self):
        if self.mode == "qp":
            return {
                "input": {"functional": "bhhlyp", "basis": "6-31g*",
                          "method": "hf", "runtype": "energy", "charge": 0},
                "scf": {"type": "rhf", "maxit": 100, "multiplicity": 1},
                "guess": {"type": "huckel"},
                "dftgrid": {"rad_type": "becke"},
            }
        if self.mode == "qmmm":
            return {
                "input": {"functional": "bhhlyp", "basis": "6-31g*",
                          "method": "tdhf", "qmmm_flag": True},
                "scf": {"type": "rohf", "maxit": 200, "multiplicity": 3,
                        "converger_type": "soscf", "conv": "1e-8"},
                "guess": {"type": "huckel"},
                "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
                "tdhf": {"type": "mrsf", "nstate": 6, "maxit": 100,
                         "maxit_zv": 200, "nvdav": 100, "zvconv": "1.0e-10"},
                "qmmm": {"pdb_file": "", "forcefield_files": "amber14-all.xml",
                         "qm_atoms": "", "cutoff": "PME",
                         "embedding": "electrostatic", "n_steps": 1000,
                         "timestep": 1.0, "temperature": 300.0},
            }
        if self.mode == "namd":
            return {
                "input": {"functional": "bhhlyp", "charge": 0,
                          "method": "tdhf", "basis": "6-31G*"},
                "guess": {"type": "huckel"},
                "scf": {"type": "rohf", "maxit": 200, "multiplicity": 3,
                        "converger_type": "soscf", "conv": "1e-8"},
                "tdhf": {"type": "mrsf", "nstate": 6, "maxit": 100,
                         "maxit_zv": 200, "nvdav": 100, "zvconv": "1.0e-10"},
                "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
                "properties": {"export": True, "nac": "nacme",
                               "back_door": True, "grad": 5},
                "nac": {},
            }
        if self.mode == "namd_qmmm":
            return {
                "input": {"functional": "bhhlyp", "charge": 0,
                          "method": "tdhf", "basis": "6-31G*",
                          "qmmm_flag": True},
                "guess": {"type": "huckel"},
                "scf": {"type": "rohf", "maxit": 200, "multiplicity": 3,
                        "converger_type": "soscf", "conv": "1e-8"},
                "tdhf": {"type": "mrsf", "nstate": 6, "maxit": 100,
                         "maxit_zv": 200, "nvdav": 100, "zvconv": "1.0e-10"},
                "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
                "properties": {"export": True, "nac": "nacme",
                               "back_door": True, "grad": 5},
                "nac": {},
                "qmmm": {"pdb_file": "", "forcefield_files": "amber14-all.xml",
                         "qm_atoms": "", "cutoff": "PME",
                         "embedding": "electrostatic", "n_steps": 1000,
                         "timestep": 1.0, "temperature": 300.0},
            }
        return {}

    def _default_namd(self):
        return {
            "control": {"qc_ncpu": 16, "gl_seed": 1, "jobtype": "md",
                        "qm": "openqp", "abinit": "openqp"},
            "molecule": {"ci": 5, "spin": 0,
                         "coupling": "1 2, 1 3, 1 4, 1 5, 2 3, 2 4, 2 5, 3 4, 3 5, 4 5"},
            "openqp": {"openqp": "", "threads": 40, "use_hpc": -1, "align_mo": 1},
            "md": {"initcond": 0, "temp": 300, "randvelo": 0,
                   "step": 80, "direct": 80, "checkpoint": 1, "buffer": 0,
                   "size": 20.67, "root": 5, "activestate": 1,
                   "sfhp": "fssh", "nactype": "dcm", "phasecheck": 0,
                   "adjust": 0, "reflect": 1, "deco": "OFF",
                   "substep": 20, "verbose": 1, "thermo": 0, "restart": ""},
        }

    # ================================================================
    #  Dependency checking
    # ================================================================
    def _check_dependencies(self):
        check_mode_dependencies(self.mode)

    @staticmethod
    def check_dependencies():
        """Print installation status of all required packages."""
        results = check_all_dependencies()
        print("Package dependencies:")
        print("-" * 60)
        for pkg, info in results.items():
            status = f"\u2713 {info['info']}" if info["available"] else "\u2717 not installed"
            print(f"  {info['label']}    {status}")
        print()
        print("Mode requirements:")
        print("-" * 60)
        for mode, pkgs in MODE_DEPENDENCIES.items():
            all_ok = all(results[p]["available"] for p in pkgs)
            mark = "\u2713" if all_ok else "\u2717"
            print(f"  {mark} {mode:12s}  needs: {', '.join(pkgs)}")

    # ================================================================
    #  Template system
    # ================================================================
    @staticmethod
    def _load_template(source):
        if isinstance(source, str) and source in BUILTIN_TEMPLATES:
            return deepcopy(BUILTIN_TEMPLATES[source])
        path = Path(source)
        if not path.is_file():
            raise FileNotFoundError(
                f"Template not found: '{source}'. "
                f"Built-in: {list(BUILTIN_TEMPLATES.keys())}"
            )
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_template(self, filepath, description=None):
        """Save current settings as a reusable JSON template."""
        tpl = {
            "description": description or f"SimulationManager template ({self.mode})",
            "mode": self.mode,
            "sections": deepcopy(self._sections),
        }
        if self._namd_params:
            tpl["namd"] = deepcopy(self._namd_params)
        if "qmmm" in tpl["sections"]:
            tpl["sections"]["qmmm"]["pdb_file"] = ""
            tpl["sections"]["qmmm"]["qm_atoms"] = ""
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(tpl, f, indent=2, ensure_ascii=False)
        print(f"Template saved: {filepath}")

    @classmethod
    def from_template(cls, template, **overrides):
        """Create a SimulationManager from a template with per-job overrides."""
        return cls(template=template, **overrides)

    @staticmethod
    def list_templates(directory=None):
        """List built-in and user templates."""
        from .templates import TEMPLATES_DIR
        print(f"Built-in templates  ({TEMPLATES_DIR}):")
        print("-" * 60)
        for name, tpl in BUILTIN_TEMPLATES.items():
            print(f"  {name:24s} {tpl.get('description', '')}")
        search_dir = Path(directory) if directory else Path(".")
        user_templates = []
        for f in sorted(search_dir.glob("*.json")):
            try:
                with open(f) as fh:
                    data = json.load(fh)
                if "mode" in data and "sections" in data:
                    user_templates.append((f.name, data.get("description", "")))
            except (json.JSONDecodeError, KeyError):
                continue
        if user_templates:
            print(f"\nUser templates in {search_dir}:")
            print("-" * 60)
            for name, desc in user_templates:
                print(f"  {name:24s} {desc}")

    @staticmethod
    def create_template(filepath, mode, description=None, sections=None, namd=None):
        """Create a template JSON file from scratch."""
        tpl = {
            "description": description or f"Custom template ({mode})",
            "mode": mode,
            "sections": sections or {},
        }
        if namd:
            tpl["namd"] = namd
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(tpl, f, indent=2, ensure_ascii=False)
        print(f"Template created: {filepath}")

    # ================================================================
    #  Project persistence
    # ================================================================
    def save_project(self, filepath=None):
        """Save the full manager state to a JSON project file."""
        return _project.save_project(self, filepath)

    @classmethod
    def load_project(cls, filepath):
        """Reload a SimulationManager from a saved project file."""
        state = _project.load_project_state(filepath)

        obj = object.__new__(cls)
        obj.mode = state["mode"]
        obj.folder = state["folder"]
        obj.title = state["title"]
        obj._template_source = state.get("template_source")
        obj._sections = state.get("sections", {})
        obj._namd_params = state.get("namd_params", {})
        obj._xyz = state.get("xyz")
        obj._pdb = state.get("pdb")
        obj._qm_atoms = state.get("qm_atoms")
        obj._generated = state.get("generated", {})
        obj.configs = state.get("configs", {})
        obj._workflow = state.get("workflow", {})

        obj.openqp = None
        obj.pyrai2md = None
        obj._setup_engines()

        print(f"Project loaded: {filepath}")
        print(f"  Mode: {obj.mode} | Folder: {obj.folder} | Title: {obj.title}")

        wf = obj._workflow
        if wf:
            if "wigner" in wf:
                w = wf["wigner"]
                print(f"  Wigner: {w.get('nsamples', '?')} samples "
                      f"\u2192 {w.get('geom_file', '?')}")
            if "absorption_batch" in wf:
                ab = wf["absorption_batch"]
                print(f"  Absorption batch: {ab.get('n_samples', '?')} samples "
                      f"in {ab.get('base_folder', '?')}")
            if "absorption_spectrum" in wf:
                sp = wf["absorption_spectrum"]
                print(f"  Spectrum: {sp.get('n_converged', '?')} converged "
                      f"\u2192 {sp.get('spectrum_file', '?')}")
            if "selected" in wf:
                sel = wf["selected"]
                print(f"  Selection: {sel.get('n_selected', '?')} samples "
                      f"in [{sel.get('e_min', '?')}, {sel.get('e_max', '?')}] eV")

        return obj

    def _save_workflow(self, key, data):
        """Store workflow step results and auto-save if project file exists."""
        self._workflow[key] = data
        proj_path = Path(self.folder) / f"{self.title}.project.json"
        if proj_path.is_file():
            self.save_project(str(proj_path))

    def show_workflow(self):
        """Print a summary of completed workflow steps."""
        _project.show_workflow(self._workflow, self.title, self.mode, self.folder)

    def get_workflow(self, step=None):
        """Retrieve saved workflow results."""
        return _project.get_workflow(self._workflow, step)

    # ================================================================
    #  Input generation
    # ================================================================
    def _generate_inputs(self, openqp_content, xyz_content, input_content):
        if openqp_content is not None:
            text, _ = resolve_content(openqp_content)
            self._generated["openqp"] = text
        else:
            xyz_for_system = self._xyz if self.mode == "qp" else None
            self._generated["openqp"] = sections_to_inp(
                self._sections, xyz=xyz_for_system,
            )

        if xyz_content is not None:
            text, _ = resolve_content(xyz_content)
            self._generated["xyz"] = text
        elif self._xyz is not None:
            text, _ = resolve_content(self._xyz)
            self._generated["xyz"] = text

        if self._pdb is not None:
            text, _ = resolve_content(self._pdb)
            self._generated["pdb"] = text

        if input_content is not None:
            text, _ = resolve_content(input_content)
            self._generated["namd"] = text
        elif self.mode in {"namd", "namd_qmmm"} and self._namd_params:
            self._generated["namd"] = namd_dict_to_control(
                self._namd_params, self.title,
            )

    # ================================================================
    #  Config inspection & modification
    # ================================================================
    def show_config(self, name=None):
        """Inspect generated inputs."""
        if name == "sections":
            return deepcopy(self._sections)
        if name == "namd_sections":
            return deepcopy(self._namd_params)
        if name is not None:
            return self._generated.get(name)

        sep = "=" * 60
        for key, content in self._generated.items():
            print(f"\n{sep}\n  {key.upper()} INPUT\n{sep}")
            print(content)
        print(f"\n{sep}")
        print(f"  MODE: {self.mode}  |  FOLDER: {self.folder}  |  TITLE: {self.title}")
        if self._template_source:
            print(f"  TEMPLATE: {self._template_source}")
        print(sep)

    def set_config(self, section, key, value):
        """Modify OpenQP sections and regenerate."""
        if section not in self._sections:
            self._sections[section] = {}
        self._sections[section][key] = value
        xyz_for_system = self._xyz if self.mode == "qp" else None
        self._generated["openqp"] = sections_to_inp(
            self._sections, xyz=xyz_for_system,
        )
        print(f"Config updated: [{section}] {key} = {value}")

    def set_namd_config(self, section, key, value):
        """Modify PyRAI2MD control sections and regenerate."""
        if section not in self._namd_params:
            self._namd_params[section] = {}
        self._namd_params[section][key] = value
        self._generated["namd"] = namd_dict_to_control(
            self._namd_params, self.title,
        )
        print(f"NAMD config updated: &{section.upper()} {key} = {value}")

    def set_xyz(self, xyz):
        text, _ = resolve_content(xyz)
        self._generated["xyz"] = text
        self._xyz = xyz
        if self.mode == "qp":
            self._generated["openqp"] = sections_to_inp(
                self._sections, xyz=xyz,
            )
        print("XYZ content updated.")

    def set_pdb(self, pdb):
        text, _ = resolve_content(pdb)
        self._generated["pdb"] = text
        self._pdb = pdb
        print("PDB content updated.")

    # ================================================================
    #  Engine setup
    # ================================================================
    def _setup_engines(self):
        if self.mode in {"qp", "qmmm", "namd", "namd_qmmm"}:
            self.openqp = OpenQPEngine(self.configs.get("openqp"))
        if self.mode in {"namd", "namd_qmmm"}:
            self.pyrai2md = PyRAI2MDEngine(
                self.mode, self.configs.get("pyrai2md"),
            )

    def _setup_files(self):
        """Write all generated inputs to disk."""
        openqp_content = self._generated.get("openqp", "")
        xyz_content = self._generated.get("xyz", "")
        pdb_content = self._generated.get("pdb", "")

        folder = Path(self.folder)
        folder.mkdir(parents=True, exist_ok=True)

        if pdb_content and self.mode in {"qmmm", "namd_qmmm"}:
            pdb_filename = self._sections.get("qmmm", {}).get(
                "pdb_file", f"{self.title}.pdb",
            )
            (folder / pdb_filename).write_text(
                pdb_content.strip() + "\n", encoding="utf-8",
            )

        if self.mode == "qp" and openqp_content:
            self.openqp.setup_input(
                folder=self.folder, openqp_content=openqp_content,
                title=self.title,
            )
        elif self.mode == "qmmm" and openqp_content:
            if xyz_content:
                self.openqp.setup_input_with_xyz(
                    folder=self.folder, openqp_content=openqp_content,
                    xyz_content=xyz_content, title=self.title,
                )
            else:
                self.openqp.setup_input(
                    folder=self.folder, openqp_content=openqp_content,
                    title=self.title,
                )

        if self.mode in {"namd", "namd_qmmm"}:
            namd_content = self._generated.get("namd", "")
            self.pyrai2md.setup_files(
                folder=self.folder, input_content=namd_content,
                openqp_content=openqp_content, xyz_content=xyz_content,
                title=self.title,
            )

    # ================================================================
    #  Prepare & Run
    # ================================================================
    def prepare(self):
        if self.openqp is not None:
            self.openqp.prepare()
        if self.pyrai2md is not None:
            self.pyrai2md.prepare()

    def run(self):
        """Rewrite files from current config, then execute."""
        self._check_dependencies()
        self._setup_files()
        self.prepare()
        if self.mode == "qp":
            return self.openqp.run()
        elif self.mode == "qmmm":
            return self.openqp.run_qmmm()
        elif self.mode in {"namd", "namd_qmmm"}:
            self.pyrai2md.run(folder=self.folder, title=self.title)

    # ================================================================
    #  Wigner sampling  (delegates to opensm.wigner)
    # ================================================================
    def wigner_sampling(
        self, molden_file=None, nsamples=1000, temp=298.15, seed=1,
        scale=1.0, refxyz=None, skip_first=0, skip_last=0,
        output_dir=None, hess_sections=None, hess_state=1, hess_dir=None,
    ):
        """Generate Wigner-sampled initial conditions.
        See opensm.wigner.wigner_sampling() for full parameter docs."""
        result = _wigner.wigner_sampling(
            molden_file=molden_file,
            nsamples=nsamples, temp=temp, seed=seed, scale=scale,
            refxyz=refxyz if refxyz is not None else self._xyz,
            skip_first=skip_first, skip_last=skip_last,
            output_dir=output_dir or self.folder,
            sections=self._sections, xyz=self._xyz, title=self.title,
            hess_sections=hess_sections, hess_state=hess_state,
            hess_dir=hess_dir,
        )
        self._save_workflow("wigner", result)
        return result

    # ================================================================
    #  Absorption spectrum  (delegates to opensm.spectrum)
    # ================================================================
    def absorption_batch(
        self, geom_file=None, n_samples=None, start_sample=0,
        base_folder=None, spectrum_sections=None,
        scheduler="slurm", submit=False, job_name=None,
        partition=None, nodes=1, ntasks=1, cpus_per_task=16,
        mem=None, time="02:00:00", account=None, modules=None,
        pre_commands=None, post_commands=None, extra_headers=None,
    ):
        """Create OpenQP single-point folders for each Wigner sample.

        ``geom_file`` is optional if ``wigner_sampling()`` has been run
        (or loaded from a project file).

        See opensm.spectrum.absorption_batch() for full parameter docs."""
        # ── Auto-fill from workflow ──
        if geom_file is None:
            wf = self.get_workflow("wigner")
            if wf is None:
                raise ValueError("geom_file not provided and no wigner_sampling in workflow")
            geom_file = wf["geom_file"]

        base = base_folder or str(Path(self.folder) / "absorption")
        result = _spectrum.absorption_batch(
            geom_file=geom_file, sections=self._sections,
            title=self.title, base_folder=base,
            n_samples=n_samples, start_sample=start_sample,
            spectrum_sections=spectrum_sections,
            job_script_func=_batch.generate_job_script,
            scheduler=scheduler, submit=submit, job_name=job_name,
            partition=partition, nodes=nodes, ntasks=ntasks,
            cpus_per_task=cpus_per_task, mem=mem, time=time,
            account=account, modules=modules,
            pre_commands=pre_commands, post_commands=post_commands,
            extra_headers=extra_headers,
            master_script_func=_batch.generate_master_script,
            submit_func=_batch.submit_jobs if submit else None,
        )
        self._save_workflow("absorption_batch", {
            "base_folder": result["base_folder"],
            "n_samples": result["n_samples"],
            "geom_file": str(geom_file),
        })
        return result

    def submit(self, workflow_step=None, script_path=None):
        """
        Submit batch jobs by running a submit_all.sh script.

        If ``script_path`` is given, runs that script directly.
        Otherwise auto-detects from the latest workflow step that
        produced a submit_all.sh.

        Parameters
        ----------
        workflow_step : str, optional
            Which workflow step to submit: 'absorption_batch' or
            'batch' / 'batch_from_wigner'.  If None, auto-detects.
        script_path : str or Path, optional
            Explicit path to submit_all.sh (overrides auto-detection).

        Example
        -------
            manager.absorption_batch(partition="compute")
            manager.submit()  # runs absorption/submit_all.sh
        """
        import subprocess

        if script_path is not None:
            submit_script = Path(script_path)
        else:
            # Auto-detect from workflow
            candidates = []
            if workflow_step:
                steps = [workflow_step]
            else:
                steps = ["absorption_batch", "selected", "wigner"]

            for step in steps:
                wf = self.get_workflow(step)
                if wf and "base_folder" in wf:
                    candidate = Path(wf["base_folder"]) / "submit_all.sh"
                    if candidate.is_file():
                        candidates.append(candidate)

            # Also check the main folder
            main_submit = Path(self.folder) / "submit_all.sh"
            if main_submit.is_file():
                candidates.append(main_submit)

            if not candidates:
                raise FileNotFoundError(
                    "No submit_all.sh found. Run absorption_batch() or "
                    "batch_from_wigner() first, or pass script_path= explicitly."
                )
            submit_script = candidates[0]

        if not submit_script.is_file():
            raise FileNotFoundError(f"Submit script not found: {submit_script}")

        print(f"Running: bash {submit_script}")
        result = subprocess.run(
            ["bash", str(submit_script)],
            capture_output=True, text=True,
        )
        if result.stdout:
            print(result.stdout)
        if result.returncode != 0:
            print(f"Error (exit code {result.returncode}):")
            print(result.stderr)
        else:
            print("All jobs submitted successfully.")

    def absorption_spectrum(
        self, base_folder=None, initial_state=1,
        e_min=None, e_max=None, npts=5000, sigma=0.15,
        output_dir=None, plot=True, figsize=(9, 5.5), x_unit="eV",
    ):
        """Read results and build state-resolved absorption spectrum.
        See opensm.spectrum.absorption_spectrum() for full parameter docs."""
        base = base_folder or str(Path(self.folder) / "absorption")
        result = _spectrum.absorption_spectrum(
            base_folder=base, initial_state=initial_state,
            e_min=e_min, e_max=e_max, npts=npts, sigma=sigma,
            output_dir=output_dir, plot=plot, figsize=figsize, x_unit=x_unit,
        )
        self._save_workflow("absorption_spectrum", {
            "all_transitions": {str(k): v for k, v in result["all_transitions"].items()},
            "n_converged": result["n_converged"],
            "n_failed": result["n_failed"],
            "spectrum_file": result["spectrum_file"],
            "transitions_file": result.get("transitions_file"),
            "plot_file": result["plot_file"],
            "e_min": float(result["e_grid"][0]),
            "e_max": float(result["e_grid"][-1]),
            "sigma": sigma,
        })
        return result

    def select_wigner_window(
        self, geom_file=None, velo_file=None, all_transitions=None,
        e_min=None, e_max=None, w_min=None, w_max=None,
        target_state=None, output_dir=None,
    ):
        """Filter Wigner samples by an energy window.

        All file arguments are optional if previous workflow steps have
        been run (or loaded from a project file).  Defaults:
            geom_file       ← wigner_sampling result
            velo_file       ← wigner_sampling result
            all_transitions ← absorption_spectrum result

        Parameters
        ----------
        geom_file, velo_file : str or Path, optional
        all_transitions : dict, optional
        e_min, e_max : float, optional
            Energy window in eV.
        w_min, w_max : float, optional
            Wavelength window in nm (alternative to e_min/e_max).
        target_state : int, optional
        output_dir : str or Path, optional

        See opensm.spectrum.select_wigner_window() for full docs.
        """
        # ── Auto-fill from workflow ──
        if geom_file is None:
            wf = self.get_workflow("wigner")
            if wf is None:
                raise ValueError("geom_file not provided and no wigner_sampling in workflow")
            geom_file = wf["geom_file"]
        if velo_file is None:
            wf = self.get_workflow("wigner")
            if wf is None:
                raise ValueError("velo_file not provided and no wigner_sampling in workflow")
            velo_file = wf["velo_file"]
        if all_transitions is None:
            wf = self.get_workflow("absorption_spectrum")
            if wf is None:
                raise ValueError("all_transitions not provided and no absorption_spectrum in workflow")
            all_transitions = wf["all_transitions"]
        if e_min is None and e_max is None and w_min is None and w_max is None:
            raise ValueError(
                "Provide an energy window: e_min/e_max (eV) or w_min/w_max (nm)."
            )

        result = _spectrum.select_wigner_window(
            geom_file=geom_file, velo_file=velo_file,
            all_transitions=all_transitions,
            e_min=e_min, e_max=e_max, w_min=w_min, w_max=w_max,
            target_state=target_state,
            output_dir=output_dir or self.folder,
        )
        self._save_workflow("selected", {
            "geom_file": result["geom_file"],
            "velo_file": result["velo_file"],
            "n_selected": result["n_selected"],
            "n_total": result["n_total"],
            "selected_indices": result["selected_indices"],
            "e_min": e_min, "e_max": e_max,
            "target_state": target_state,
        })
        return result

    # ================================================================
    #  Batch NAMD from Wigner samples
    # ================================================================
    def batch_from_wigner(
        self, geom_file=None, velo_file=None, n_traj=None, base_folder=None,
        start_sample=0, scheduler="slurm", submit=False,
        job_name=None, partition=None, nodes=1, ntasks=1,
        cpus_per_task=16, mem=None, time="24:00:00", account=None,
        modules=None, pre_commands=None, post_commands=None,
        python_cmd="python", extra_headers=None,
    ):
        """Create NAMD trajectory folders from Wigner-sampled initial conditions.

        File arguments are optional if previous workflow steps have been
        run.  Defaults:
            geom_file ← select_wigner_window result (or wigner_sampling)
            velo_file ← select_wigner_window result (or wigner_sampling)
        """
        if self.mode not in {"namd", "namd_qmmm"}:
            raise ValueError(f"batch_from_wigner() is for NAMD modes, not '{self.mode}'")

        # ── Auto-fill from workflow: prefer selected, fallback to wigner ──
        if geom_file is None or velo_file is None:
            sel = self.get_workflow("selected")
            wig = self.get_workflow("wigner")
            source = sel or wig
            if source is None:
                raise ValueError(
                    "geom_file/velo_file not provided and no wigner or selection in workflow"
                )
            if geom_file is None:
                geom_file = source["geom_file"]
            if velo_file is None:
                velo_file = source["velo_file"]

        from .parsers import read_xyz_frames, read_velo_frames

        geom_frames = read_xyz_frames(geom_file)
        velo_frames = read_velo_frames(velo_file)
        n_available = min(len(geom_frames), len(velo_frames))

        if n_traj is None:
            n_traj = n_available - start_sample
        if start_sample + n_traj > n_available:
            raise ValueError(
                f"Requested {n_traj} trajectories starting at {start_sample}, "
                f"but only {n_available} available."
            )

        print(f"Creating {n_traj} NAMD trajectories from Wigner samples "
              f"[{start_sample}:{start_sample + n_traj}]")

        base = Path(base_folder) if base_folder else Path(self.folder)
        base.mkdir(parents=True, exist_ok=True)

        job_prefix = job_name or self.title
        modules = modules or []
        pre_commands = pre_commands or []
        post_commands = post_commands or []
        extra_headers = extra_headers or []

        folders = []
        scripts = []

        for i in range(n_traj):
            sample_idx = start_sample + i
            traj_idx = i + 1
            traj_name = f"{self.title}_traj_{traj_idx:04d}"
            traj_folder = base / traj_name

            atoms, coords_ang = geom_frames[sample_idx]
            velo = velo_frames[sample_idx]
            natoms = len(atoms)

            xyz_lines = [f"{natoms}", ""]
            for a in range(natoms):
                xyz_lines.append(
                    f"{atoms[a]:<2s}  {coords_ang[a][0]:20.12f}  "
                    f"{coords_ang[a][1]:20.12f}  {coords_ang[a][2]:20.12f}"
                )
            xyz_content = "\n".join(xyz_lines)

            velo_lines = []
            for a in range(natoms):
                velo_lines.append(
                    f"  {velo[a][0]:24.16e}  {velo[a][1]:24.16e}  "
                    f"{velo[a][2]:24.16e}"
                )
            velo_content = "\n".join(velo_lines)

            self.set_namd_config("control", "gl_seed", sample_idx + 1)
#            self.set_namd_config("md", "initcond", 1)
            self.folder = str(traj_folder)
            self._generated["xyz"] = xyz_content
            self._setup_engines()
            self._setup_files()

            traj_folder = Path(traj_folder)
            (traj_folder / f"{self.title}.velo").write_text(
                velo_content + "\n", encoding="utf-8",
            )

            run_script = _batch.generate_run_script(self.title)
            (traj_folder / "run.py").write_text(run_script, encoding="utf-8")

            script_content = _batch.generate_job_script(
                traj_folder=traj_folder, traj_name=traj_name,
                traj_idx=traj_idx, job_prefix=job_prefix,
                scheduler=scheduler, partition=partition,
                nodes=nodes, ntasks=ntasks, cpus_per_task=cpus_per_task,
                mem=mem, time=time, account=account, modules=modules,
                pre_commands=pre_commands, post_commands=post_commands,
                python_cmd=python_cmd, extra_headers=extra_headers,
            )
            script_path = traj_folder / "submit.sh"
            script_path.write_text(script_content, encoding="utf-8")
            os.chmod(script_path, 0o755)

            folders.append(traj_folder)
            scripts.append(script_path)

        self.folder = str(base)

        print(f"\nCreated {n_traj} trajectory folders in: {base}")
        print(f"  Folders: {self.title}_traj_0001 \u2192 {self.title}_traj_{n_traj:04d}")
        print(f"  Each folder contains: {self.title} (control), "
              f"{self.title}.openqp, {self.title}.xyz, {self.title}.velo")

        if submit:
            _batch.submit_jobs(scripts, scheduler)

        master = _batch.generate_master_script(base, scripts, scheduler)
        master_path = base / "submit_all.sh"
        master_path.write_text(master, encoding="utf-8")
        os.chmod(master_path, 0o755)
        print(f"  Master: {master_path}")

        return folders

    # ================================================================
    #  Trajectory analysis  (delegates to opensm.analysis)
    # ================================================================
    def analyze_trajectories(self, base_folder=None, title=None, load=True):
        """
        Create a TrajectoryAnalysis for an ensemble of NAMD runs.

        Parameters
        ----------
        base_folder : str or Path, optional
            Directory containing trajectory subfolders.
            Defaults to the folder set by batch_from_wigner() in the
            workflow, or self.folder if unavailable.
        title : str, optional
            Molecule title used in output filenames.
            Defaults to self.title.
        load : bool
            Call load() automatically (default True).

        Returns
        -------
        TrajectoryAnalysis

        Examples
        --------
            ana = manager.analyze_trajectories()
            ana.summary()
            ana.save()
            ana.plot_population("population.png")
        """
        if base_folder is None:
            wf = self.get_workflow("batch_from_wigner")
            if wf and "base_folder" in wf:
                base_folder = wf["base_folder"]
            else:
                base_folder = self.folder

        ana = _analysis.TrajectoryAnalysis(
            base_folder=base_folder,
            title=title or self.title,
        )
        if load:
            ana.load()
        return ana

    # ================================================================
    #  Batch NAMD (random velocity)
    # ================================================================
    def batch(
        self, n_traj, base_folder=None, seed_start=1,
        scheduler="slurm", submit=False, job_name=None,
        partition=None, nodes=1, ntasks=1, cpus_per_task=16,
        mem=None, time="24:00:00", account=None, modules=None,
        pre_commands=None, post_commands=None, python_cmd="python",
        extra_headers=None,
    ):
        """Create N trajectory folders with random seeds."""
        if self.mode not in {"namd", "namd_qmmm"}:
            raise ValueError(f"batch() is for NAMD modes, not '{self.mode}'")

        base = Path(base_folder) if base_folder else Path(self.folder)
        base.mkdir(parents=True, exist_ok=True)

        job_prefix = job_name or self.title
        modules = modules or []
        pre_commands = pre_commands or []
        post_commands = post_commands or []
        extra_headers = extra_headers or []

        folders = []
        scripts = []

        for i in range(n_traj):
            traj_idx = i + 1
            seed = seed_start + i
            traj_name = f"{self.title}_traj_{traj_idx:04d}"
            traj_folder = base / traj_name

            self.set_namd_config("control", "gl_seed", seed)
            self.folder = str(traj_folder)
            self._setup_engines()
            self._setup_files()

            run_script = _batch.generate_run_script(self.title)
            (traj_folder / "run.py").write_text(run_script, encoding="utf-8")

            script_content = _batch.generate_job_script(
                traj_folder=traj_folder, traj_name=traj_name,
                traj_idx=traj_idx, job_prefix=job_prefix,
                scheduler=scheduler, partition=partition,
                nodes=nodes, ntasks=ntasks, cpus_per_task=cpus_per_task,
                mem=mem, time=time, account=account, modules=modules,
                pre_commands=pre_commands, post_commands=post_commands,
                python_cmd=python_cmd, extra_headers=extra_headers,
            )
            script_path = traj_folder / "submit.sh"
            script_path.write_text(script_content, encoding="utf-8")
            os.chmod(script_path, 0o755)

            folders.append(traj_folder)
            scripts.append(script_path)

        self.folder = str(base)

        print(f"Created {n_traj} trajectory folders in: {base}")
        print(f"  Seed range: {seed_start} \u2192 {seed_start + n_traj - 1}")
        print(f"  Folders: {self.title}_traj_0001 \u2192 {self.title}_traj_{n_traj:04d}")

        if submit:
            _batch.submit_jobs(scripts, scheduler)

        master = _batch.generate_master_script(base, scripts, scheduler)
        master_path = base / "submit_all.sh"
        master_path.write_text(master, encoding="utf-8")
        os.chmod(master_path, 0o755)
        print(f"  Master: {master_path}")

        return folders
