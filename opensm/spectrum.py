"""
Absorption spectrum workflow:
    - build_spectrum_input()  — build OpenQP single-point sections
    - absorption_batch()      — create calculation folders per Wigner sample
    - absorption_spectrum()   — read results and build broadened spectrum
    - select_wigner_window()  — filter samples by energy window
"""

import os
import numpy as np
from copy import deepcopy
from pathlib import Path

from .builders import sections_to_inp
from .engines import OpenQPEngine
from .parsers import read_xyz_frames, read_velo_frames, read_openqp_transitions


def build_spectrum_input(sections, spectrum_sections=None):
    """
    Build OpenQP sections for a single-point absorption calculation.

    Reuses electronic-structure sections, strips dynamics-only ones
    (properties, nac, qmmm, hess), forces runtype=energy.

    Parameters
    ----------
    sections : dict
        The manager's OpenQP section dict.
    spectrum_sections : dict, optional
        Extra overrides.

    Returns
    -------
    dict
        OpenQP sections dict with runtype=energy.
    """
    spec_secs = {}
    _skip = {"properties", "nac", "qmmm", "hess"}
    for sec_name, sec_params in sections.items():
        if sec_name in _skip:
            continue
        spec_secs[sec_name] = deepcopy(sec_params)

    spec_secs.setdefault("input", {})
    spec_secs["input"]["runtype"] = "energy"

    if spectrum_sections:
        for sec_name, sec_params in spectrum_sections.items():
            if sec_name in spec_secs:
                spec_secs[sec_name].update(sec_params)
            else:
                spec_secs[sec_name] = deepcopy(sec_params)

    return spec_secs


def absorption_batch(
    geom_file,
    sections,
    title,
    base_folder,
    n_samples=None,
    start_sample=0,
    spectrum_sections=None,
    job_script_func=None,
    scheduler="slurm",
    submit=False,
    job_name=None,
    partition=None,
    nodes=1,
    ntasks=1,
    cpus_per_task=16,
    mem=None,
    time="02:00:00",
    account=None,
    modules=None,
    pre_commands=None,
    post_commands=None,
    extra_headers=None,
    master_script_func=None,
    submit_func=None,
):
    """
    Create OpenQP single-point folders for each Wigner sample.

    Parameters
    ----------
    geom_file : str or Path
        Wigner geometry file (wigner_geom.xyz).
    sections : dict
        Manager's OpenQP sections.
    title : str
        Project name.
    base_folder : str or Path
        Parent directory for calculation folders.
    n_samples : int, optional
        Number of samples (default: all available).
    start_sample : int
        First sample index (default 0).
    spectrum_sections : dict, optional
        Extra OpenQP overrides for the single-point input.
    job_script_func : callable, optional
        Function(traj_folder, ...) → str to generate job scripts.
    scheduler, submit, job_name, partition, ... :
        Scheduler options (passed to job_script_func).
    master_script_func : callable, optional
        Function(base, scripts, scheduler) → str for master script.
    submit_func : callable, optional
        Function(scripts, scheduler) to submit jobs.

    Returns
    -------
    dict
        base_folder, folders, n_samples, scripts.
    """
    geom_frames = read_xyz_frames(geom_file)
    n_available = len(geom_frames)

    if n_samples is None:
        n_samples = n_available - start_sample
    if start_sample + n_samples > n_available:
        raise ValueError(
            f"Requested {n_samples} samples starting at {start_sample}, "
            f"but only {n_available} available."
        )

    base = Path(base_folder)
    base.mkdir(parents=True, exist_ok=True)

    job_prefix = job_name or f"{title}_abs"
    modules = modules or []
    pre_commands = pre_commands or []
    post_commands = post_commands or []
    extra_headers = extra_headers or []

    spec_secs = build_spectrum_input(sections, spectrum_sections)

    print(f"Creating {n_samples} absorption single-point calculations "
          f"[{start_sample}:{start_sample + n_samples}]")

    folders = []
    scripts = []

    for i in range(n_samples):
        sample_idx = start_sample + i
        sample_name = f"{title}_abs_{sample_idx + 1:04d}"
        sample_folder = base / sample_name
        sample_folder.mkdir(parents=True, exist_ok=True)

        atoms, coords_ang = geom_frames[sample_idx]
        natoms = len(atoms)

        xyz_lines = [f"{natoms}", ""]
        for a in range(natoms):
            xyz_lines.append(
                f"{atoms[a]:<2s}  {coords_ang[a][0]:20.12f}  "
                f"{coords_ang[a][1]:20.12f}  {coords_ang[a][2]:20.12f}"
            )
        xyz_content = "\n".join(xyz_lines)

        inp_content = sections_to_inp(spec_secs, xyz=xyz_content)

        engine = OpenQPEngine()
        engine.setup_input_with_xyz(
            folder=str(sample_folder),
            openqp_content=inp_content,
            xyz_content=xyz_content,
            title=sample_name,
        )

        # Run script
        run_script = (
            f"#!/usr/bin/env python\n"
            f"import os, sys\n"
            f"os.chdir(os.path.dirname(os.path.abspath(__file__)))\n\n"
            f"from oqp.pyoqp import Runner\n\n"
            f'title = "{sample_name}"\n'
            f"pyoqp = Runner(\n"
            f"    project=title,\n"
            f'    input_file=f"{{title}}.inp",\n'
            f'    log=f"{{title}}.log",\n'
            f"    silent=1, usempi=False,\n"
            f")\n"
            f"pyoqp.run()\n"
        )
        (sample_folder / "run.py").write_text(run_script, encoding="utf-8")

        # Job script
        if job_script_func:
            script_content = job_script_func(
                traj_folder=sample_folder, traj_name=sample_name,
                traj_idx=sample_idx + 1, job_prefix=job_prefix,
                scheduler=scheduler, partition=partition,
                nodes=nodes, ntasks=ntasks, cpus_per_task=cpus_per_task,
                mem=mem, time=time, account=account, modules=modules,
                pre_commands=pre_commands, post_commands=post_commands,
                python_cmd="python", extra_headers=extra_headers,
            )
            script_path = sample_folder / "submit.sh"
            script_path.write_text(script_content, encoding="utf-8")
            os.chmod(script_path, 0o755)
            scripts.append(script_path)

        folders.append(sample_folder)

    # Master script
    if master_script_func and scripts:
        master = master_script_func(base, scripts, scheduler)
        master_path = base / "submit_all.sh"
        master_path.write_text(master, encoding="utf-8")
        os.chmod(master_path, 0o755)
        print(f"  Master: {master_path}")

    print(f"\nCreated {n_samples} absorption folders in: {base}")
    print(f"  Folders: {title}_abs_0001 → {title}_abs_{n_samples:04d}")

    if submit and submit_func and scripts:
        submit_func(scripts, scheduler)

    return {
        "base_folder": str(base),
        "folders": folders,
        "n_samples": n_samples,
        "scripts": scripts,
    }


def absorption_spectrum(
    base_folder,
    initial_state=1,
    e_min=None,
    e_max=None,
    npts=5000,
    sigma=0.15,
    output_dir=None,
    plot=True,
    figsize=(9, 5.5),
):
    """
    Read OpenQP log files and build a state-resolved, Gaussian-broadened
    absorption spectrum.

    The spectrum is decomposed by final state (e.g. S₁, S₂, …) so you
    can identify which electronic band corresponds to which energy window
    before selecting samples for NAMD.

    Parameters
    ----------
    base_folder : str or Path
        Directory containing absorption calculation subfolders.
    initial_state : int
        Transitions originating from this state (default 1).
    e_min, e_max : float or None
        Energy range in eV.  If None, auto-detected from the data
        with a 1 eV padding on each side.
    npts : int
        Grid points (default 5000).
    sigma : float
        Gaussian broadening width in eV (default 0.15).
    output_dir : str or Path, optional
        Save directory (default: base_folder).
    plot : bool
        Save a PNG plot (default True).
    figsize : tuple
        Figure size (default (9, 5.5)).

    Returns
    -------
    dict
        - 'e_grid'           : ndarray — energy grid in eV
        - 'spectrum'         : ndarray — normalized total absorption
        - 'raw_spectrum'     : ndarray — unnormalized total absorption
        - 'state_spectra'    : dict    — {final_state: ndarray} per-state curves
        - 'all_transitions'  : dict    — {sample_idx: [(E, osc, fstate), ...]}
        - 'n_converged'      : int
        - 'n_failed'         : int
        - 'spectrum_file'    : str     — path to total spectrum .dat
        - 'transitions_file' : str     — path to per-sample transitions .dat
        - 'plot_file'        : str or None
    """
    base = Path(base_folder)
    if not base.is_dir():
        raise FileNotFoundError(f"Absorption folder not found: {base}")

    out_dir = Path(output_dir) if output_dir else base
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Find log files ──
    log_files = []
    for subfolder in sorted(base.iterdir()):
        if not subfolder.is_dir():
            continue
        for logfile in subfolder.glob("*.log"):
            log_files.append(logfile)
        for logfile in subfolder.glob("*.out"):
            log_files.append(logfile)

    if not log_files:
        raise FileNotFoundError(f"No .log/.out files in subfolders of {base}")

    print(f"Reading {len(log_files)} log files from: {base}")

    # ── Parse all transitions ──
    all_transitions = {}
    n_converged = 0
    n_failed = 0

    for logfile in sorted(log_files):
        transitions = read_openqp_transitions(str(logfile), initial_state)

        folder_name = logfile.parent.name
        try:
            idx_str = folder_name.rsplit("_", 1)[-1]
            sample_idx = int(idx_str) - 1
        except (ValueError, IndexError):
            sample_idx = n_converged + n_failed

        if not transitions:
            n_failed += 1
            print(f"  Warning: no transitions in {logfile.name}")
            continue

        n_converged += 1
        all_transitions[sample_idx] = transitions

    if not all_transitions:
        raise RuntimeError("No valid transitions found in any log file.")

    # ── Collect all energies to auto-detect range ──
    all_energies = []
    all_final_states = set()
    for transitions in all_transitions.values():
        for energy_ev, osc, fstate in transitions:
            all_energies.append(energy_ev)
            all_final_states.add(fstate)

    if e_min is None:
        e_min = max(0.0, min(all_energies) - 1.0)
    if e_max is None:
        e_max = max(all_energies) + 1.0

    e_grid = np.linspace(e_min, e_max, npts)
    final_states = sorted(all_final_states)

    print(f"  Converged: {n_converged} | Failed: {n_failed}")
    print(f"  Final states: {final_states}")
    print(f"  Energy range: [{e_min:.2f}, {e_max:.2f}] eV (auto)" if e_min is not None else "")

    # ── Build state-resolved spectra ──
    # For each final state, average the per-sample spectra
    state_spectra = {}

    for fstate in final_states:
        sample_spectra = []
        for sample_idx, transitions in all_transitions.items():
            y = np.zeros_like(e_grid)
            for energy_ev, osc, fs in transitions:
                if fs == fstate:
                    y += osc * np.exp(-0.5 * ((e_grid - energy_ev) / sigma) ** 2)
            sample_spectra.append(y)

        state_spectra[fstate] = np.mean(sample_spectra, axis=0)

    # ── Total spectrum ──
    total_spectrum = np.zeros_like(e_grid)
    for y in state_spectra.values():
        total_spectrum += y

    raw_spectrum = total_spectrum.copy()

    # ── Normalize all curves by total maximum ──
    max_val = np.max(total_spectrum)
    if max_val > 0:
        total_spectrum /= max_val
        for fstate in state_spectra:
            state_spectra[fstate] = state_spectra[fstate] / max_val

    # ── Save total spectrum data ──
    spec_file = out_dir / "absorption_spectrum.dat"
    # Build column stack: e_grid, total, then each state
    cols = [e_grid, total_spectrum]
    col_names = ["Energy(eV)", "Total"]
    for fstate in final_states:
        cols.append(state_spectra[fstate])
        col_names.append(f"State_{initial_state}->{fstate}")

    header = (
        f"# State-resolved absorption spectrum from {n_converged} Wigner samples\n"
        f"# sigma = {sigma} eV, initial_state = {initial_state}\n"
        f"# Columns: {', '.join(col_names)}"
    )
    np.savetxt(
        spec_file, np.column_stack(cols),
        header=header, fmt="%12.6f" + "  %18.10e" * (len(cols) - 1),
    )
    print(f"  Spectrum data  → {spec_file}")

    # ── Save per-sample transition data ──
    trans_file = out_dir / "absorption_transitions.dat"
    with open(trans_file, "w") as f:
        f.write("# sample_idx  initial_state  final_state  "
                "energy_eV  osc_strength\n")
        for sample_idx in sorted(all_transitions.keys()):
            for energy_ev, osc, fstate in all_transitions[sample_idx]:
                f.write(f"  {sample_idx:6d}  {initial_state:6d}  "
                        f"{fstate:6d}  {energy_ev:12.6f}  {osc:12.6f}\n")
    print(f"  Transitions    → {trans_file}")

    # ── State-resolved plot ──
    plot_file = None
    if plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            # Color palette for up to 10 states
            colors = [
                "#e63946", "#457b9d", "#2a9d8f", "#e9c46a", "#f4a261",
                "#264653", "#a855f7", "#ef4444", "#06b6d4", "#84cc16",
            ]

            fig, ax = plt.subplots(figsize=figsize)

            for i, fstate in enumerate(final_states):
                color = colors[i % len(colors)]
                ax.plot(
                    e_grid, state_spectra[fstate],
                    linewidth=1.6, color=color,
                    label=f"State {initial_state} \u2192 {fstate}",
                )
                ax.fill_between(
                    e_grid, state_spectra[fstate],
                    alpha=0.08, color=color,
                )

            ax.plot(
                e_grid, total_spectrum,
                linewidth=2.8, linestyle="--", color="#1e293b",
                label="Total",
            )

            ax.set_xlabel("Excitation energy (eV)", fontsize=12)
            ax.set_ylabel("Normalized absorption intensity", fontsize=12)
            ax.set_title(
                f"State-resolved absorption spectrum "
                f"({n_converged} samples, \u03c3 = {sigma} eV)",
                fontsize=13,
            )
            ax.set_xlim(e_min, e_max)
            ax.set_ylim(bottom=0)
            ax.legend(frameon=False, fontsize=10)
            fig.tight_layout()

            plot_file = str(out_dir / "absorption_spectrum.png")
            fig.savefig(plot_file, dpi=300)
            plt.close(fig)
            print(f"  Plot           → {plot_file}")
        except ImportError:
            print("  Warning: matplotlib not available, skipping plot.")

    return {
        "e_grid": e_grid,
        "spectrum": total_spectrum,
        "raw_spectrum": raw_spectrum,
        "state_spectra": state_spectra,
        "all_transitions": all_transitions,
        "n_converged": n_converged,
        "n_failed": n_failed,
        "spectrum_file": str(spec_file),
        "transitions_file": str(trans_file),
        "plot_file": plot_file,
    }


def select_wigner_window(
    geom_file,
    velo_file,
    all_transitions,
    e_min,
    e_max,
    target_state=None,
    output_dir=".",
):
    """
    Filter Wigner samples by an energy window in the absorption spectrum.

    Parameters
    ----------
    geom_file, velo_file : str or Path
        Wigner geometry and velocity files.
    all_transitions : dict
        {sample_idx: [(E_eV, osc, fstate), ...]}.
    e_min, e_max : float
        Energy window in eV.
    target_state : int, optional
        Only consider transitions to this final state.
    output_dir : str or Path
        Directory for filtered output files.

    Returns
    -------
    dict
        geom_file, velo_file, n_selected, n_total,
        selected_indices, selected_transitions.
    """
    geom_frames = read_xyz_frames(geom_file)
    velo_frames = read_velo_frames(velo_file)

    selected_indices = []
    selected_transitions = {}

    for sample_idx, transitions in sorted(all_transitions.items()):
        for energy_ev, osc, fstate in transitions:
            if target_state is not None and fstate != target_state:
                continue
            if e_min <= energy_ev <= e_max:
                selected_indices.append(sample_idx)
                selected_transitions[sample_idx] = transitions
                break

    n_selected = len(selected_indices)
    n_total = len(all_transitions)

    print(f"Energy window: [{e_min:.2f}, {e_max:.2f}] eV"
          + (f" (state {target_state})" if target_state else ""))
    print(f"  Selected {n_selected} / {n_total} samples")

    if n_selected == 0:
        print("  Warning: no samples found in the specified window.")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    filtered_geom = out_dir / "wigner_geom_selected.xyz"
    filtered_velo = out_dir / "wigner_velo_selected.xyz"

    with open(filtered_geom, "w") as f:
        for new_idx, sample_idx in enumerate(selected_indices):
            if sample_idx >= len(geom_frames):
                continue
            atoms, coords = geom_frames[sample_idx]
            natoms = len(atoms)
            f.write(f"{natoms}\n [Angstrom]\n")
            for a in range(natoms):
                f.write(f"{atoms[a]:<2s}  {coords[a,0]:24.16e}  "
                        f"{coords[a,1]:24.16e}  {coords[a,2]:24.16e}\n")

    with open(filtered_velo, "w") as f:
        for new_idx, sample_idx in enumerate(selected_indices):
            if sample_idx >= len(velo_frames):
                continue
            vel = velo_frames[sample_idx]
            natoms = vel.shape[0]
            f.write(f"{new_idx + 1} [Bohr / time_au]\n")
            for a in range(natoms):
                f.write(f"  {vel[a,0]:24.16e}  {vel[a,1]:24.16e}  "
                        f"{vel[a,2]:24.16e}\n")

    # Summary file
    summary_file = out_dir / "wigner_selection.dat"
    with open(summary_file, "w") as f:
        f.write(f"# Energy window: [{e_min:.4f}, {e_max:.4f}] eV\n")
        if target_state is not None:
            f.write(f"# Target state: {target_state}\n")
        f.write(f"# Selected {n_selected} / {n_total} samples\n#\n")
        f.write("# new_idx  original_idx  E(eV)_brightest  "
                "f_brightest  fstate\n")
        for new_idx, sample_idx in enumerate(selected_indices):
            trans = selected_transitions[sample_idx]
            best = None
            for energy_ev, osc, fstate in trans:
                if target_state is not None and fstate != target_state:
                    continue
                if e_min <= energy_ev <= e_max:
                    if best is None or osc > best[1]:
                        best = (energy_ev, osc, fstate)
            if best:
                f.write(f"  {new_idx:6d}  {sample_idx:6d}  "
                        f"{best[0]:12.6f}  {best[1]:12.6f}  {best[2]:4d}\n")

    print(f"  Filtered geometries → {filtered_geom}")
    print(f"  Filtered velocities → {filtered_velo}")
    print(f"  Selection summary   → {summary_file}")

    return {
        "geom_file": str(filtered_geom),
        "velo_file": str(filtered_velo),
        "n_selected": n_selected,
        "n_total": n_total,
        "selected_indices": selected_indices,
        "selected_transitions": selected_transitions,
    }
