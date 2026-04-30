"""
Wigner sampling from Molden frequency files.

Provides:
    - build_hess_input()  — build an OpenQP Hessian input from sections
    - run_hessian()       — run OpenQP to produce the .freq.molden file
    - wigner_sampling()   — sample displaced geometries and velocities
"""

import numpy as np
from copy import deepcopy
from math import cosh, sinh
from pathlib import Path

from .constants import (
    CM1_TO_AU, AMU_TO_AU, KB_AU, BOHR_TO_ANG,
    ELEMENT_TO_MASS,
)
from .parsers import parse_molden, parse_refxyz
from .builders import sections_to_inp
from .engines import OpenQPEngine
from .templates import WIGNER_HESS_TEMPLATE
from .utils import resolve_content


def build_hess_input(sections, xyz, hess_sections=None, hess_state=1):
    """
    Build an OpenQP Hessian input from electronic-structure sections.

    Merges the caller's sections on top of WIGNER_HESS_TEMPLATE, strips
    dynamics-only sections (properties, nac, qmmm), forces runtype=hess.

    Parameters
    ----------
    sections : dict
        The manager's OpenQP section dict.
    xyz : str or Path
        XYZ content/file for the system= block.
    hess_sections : dict, optional
        Additional overrides for the Hessian input.
    hess_state : int
        Electronic state for the Hessian (default 1).

    Returns
    -------
    str
        Complete OpenQP .inp file content.
    """
    hess_secs = deepcopy(WIGNER_HESS_TEMPLATE)

    _skip = {"properties", "nac", "qmmm", "hess"}
    for sec_name, sec_params in sections.items():
        if sec_name in _skip:
            continue
        if sec_name in hess_secs:
            hess_secs[sec_name].update(deepcopy(sec_params))
        else:
            hess_secs[sec_name] = deepcopy(sec_params)

    hess_secs.setdefault("input", {})
    hess_secs["input"]["runtype"] = "hess"
    hess_secs["hess"] = {"state": hess_state}

    if hess_sections:
        for sec_name, sec_params in hess_sections.items():
            if sec_name in hess_secs:
                hess_secs[sec_name].update(sec_params)
            else:
                hess_secs[sec_name] = deepcopy(sec_params)

    return sections_to_inp(hess_secs, xyz=xyz)


def run_hessian(sections, xyz, title, folder,
                hess_sections=None, hess_state=1):
    """
    Run an OpenQP Hessian/frequency calculation.

    Parameters
    ----------
    sections : dict
        The manager's OpenQP sections.
    xyz : str or Path
        XYZ coordinates.
    title : str
        Project name.
    folder : str or Path
        Working directory for the Hessian.
    hess_sections : dict, optional
        Extra section overrides.
    hess_state : int
        Electronic state for the Hessian.

    Returns
    -------
    Path
        Path to the generated .freq.molden file.
    """
    hess_folder = Path(folder)
    hess_folder.mkdir(parents=True, exist_ok=True)
    hess_title = f"{title}_hess"

    inp_content = build_hess_input(sections, xyz, hess_sections, hess_state)

    print(f"Running OpenQP Hessian calculation in: {hess_folder}")
    print(f"  Input: {hess_title}.inp")

    engine = OpenQPEngine()
    xyz_text, _ = resolve_content(xyz)
    engine.setup_input_with_xyz(
        folder=str(hess_folder),
        openqp_content=inp_content,
        xyz_content=xyz_text,
        title=hess_title,
    )
    engine.prepare()
    engine.run()

    molden_path = hess_folder / f"{hess_title}.freq.molden"
    if not molden_path.is_file():
        candidates = list(hess_folder.glob("*.freq.molden"))
        if candidates:
            molden_path = candidates[0]
        else:
            raise FileNotFoundError(
                f"Hessian completed but no .freq.molden found in {hess_folder}. "
                f"Check {hess_title}.log for errors."
            )

    print(f"  Hessian complete → {molden_path}")
    return molden_path


def wigner_sampling(
    molden_file=None,
    nsamples=1000,
    temp=298.15,
    seed=1,
    scale=1.0,
    refxyz=None,
    skip_first=0,
    skip_last=0,
    output_dir=".",
    sections=None,
    xyz=None,
    title="mol",
    hess_sections=None,
    hess_state=1,
    hess_dir=None,
):
    """
    Generate Wigner-sampled initial conditions.

    If ``molden_file`` is given it is used directly; otherwise an
    OpenQP Hessian is run using ``sections`` and ``xyz``.

    Parameters
    ----------
    molden_file : str or Path, optional
        Pre-computed Molden frequency file.
    nsamples : int
        Number of samples (default 1000).
    temp : float
        Temperature in Kelvin (default 298.15).
    seed : int
        Random seed (default 1).
    scale : float
        Frequency scaling factor (default 1.0).
    refxyz : str or Path, optional
        Reference XYZ for element symbols (defaults to ``xyz``).
    skip_first, skip_last : int
        Modes to skip from the beginning / end.
    output_dir : str or Path
        Directory for output files.
    sections : dict, optional
        OpenQP sections (for auto-Hessian).
    xyz : str or Path, optional
        XYZ geometry (for auto-Hessian and default refxyz).
    title : str
        Project name (for auto-Hessian).
    hess_sections : dict, optional
        Extra overrides for the Hessian input.
    hess_state : int
        Electronic state for the Hessian (default 1).
    hess_dir : str or Path, optional
        Directory for the Hessian run.

    Returns
    -------
    dict
        Keys: geom_file, velo_file, molden_file, nsamples, natoms, nmodes.
    """
    # ── Obtain the Molden file ──
    if molden_file is not None:
        molden_path = Path(molden_file)
        if not molden_path.is_file():
            raise FileNotFoundError(f"Molden file not found: {molden_file}")
    else:
        if sections is None or xyz is None:
            raise ValueError(
                "Either provide molden_file, or both sections and xyz "
                "so a Hessian can be computed automatically."
            )
        hess_folder = Path(hess_dir) if hess_dir else Path(output_dir) / "hessian"
        molden_path = run_hessian(
            sections, xyz, title, hess_folder, hess_sections, hess_state,
        )

    # ── Parse Molden ──
    atoms, coords_bohr, freqs_cm1, modes_cart = parse_molden(str(molden_path))
    natoms = len(atoms)
    nmodes = len(freqs_cm1)

    print(f"Wigner sampling: {natoms} atoms, {nmodes} modes from {molden_path.name}")
    print(f"  Temperature: {temp} K | Samples: {nsamples} | Scale: {scale} | Seed: {seed}")

    # ── Element symbols from reference XYZ ──
    ref = refxyz if refxyz is not None else xyz
    if ref is not None:
        ref_path = Path(str(ref))
        if ref_path.is_file():
            ref_atoms = parse_refxyz(str(ref_path))
            if len(ref_atoms) == natoms:
                atoms = ref_atoms
            else:
                print(f"  Warning: reference XYZ has {len(ref_atoms)} atoms, "
                      f"Molden has {natoms}. Using Molden symbols.")

    # ── Masses in atomic units ──
    masses_au = np.array(
        [ELEMENT_TO_MASS.get(a, 1.0) * AMU_TO_AU for a in atoms], dtype=float,
    )

    # ── Frequencies ──
    omega_au = freqs_cm1 * scale * CM1_TO_AU
    kbt_au = KB_AU * temp if temp > 0 else 0.0

    # ── Active modes ──
    end_idx = nmodes - skip_last if skip_last > 0 else nmodes
    active_modes = [i for i in range(skip_first, end_idx) if omega_au[i] > 1e-10]

    print(f"  Active modes: {len(active_modes)} "
          f"(skipping first {skip_first}, last {skip_last}, "
          f"zero-freq modes filtered)")

    # ── Mass-weight and normalise ──
    modes_mw = np.copy(modes_cart)
    for i in range(nmodes):
        for a in range(natoms):
            modes_mw[i, a, :] *= np.sqrt(masses_au[a])
        norm = np.linalg.norm(modes_mw[i])
        if norm > 1e-12:
            modes_mw[i] /= norm

    # ── Wigner widths ──
    sigma_Q = np.zeros(nmodes)
    sigma_P = np.zeros(nmodes)
    for i in active_modes:
        w = omega_au[i]
        if kbt_au > 0:
            x = w / (2.0 * kbt_au)
            coth_x = 1.0 if x > 500.0 else cosh(x) / sinh(x)
        else:
            coth_x = 1.0
        sigma_Q[i] = np.sqrt(0.5 / w * coth_x)
        sigma_P[i] = np.sqrt(0.5 * w * coth_x)

    # ── Sample ──
    rng = np.random.RandomState(int(seed))
    all_geoms = []
    all_velos = []

    for _ in range(nsamples):
        disp = np.zeros((natoms, 3))
        vel = np.zeros((natoms, 3))
        for i in active_modes:
            Qi = rng.normal(0.0, sigma_Q[i])
            Pi = rng.normal(0.0, sigma_P[i])
            for a in range(natoms):
                inv_sqrt_m = 1.0 / np.sqrt(masses_au[a])
                disp[a] += Qi * modes_mw[i, a] * inv_sqrt_m
                vel[a]  += Pi * modes_mw[i, a] * inv_sqrt_m

        all_geoms.append((coords_bohr + disp) * BOHR_TO_ANG)
        all_velos.append(vel)

    # ── Write outputs ──
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    geom_file = out_dir / "wigner_geom.xyz"
    velo_file = out_dir / "wigner_velo.xyz"

    with open(geom_file, "w") as f:
        for geom in all_geoms:
            f.write(f"{natoms}\n [Angstrom]\n")
            for a in range(natoms):
                f.write(f"{atoms[a]:<2s}  {geom[a,0]:24.16e}  "
                        f"{geom[a,1]:24.16e}  {geom[a,2]:24.16e}\n")

    with open(velo_file, "w") as f:
        for idx, vel in enumerate(all_velos):
            f.write(f"{idx + 1} [Bohr / time_au]\n")
            for a in range(natoms):
                f.write(f"  {vel[a,0]:24.16e}  {vel[a,1]:24.16e}  "
                        f"{vel[a,2]:24.16e}\n")

    print(f"\n  Wigner sampling complete: {nsamples} samples generated")
    print(f"  Geometries → {geom_file}")
    print(f"  Velocities → {velo_file}")

    return {
        "geom_file": str(geom_file),
        "velo_file": str(velo_file),
        "molden_file": str(molden_path),
        "nsamples": nsamples,
        "natoms": natoms,
        "nmodes": len(active_modes),
    }
