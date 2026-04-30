"""
Trajectory analysis for PyRAI2MD NAMD ensembles.

File formats parsed
-------------------
title.md.energies  — time(au)  Epot  Ekin  Etot  E_s1  E_s2  ... (Hartree)
title.md.xyz       — multi-frame XYZ; comment: "title coord N state S"
title.sh.xyz       — hop-event frames; comment: "title coord N state F to T CI"
"""

import re
import numpy as np
from pathlib import Path

AU_TO_FS = 0.02418884      # 1 a.u. of time → femtoseconds
HARTREE_TO_EV = 27.211386  # 1 Hartree → eV


class TrajectoryAnalysis:
    """
    Analyse an ensemble of PyRAI2MD NAMD trajectories.

    Parameters
    ----------
    base_folder : str or Path
        Directory containing trajectory subfolders.
    title : str, optional
        Molecule title used in filenames (e.g. 'test').
        Auto-detected from folder contents if omitted.

    Examples
    --------
    Direct use::

        ana = TrajectoryAnalysis("namd_trajectories", title="mol")
        ana.load()
        ana.summary()
        ana.save()
        ana.plot_population("population.png")

    Via manager (preferred)::

        ana = manager.analyze_trajectories()
        ana.summary()
    """

    def __init__(self, base_folder, title=None):
        self.base_folder = Path(base_folder)
        self.title = title
        self._trajs = {}
        self._loaded = False

    # ------------------------------------------------------------------ #
    #  Loading                                                             #
    # ------------------------------------------------------------------ #

    def load(self):
        """
        Load all trajectories from *base_folder*.

        Returns
        -------
        dict
            n_loaded, n_failed, n_steps_min, n_steps_max, states_found.
        """
        self._trajs = {}
        subdirs = sorted(d for d in self.base_folder.iterdir() if d.is_dir())
        if not subdirs:
            raise FileNotFoundError(f"No subdirectories in {self.base_folder}")

        if self.title is None:
            for d in subdirs:
                for f in d.glob("*.md.energies"):
                    self.title = f.name.replace(".md.energies", "")
                    break
                if self.title:
                    break
            if self.title is None:
                raise FileNotFoundError("No .md.energies files found; set title= explicitly")

        n_failed = 0
        all_states = set()

        for d in subdirs:
            energy_file = d / f"{self.title}.md.energies"
            xyz_file    = d / f"{self.title}.md.xyz"
            sh_file     = d / f"{self.title}.sh.xyz"

            if not energy_file.is_file():
                continue

            time_au, epot, ekin, etot, state_e = _parse_energies(energy_file)
            if len(time_au) == 0:
                n_failed += 1
                continue

            states = _parse_states_xyz(xyz_file)
            if len(states) == 0:
                n_failed += 1
                continue

            hops = _parse_hops(sh_file)

            n = min(len(time_au), len(states))
            entry = {
                "time_au": time_au[:n],
                "time_fs": time_au[:n] * AU_TO_FS,
                "state":   states[:n],
                "epot":    epot[:n],
                "ekin":    ekin[:n],
                "etot":    etot[:n],
                "state_e": state_e[:n] if state_e is not None else None,
                "hops":    hops,
                "n_steps": n,
            }
            all_states.update(entry["state"].tolist())
            self._trajs[d.name] = entry

        n_loaded = len(self._trajs)
        if n_loaded == 0:
            raise RuntimeError("No valid trajectories loaded.")

        steps = [d["n_steps"] for d in self._trajs.values()]
        self._loaded = True

        summary = {
            "n_loaded":    n_loaded,
            "n_failed":    n_failed,
            "n_steps_min": min(steps),
            "n_steps_max": max(steps),
            "states_found": sorted(all_states),
        }
        print(f"Loaded {n_loaded} trajectories from: {self.base_folder}")
        print(f"  Failed/empty:  {n_failed}")
        print(f"  Steps range:   {min(steps)}–{max(steps)}")
        print(f"  States found:  {sorted(all_states)}")
        return summary

    def _require_loaded(self):
        if not self._loaded:
            raise RuntimeError("Call load() before running analysis.")

    # ------------------------------------------------------------------ #
    #  Properties                                                          #
    # ------------------------------------------------------------------ #

    @property
    def n_traj(self):
        return len(self._trajs)

    @property
    def trajectory_names(self):
        return list(self._trajs.keys())

    # ------------------------------------------------------------------ #
    #  Analysis methods                                                    #
    # ------------------------------------------------------------------ #

    def population(self, time_unit="fs"):
        """
        Ensemble-average state populations vs time.

        At step *i*, population of state S = (number of trajectories in S
        at step i) / (number of trajectories that have reached step i).

        Parameters
        ----------
        time_unit : str
            'fs' (default) or 'au'.

        Returns
        -------
        dict
            time, populations {state: array}, n_alive, states.
        """
        self._require_loaded()
        max_steps = max(d["n_steps"] for d in self._trajs.values())
        ref = max(self._trajs.values(), key=lambda d: d["n_steps"])
        time = ref["time_fs"] if time_unit == "fs" else ref["time_au"]

        all_states = sorted({s for d in self._trajs.values() for s in d["state"].tolist()})
        counts  = {s: np.zeros(max_steps) for s in all_states}
        n_alive = np.zeros(max_steps)

        for d in self._trajs.values():
            n = d["n_steps"]
            n_alive[:n] += 1
            for i, s in enumerate(d["state"]):
                counts[s][i] += 1

        with np.errstate(invalid="ignore"):
            populations = {
                s: np.where(n_alive > 0, counts[s] / n_alive, np.nan)
                for s in all_states
            }

        return {"time": time, "populations": populations,
                "n_alive": n_alive, "states": all_states}

    def average_energy(self, unit="eV", time_unit="fs"):
        """
        Ensemble-average energies vs time.

        Parameters
        ----------
        unit : str
            'eV' (default) or 'au'.
        time_unit : str
            'fs' (default) or 'au'.

        Returns
        -------
        dict
            time, epot, ekin, etot, state_energies {1-idx: array}, n_alive, unit.
        """
        self._require_loaded()
        max_steps = max(d["n_steps"] for d in self._trajs.values())
        ref = max(self._trajs.values(), key=lambda d: d["n_steps"])
        time = ref["time_fs"] if time_unit == "fs" else ref["time_au"]
        conv = HARTREE_TO_EV if unit == "eV" else 1.0

        epot_s = np.zeros(max_steps)
        ekin_s = np.zeros(max_steps)
        etot_s = np.zeros(max_steps)
        n_alive = np.zeros(max_steps)

        n_states = max(
            (d["state_e"].shape[1] for d in self._trajs.values() if d["state_e"] is not None),
            default=0,
        )
        se_sum = np.zeros((max_steps, n_states))

        for d in self._trajs.values():
            n = d["n_steps"]
            n_alive[:n] += 1
            epot_s[:n]  += d["epot"]
            ekin_s[:n]  += d["ekin"]
            etot_s[:n]  += d["etot"]
            if d["state_e"] is not None:
                ns = min(d["state_e"].shape[1], n_states)
                se_sum[:n, :ns] += d["state_e"][:, :ns]

        with np.errstate(invalid="ignore"):
            safe = n_alive > 0
            epot = np.where(safe, conv * epot_s / n_alive, np.nan)
            ekin = np.where(safe, conv * ekin_s / n_alive, np.nan)
            etot = np.where(safe, conv * etot_s / n_alive, np.nan)
            state_energies = {
                i + 1: np.where(safe, conv * se_sum[:, i] / n_alive, np.nan)
                for i in range(n_states)
            }

        return {"time": time, "epot": epot, "ekin": ekin, "etot": etot,
                "state_energies": state_energies, "n_alive": n_alive, "unit": unit}

    def hop_statistics(self):
        """
        Surface hopping event summary.

        Returns
        -------
        dict
            total_hops, hops_per_traj, from_to_counts,
            first_hop_times (fs), hop_events list.
        """
        self._require_loaded()
        total = 0
        per_traj = {}
        from_to = {}
        first_hop_times = []
        events = []

        for name, d in self._trajs.items():
            hops = d["hops"]
            per_traj[name] = len(hops)
            total += len(hops)

            if hops:
                step_i = hops[0]["step"] - 1
                if step_i < len(d["time_fs"]):
                    first_hop_times.append(float(d["time_fs"][step_i]))

            for h in hops:
                key = (h["from"], h["to"])
                from_to[key] = from_to.get(key, 0) + 1
                step_i = h["step"] - 1
                t_fs = float(d["time_fs"][step_i]) if step_i < len(d["time_fs"]) else np.nan
                events.append({
                    "traj": name, "step": h["step"],
                    "from": h["from"], "to": h["to"], "time_fs": t_fs,
                })

        return {
            "total_hops":      total,
            "hops_per_traj":   per_traj,
            "from_to_counts":  from_to,
            "first_hop_times": first_hop_times,
            "hop_events":      events,
        }

    def survival_probability(self, initial_state=None, time_unit="fs"):
        """
        Fraction of trajectories still in *initial_state* vs time.

        Parameters
        ----------
        initial_state : int, optional
            Defaults to the state at step 0 of the first trajectory.
        time_unit : str
            'fs' (default) or 'au'.

        Returns
        -------
        dict
            time, survival, n_alive, initial_state.
        """
        self._require_loaded()
        if initial_state is None:
            initial_state = int(next(iter(self._trajs.values()))["state"][0])

        max_steps = max(d["n_steps"] for d in self._trajs.values())
        ref = max(self._trajs.values(), key=lambda d: d["n_steps"])
        time = ref["time_fs"] if time_unit == "fs" else ref["time_au"]

        in_s0   = np.zeros(max_steps)
        n_alive = np.zeros(max_steps)

        for d in self._trajs.values():
            n = d["n_steps"]
            n_alive[:n] += 1
            in_s0[:n] += (d["state"] == initial_state).astype(float)

        with np.errstate(invalid="ignore"):
            survival = np.where(n_alive > 0, in_s0 / n_alive, np.nan)

        return {"time": time, "survival": survival,
                "n_alive": n_alive, "initial_state": initial_state}

    # ------------------------------------------------------------------ #
    #  Plotting                                                            #
    # ------------------------------------------------------------------ #

    def plot_population(self, output_file=None, figsize=(9, 5), time_unit="fs", title=None):
        """
        Plot state populations vs time.

        Parameters
        ----------
        output_file : str or Path, optional
            If given, saves PNG; otherwise shows interactively.
        """
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        res = self.population(time_unit=time_unit)
        colors = [
            "#e63946", "#457b9d", "#2a9d8f", "#e9c46a", "#f4a261",
            "#264653", "#a855f7", "#ef4444", "#06b6d4", "#84cc16",
        ]
        fig, ax = plt.subplots(figsize=figsize)
        for i, s in enumerate(res["states"]):
            ax.plot(res["time"], res["populations"][s],
                    color=colors[i % len(colors)], linewidth=2, label=f"State {s}")
        ax.set_xlabel(f"Time ({time_unit})", fontsize=12)
        ax.set_ylabel("Population", fontsize=12)
        ax.set_ylim(0, 1.05)
        ax.set_title(title or f"State populations ({self.n_traj} trajectories)", fontsize=13)
        ax.legend(frameon=False, fontsize=10)
        fig.tight_layout()
        _save_or_show(fig, output_file, "Population plot")

    def plot_energy(self, output_file=None, figsize=(9, 5), unit="eV",
                    time_unit="fs", show_states=True, title=None):
        """
        Plot ensemble-average energies vs time.

        Parameters
        ----------
        output_file : str or Path, optional
        unit : str
            'eV' (default) or 'au'.
        show_states : bool
            Include individual state energy curves.
        """
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        res = self.average_energy(unit=unit, time_unit=time_unit)
        colors = [
            "#e63946", "#457b9d", "#2a9d8f", "#e9c46a", "#f4a261",
            "#264653", "#a855f7", "#ef4444", "#06b6d4", "#84cc16",
        ]
        fig, ax = plt.subplots(figsize=figsize)

        if show_states:
            for i, (s_idx, s_e) in enumerate(res["state_energies"].items()):
                ax.plot(res["time"], s_e, color=colors[i % len(colors)],
                        linewidth=1.5, alpha=0.7, label=f"E State {s_idx}")

        ax.plot(res["time"], res["epot"], "k--", linewidth=2, label="E_pot (active)", alpha=0.85)
        ax.plot(res["time"], res["ekin"], "b:",  linewidth=1.5, label="E_kin", alpha=0.85)
        ax.set_xlabel(f"Time ({time_unit})", fontsize=12)
        ax.set_ylabel(f"Energy ({unit})", fontsize=12)
        ax.set_title(title or f"Average energy ({self.n_traj} trajectories)", fontsize=13)
        ax.legend(frameon=False, fontsize=9)
        fig.tight_layout()
        _save_or_show(fig, output_file, "Energy plot")

    def plot_survival(self, output_file=None, figsize=(9, 5), time_unit="fs",
                      initial_state=None, title=None):
        """Plot survival probability of the initial state vs time."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        res = self.survival_probability(initial_state=initial_state, time_unit=time_unit)
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(res["time"], res["survival"], "#e63946", linewidth=2)
        ax.set_xlabel(f"Time ({time_unit})", fontsize=12)
        ax.set_ylabel("Survival probability", fontsize=12)
        ax.set_ylim(0, 1.05)
        ax.set_title(
            title or f"Survival of S{res['initial_state']} ({self.n_traj} trajectories)",
            fontsize=13,
        )
        fig.tight_layout()
        _save_or_show(fig, output_file, "Survival plot")

    # ------------------------------------------------------------------ #
    #  Save & summary                                                      #
    # ------------------------------------------------------------------ #

    def save(self, output_dir=None, time_unit="fs", energy_unit="eV"):
        """
        Save all analysis to .dat files.

        Writes: population.dat, energy.dat, hop_statistics.dat, survival.dat

        Parameters
        ----------
        output_dir : str or Path, optional
            Default: base_folder.
        """
        self._require_loaded()
        out = Path(output_dir) if output_dir else self.base_folder
        out.mkdir(parents=True, exist_ok=True)

        # ── population ──
        pop = self.population(time_unit=time_unit)
        states = pop["states"]
        cols = [pop["time"]] + [pop["populations"][s] for s in states] + [pop["n_alive"]]
        names = ([f"Time({time_unit})"] + [f"Pop_S{s}" for s in states] + ["N_alive"])
        np.savetxt(
            out / "population.dat", np.column_stack(cols),
            header=f"# {self.n_traj} trajectories\n# Columns: {', '.join(names)}",
            fmt="%14.6f",
        )
        print(f"  Population   → {out / 'population.dat'}")

        # ── energy ──
        eng = self.average_energy(unit=energy_unit, time_unit=time_unit)
        e_cols  = [eng["time"], eng["epot"], eng["ekin"], eng["etot"]]
        e_names = [f"Time({time_unit})", f"Epot({energy_unit})",
                   f"Ekin({energy_unit})", f"Etot({energy_unit})"]
        for s_idx, s_e in eng["state_energies"].items():
            e_cols.append(s_e)
            e_names.append(f"E_S{s_idx}({energy_unit})")
        np.savetxt(
            out / "energy.dat", np.column_stack(e_cols),
            header=f"# Average energies\n# Columns: {', '.join(e_names)}",
            fmt="%18.8f",
        )
        print(f"  Energy       → {out / 'energy.dat'}")

        # ── hop statistics ──
        hops = self.hop_statistics()
        with open(out / "hop_statistics.dat", "w") as f:
            f.write(f"# Hop statistics — {self.n_traj} trajectories\n")
            f.write(f"# Total hops: {hops['total_hops']}\n#\n")
            f.write("# Channel counts\n# from  to  count\n")
            for (fr, to), cnt in sorted(hops["from_to_counts"].items()):
                f.write(f"  {fr:4d}  {to:4d}  {cnt:6d}\n")
            f.write("\n# Individual events\n# traj  step  from  to  time_fs\n")
            for ev in hops["hop_events"]:
                f.write(f"  {ev['traj']}  {ev['step']:5d}  "
                        f"{ev['from']:4d}  {ev['to']:4d}  {ev['time_fs']:12.4f}\n")
        print(f"  Hop stats    → {out / 'hop_statistics.dat'}")

        # ── survival ──
        surv = self.survival_probability(time_unit=time_unit)
        np.savetxt(
            out / "survival.dat",
            np.column_stack([surv["time"], surv["survival"], surv["n_alive"]]),
            header=(f"# Survival of S{surv['initial_state']}\n"
                    f"# Columns: Time({time_unit}), Survival, N_alive"),
            fmt="%14.6f",
        )
        print(f"  Survival     → {out / 'survival.dat'}")

        return {
            "population_file":    str(out / "population.dat"),
            "energy_file":        str(out / "energy.dat"),
            "hop_statistics_file":str(out / "hop_statistics.dat"),
            "survival_file":      str(out / "survival.dat"),
        }

    def summary(self):
        """Print a concise text summary of the ensemble."""
        self._require_loaded()
        hops = self.hop_statistics()
        surv = self.survival_probability()
        pop  = self.population()

        final_pop = {s: float(np.nanmean(v[-5:])) for s, v in pop["populations"].items()}
        final_surv = float(np.nanmean(surv["survival"][-5:]))
        t_max = float(pop["time"][-1])

        print(f"\nEnsemble: {self.base_folder}  ({self.n_traj} trajectories)")
        print(f"  Time range:       0 – {t_max:.1f} fs")
        print(f"  Total hops:       {hops['total_hops']}")
        if hops["from_to_counts"]:
            for (fr, to), cnt in sorted(hops["from_to_counts"].items()):
                print(f"    S{fr} → S{to}: {cnt}")
        if hops["first_hop_times"]:
            print(f"  Mean first-hop:   {np.mean(hops['first_hop_times']):.2f} fs")
        print(f"  Final populations:")
        for s, p in sorted(final_pop.items()):
            print(f"    S{s}: {p:.3f}")
        print(f"  Survival S{surv['initial_state']}:      {final_surv:.3f}")


# ------------------------------------------------------------------ #
#  Module-level file parsers                                          #
# ------------------------------------------------------------------ #

def _parse_energies(filepath):
    """Return (time_au, epot, ekin, etot, state_e_or_None) as ndarrays."""
    time_l, epot_l, ekin_l, etot_l, se_l = [], [], [], [], []
    with open(filepath, "r") as f:
        for line in f:
            if not line.strip() or "time" in line.lower():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                time_l.append(float(parts[0]))
                epot_l.append(float(parts[1]))
                ekin_l.append(float(parts[2]))
                etot_l.append(float(parts[3]))
                if len(parts) > 4:
                    se_l.append([float(x) for x in parts[4:]])
            except ValueError:
                continue

    time_au = np.array(time_l)
    epot    = np.array(epot_l)
    ekin    = np.array(ekin_l)
    etot    = np.array(etot_l)

    if se_l:
        max_n  = max(len(r) for r in se_l)
        state_e = np.array([r + [np.nan] * (max_n - len(r)) for r in se_l])
    else:
        state_e = None

    return time_au, epot, ekin, etot, state_e


_STATE_RE = re.compile(r"state\s+(\d+)")


def _parse_states_xyz(xyz_file):
    """
    Extract current state at each step from .md.xyz comment lines.

    Comment format: 'title coord N state S'
    Returns int ndarray of length = number of frames.
    """
    if not xyz_file.is_file():
        return np.array([], dtype=int)

    states = []
    with open(xyz_file, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        try:
            natoms = int(stripped)
        except ValueError:
            i += 1
            continue
        i += 1
        if i >= len(lines):
            break
        comment = lines[i]
        m = _STATE_RE.search(comment)
        states.append(int(m.group(1)) if m else 0)
        i += 1 + natoms

    return np.array(states, dtype=int)


_HOP_RE = re.compile(r"coord\s+(\d+)\s+state\s+(\d+)\s+to\s+(\d+)")


def _parse_hops(sh_file):
    """
    Parse hop events from .sh.xyz.

    Returns list of dicts: {step, from, to}.
    """
    if not sh_file.is_file():
        return []
    hops = []
    with open(sh_file, "r") as f:
        for line in f:
            m = _HOP_RE.search(line)
            if m:
                hops.append({
                    "step": int(m.group(1)),
                    "from": int(m.group(2)),
                    "to":   int(m.group(3)),
                })
    return hops


def _save_or_show(fig, output_file, label):
    import matplotlib.pyplot as plt
    if output_file:
        fig.savefig(output_file, dpi=300)
        plt.close(fig)
        print(f"{label} → {output_file}")
    else:
        plt.show()
