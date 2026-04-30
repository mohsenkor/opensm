"""
Batch job generation: SLURM/PBS script builders and submission helpers.
"""

import os
import subprocess
from pathlib import Path


def generate_run_script(title):
    """Generate a Python run.py script for a single PyRAI2MD trajectory."""
    return (
        f"#!/usr/bin/env python\n"
        f"import os, sys\n"
        f"os.chdir(os.path.dirname(os.path.abspath(__file__)))\n\n"
        f"from PyRAI2MD.pyrai2md import PYRAI2MD\n\n"
        f'title = "{title}"\n'
        f"pmd = PYRAI2MD(title)\n"
        f"pmd.run()\n"
    )


def generate_job_script(
    traj_folder, traj_name, traj_idx,
    job_prefix, scheduler, partition, nodes, ntasks,
    cpus_per_task, mem, time, account, modules,
    pre_commands, post_commands, python_cmd, extra_headers,
):
    """Generate a SLURM or PBS submission script."""
    abs_folder = Path(traj_folder).resolve()
    lines = ["#!/bin/bash"]

    if scheduler == "slurm":
        lines.append(f"#SBATCH --job-name={job_prefix}_{traj_idx:04d}")
        lines.append(f"#SBATCH --output={abs_folder}/slurm_%j.out")
        lines.append(f"#SBATCH --error={abs_folder}/slurm_%j.err")
        lines.append(f"#SBATCH --nodes={nodes}")
        lines.append(f"#SBATCH --ntasks={ntasks}")
        lines.append(f"#SBATCH --cpus-per-task={cpus_per_task}")
        lines.append(f"#SBATCH --time={time}")
        if partition:
            lines.append(f"#SBATCH --partition={partition}")
        if mem:
            lines.append(f"#SBATCH --mem={mem}")
        if account:
            lines.append(f"#SBATCH --account={account}")
        for h in extra_headers:
            lines.append(f"#SBATCH {h}")

    elif scheduler == "pbs":
        lines.append(f"#PBS -N {job_prefix}_{traj_idx:04d}")
        lines.append(f"#PBS -o {abs_folder}/pbs.out")
        lines.append(f"#PBS -e {abs_folder}/pbs.err")
        lines.append(f"#PBS -l nodes={nodes}:ppn={cpus_per_task}")
        lines.append(f"#PBS -l walltime={time}")
        if partition:
            lines.append(f"#PBS -q {partition}")
        if mem:
            lines.append(f"#PBS -l mem={mem}")
        if account:
            lines.append(f"#PBS -A {account}")
        for h in extra_headers:
            lines.append(f"#PBS {h}")

    lines.append("")

    for mod in modules:
        lines.append(f"module load {mod}")
    if modules:
        lines.append("")

    for cmd in pre_commands:
        lines.append(cmd)
    if pre_commands:
        lines.append("")

    lines.append(f"cd {abs_folder}")
    lines.append(f"{python_cmd} run.py")
    lines.append("")

    for cmd in post_commands:
        lines.append(cmd)

    return "\n".join(lines) + "\n"


def generate_master_script(base, scripts, scheduler):
    """Generate a submit_all.sh that submits every job."""
    cmd = "sbatch" if scheduler == "slurm" else "qsub"
    lines = ["#!/bin/bash", "# Submit all jobs", ""]
    for script in scripts:
        lines.append(f"{cmd} {script.resolve()}")
    lines.append("")
    lines.append(f'echo "Submitted {len(scripts)} jobs"')
    return "\n".join(lines) + "\n"


def submit_jobs(scripts, scheduler):
    """Actually submit all job scripts via sbatch/qsub."""
    cmd = "sbatch" if scheduler == "slurm" else "qsub"
    for script in scripts:
        try:
            result = subprocess.run(
                [cmd, str(script)],
                capture_output=True, text=True,
            )
            print(f"  Submitted: {script.name} → {result.stdout.strip()}")
            if result.returncode != 0:
                print(f"    Error: {result.stderr.strip()}")
        except FileNotFoundError:
            print(f"  Error: '{cmd}' not found. Are you on a login node?")
            break
