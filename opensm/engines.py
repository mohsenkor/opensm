"""
Compute engines: OpenQPEngine (QM) and PyRAI2MDEngine (NAMD).
"""

import os
from pathlib import Path

from .utils import resolve_content


class OpenQPEngine:
    """Wrapper around the OpenQP Runner for QM and QM/MM calculations."""

    Runner = None

    def __init__(self, config=None):
        self.config = config or {}
        self.project = self.config.get("project", "default")
        self.workdir = self.config.get("workdir", ".")
        self.input_dict = self.config.get("input_dict", {})
        self.back_door_data = self.config.get("back_door_data", {})
        self.results = None

    def setup_input(self, folder, openqp_content, title):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.workdir = str(folder)
        self.project = title
        text, _ = resolve_content(openqp_content)
        (folder / f"{title}.inp").write_text(text.strip() + "\n", encoding="utf-8")
        return folder

    def setup_input_with_xyz(self, folder, openqp_content, xyz_content, title):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.workdir = str(folder)
        self.project = title
        oqp_text, _ = resolve_content(openqp_content)
        (folder / f"{title}.inp").write_text(oqp_text.strip() + "\n", encoding="utf-8")
        xyz_text, _ = resolve_content(xyz_content)
        (folder / f"{title}.xyz").write_text(xyz_text.strip() + "\n", encoding="utf-8")
        return folder

    def prepare(self):
        if OpenQPEngine.Runner is None:
            from oqp.pyoqp import Runner
            OpenQPEngine.Runner = Runner
        Path(self.workdir).mkdir(parents=True, exist_ok=True)
        print(f"OpenQP engine prepared  (project={self.project}, workdir={self.workdir})")

    def run(self):
        if OpenQPEngine.Runner is None:
            self.prepare()
        pyoqp = OpenQPEngine.Runner(
            project=self.project,
            input_file=f"{self.workdir}/{self.project}.inp",
            input_dict=self.input_dict,
            log=f"{self.workdir}/{self.project}.log",
            silent=1, usempi=False,
        )
        original_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            pyoqp.back_door(self.back_door_data)
            pyoqp.run()
            self.results = pyoqp.results()
        finally:
            os.chdir(original_dir)
        print(f"OpenQP finished  (project={self.project})")
        return self.results

    def run_qmmm(self):
        from oqp.library.qmmm_md import QMMM_MD
        input_file = f"{self.project}.inp"
        original_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            print(f"Running QM/MM in: {self.workdir} (cfg={input_file})")
            md = QMMM_MD(oqp_cfg=input_file)
            md.run()
        finally:
            os.chdir(original_dir)
        print(f"QM/MM finished  (project={self.project})")


class PyRAI2MDEngine:
    """Wrapper around PyRAI2MD for non-adiabatic molecular dynamics."""

    def __init__(self, mode, config=None):
        self.mode = mode
        self.config = config or {}

    def prepare(self):
        print("Preparing PyRAI2MD engine")

    def _set_runtype(self, openqp_content):
        runtype = "prop" if self.mode == "namd" else "qmmm"
        lines, found = [], False
        for line in openqp_content.splitlines():
            if line.strip().startswith("runtype"):
                lines.append(f"runtype={runtype}")
                found = True
            else:
                lines.append(line)
        if not found:
            new_lines = []
            for line in lines:
                new_lines.append(line)
                if line.strip().lower() == "[input]":
                    new_lines.append(f"runtype={runtype}")
            lines = new_lines
        return "\n".join(lines)

    def setup_files(self, folder, input_content, openqp_content, xyz_content, title):
        """
        Write the three NAMD files into folder:
            title           — PyRAI2MD control file
            title.openqp    — OpenQP electronic structure input
            title.xyz       — molecular coordinates
        """
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)

        inp_text, _ = resolve_content(input_content)
        (folder / title).write_text(inp_text.strip() + "\n", encoding="utf-8")

        oqp_text, _ = resolve_content(openqp_content)
        oqp_text = self._set_runtype(oqp_text)
        (folder / f"{title}.openqp").write_text(oqp_text.strip() + "\n", encoding="utf-8")

        xyz_text, _ = resolve_content(xyz_content)
        (folder / f"{title}.xyz").write_text(xyz_text.strip() + "\n", encoding="utf-8")

        return folder

    def run(self, folder, title):
        from PyRAI2MD.pyrai2md import PYRAI2MD
        folder = Path(folder).resolve()
        original_dir = os.getcwd()
        try:
            os.chdir(folder)
            print(f"Running PyRAI2MD in: {folder}")
            pmd = PYRAI2MD(title)
            pmd.run()
        finally:
            os.chdir(original_dir)
