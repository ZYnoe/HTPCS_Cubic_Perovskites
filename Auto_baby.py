import subprocess
import time
from pathlib import Path
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter

from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Kpoints

class Experiment:

    def __init__(self, material_dir: str = "") -> None:
        if not material_dir: raise ValueError("material_dir must be provided.")
        self.base_dir = Path(material_dir)
        self.scf_dir = self.base_dir / "1_scf"
        self.dos_dir = self.base_dir / "2_dos"
        self.band_dir = self.base_dir / "3_bs"
        self.job_script = "job.sh"
        self.submit_command = "qsub"

    def ensure_directories(self) -> None:
        """Ensure key calculation directories and the base directory exist."""
        for path in (self.base_dir, self.scf_dir, self.dos_dir, self.band_dir):
            path.mkdir(parents=True, exist_ok=True)
        print(f"All calculations are organized in directories under {self.base_dir} based on their names.")

    def ensure_required_files(self) -> None:
        """Check each calculation directory contains its required input files."""
        required_files_by_dir = {
            self.scf_dir: ("POSCAR", "POTCAR", "INCAR", "job.sh"),
            self.dos_dir: ("INCAR",),
            self.band_dir: ("INCAR","KPOINTS"),
        }
        missing_by_dir = {}

        for directory, files in required_files_by_dir.items():
            missing = [fname for fname in files if not (directory / fname).is_file()]
            if missing:
                missing_by_dir[directory] = missing

        if missing_by_dir:
            problems = "; ".join(f"{directory}: {', '.join(files)}" for directory, files in missing_by_dir.items())
            raise FileNotFoundError(f"Required files not found -> {problems}")

        print("All required calculation inputs and job scripts are present where needed.")

    def scf_calculation(self, interval: int = 1) -> None:

        # SCF calculation
        process = subprocess.Popen(
            f"cd {self.scf_dir} && {self.submit_command} {self.job_script} > job_id.txt",
            shell=True)
        process.wait(5)
        with open(f"{self.scf_dir}/job_id.txt") as f:
            job_id = f.read().strip()
        while True:
            ret = subprocess.check_output(
                f"sacct -j {job_id} --format=State --noheader", shell=True
            ).decode("utf-8").strip()
            main_state = ret.splitlines()[0].strip() if ret else "UNKNOWN"
            print(f"Current SLURM job state: {main_state}")
            if "COMPLETED" in main_state or "FAILED" in main_state or "CANCELLED" in main_state or "TIMEOUT" in main_state:
                print(f"Final SLURM job state: {main_state}")
                break

            time.sleep(interval)

    def dos(self, interval: int = 1) -> None:

        # DOS calculation
        process = subprocess.Popen(
            f"cd {self.dos_dir} && {self.submit_command} {self.job_script} > job_id.txt",
            shell=True)
        process.wait(5)
        with open(f"{self.dos_dir}/job_id.txt") as f:
            job_id = f.read().strip()
        while True:
            ret = subprocess.check_output(
                f"sacct -j {job_id} --format=State --noheader", shell=True
            ).decode("utf-8").strip()
            main_state = ret.splitlines()[0].strip() if ret else "UNKNOWN"
            print(f"Current SLURM job state: {main_state}")
            if "COMPLETED" in main_state or "FAILED" in main_state or "CANCELLED" in main_state or "TIMEOUT" in main_state:
                print(f"Final SLURM job state: {main_state}")
                break
            
            time.sleep(interval)

    def band_structure(self, interval: int = 1) -> None:

        # Band structure calculation
        process = subprocess.Popen(
            f"cd {self.band_dir} && {self.submit_command} {self.job_script} > job_id.txt",
            shell=True)
        process.wait(5)
        with open(f"{self.band_dir}/job_id.txt") as f:
            job_id = f.read().strip()
        while True:
            ret = subprocess.check_output(
                f"sacct -j {job_id} --format=State --noheader", shell=True
            ).decode("utf-8").strip()
            main_state = ret.splitlines()[0].strip() if ret else "UNKNOWN"
            print(f"Current SLURM job state: {main_state}")
            if "COMPLETED" in main_state or "FAILED" in main_state or "CANCELLED" in main_state or "TIMEOUT" in main_state:
                print(f"Final SLURM job state: {main_state}")
                break

            time.sleep(interval)

    def plot_dos(self, material_path: str = ''):
        if not material_path:
            raise ValueError("material_path must be provided.")
        # load data
        result = Vasprun(f'{material_path}/2_dos/vasprun.xml', parse_potcar_file=False)
        # complete_dos = result.complete_dos
        # pdos_Si = complete_dos.get_element_spd_dos('Si')

        plotter = DosPlotter()
        plotter.add_dos('Total DOS', result.tdos)
        # plotter.add_dos('Si(s)', pdos_Si[OrbitalType.s])
        # plotter.add_dos('Si(p)', pdos_Si[OrbitalType.p])
        plotter.get_plot(xlim=(-13, 15), ylim=(0, 3))
        plotter.save_plot(f"{material_path}/dos.png")
    
    def plot_band_structure(self, material_path: str = ''):
        if not material_path:
            raise ValueError("material_path must be provided.")
        vaspout = Vasprun(f"{material_path}/3_bs/vasprun.xml")
        band_str = vaspout.get_band_structure(line_mode=True)
        plt = BSPlotter(band_str)
        plt.get_plot(ylim=[-12,10])
        plt.save_plot(f"{material_path}/band_structure.png")
        
        


def scf_ready(material: str):
    experiment = Experiment(material_dir=material)
    experiment.ensure_required_files()
    experiment.scf_calculation()

def dos_ready(material: str):
    experiment = Experiment(material_dir=material)
    experiment.dos()
    # experiment.plot_dos(material_path=material)

def bs_ready(material: str):
    structure = Structure.from_file(f'{material}/1_scf/POSCAR')
    kpath = HighSymmKpath(structure)
    kpts = Kpoints.automatic_linemode(divisions=40, ibz=kpath)
    kpts.write_file(f"{material}/3_bs/KPOINTS")

    experiment = Experiment(material_dir=material)
    experiment.band_structure()
    # experiment.plot_band_structure(material_path=material)