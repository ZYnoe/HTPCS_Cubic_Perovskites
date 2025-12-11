INCAR_SCF = """PREC=Accurate
ENCUT=500
NELMIN=10
NELM=100          # Increase max electronic steps for better convergence
EDIFF=1E-6        # Tighten electronic energy convergence criterion
NSW=0
IBRION=-1
ISMEAR=0
SIGMA=0.01
ADDGRID=TRUE
LWAVE=TRUE
LCHARG=TRUE
LREAL=FALSE
LASPH=TRUE
LORBIT=11          # Use 11 to get detailed DOS and projected info

KGAMMA=TRUE
KSPACING=0.15

NPAR=32
"""

INCAR_DOS = """PREC=Accurate
ENCUT=500
NELMIN=10
NELM=100           # Same as SCF to ensure good electronic convergence
EDIFF=1E-6         # Keep tighter convergence
NSW=0
IBRION=-1
ISMEAR=-5          # Tetrahedron method with Blöchl corrections for DOS
SIGMA=0.01
ADDGRID=TRUE
LCHARG=FALSE
LWAVE=FALSE
LREAL=FALSE
LASPH=TRUE
LORBIT=11

NEDOS=2001
EMAX=20.00
EMIN=-20.00

KGAMMA=TRUE
KSPACING=0.15

NPAR=32
"""

INCAR_BS = """PREC=Accurate
ENCUT=500
NELMIN=10
NELM=100           # Increase to maintain SCF quality
EDIFF=1E-6
NSW=0
IBRION=-1
ISMEAR=0
SIGMA=0.01
ADDGRID=TRUE
LWAVE=FALSE
LCHARG=FALSE
LREAL=FALSE
LASPH=TRUE
LORBIT=11

ICHARG=11          # Use charge density from previous SCF for band structure calc

KGAMMA=TRUE
KSPACING=0.15

NPAR=32

"""

job_sh_template = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=g1,g2
##
#SBATCH --job-name=zylu_demo
#SBATCH --time=99-24:59          # Runtime: Day-HH:MM
#SBATCH -o log.txt         # STDOUT
#SBATCH -e error.txt         # STDERR
#################################################################################

#################################################################################
VASP_BIN=/TGM/Apps/VASP/VASP_BIN/6.3.2/vasp.6.3.2.std.x
#################################################################################

hostname
date

#################################################################################
### Module Env
module purge
module add compiler/2022.1.0
module add mkl/2022.1.0
module add mpi/2021.6.0
export OMP_NUM_THREADS=1
export RSH_COMMAND="/usr/bin/ssh -x"
scontrol show hostname ${SLURM_JOB_NODELIST} | perl -ne 'chomb; print "$_"x1'> ${SLURM_SUBMIT_DIR}/nodelist
#################################################################################

mpirun -np $SLURM_NTASKS $VASP_BIN > stdout
"""

KPOINTS = """Automatic mesh
0
Auto
  10

"""

import argparse
import multiprocessing as mp
import os
import sys
from pathlib import Path
import shutil
from fireworks.core.firework import FWorker
from mp_api.client import MPRester
from emmet.core.symmetry import CrystalSystem
from fireworks import Firework, LaunchPad, FileTransferTask, PyTask, Workflow, ScriptTask
from fireworks.core.rocket_launcher import rapidfire
import Auto_baby

SCRIPT_DIR = Path(__file__).resolve().parent
_REMOVED_SYS_PATH_ENTRIES = []
_SCRIPT_DIR_STR = str(SCRIPT_DIR)
for idx in reversed(range(len(sys.path))):
    entry = sys.path[idx]
    if entry in ("", _SCRIPT_DIR_STR):
        _REMOVED_SYS_PATH_ENTRIES.append((idx, entry))
        sys.path.pop(idx)

def generate_potcar_from_poscar_flat(base_path):
    """
    Concatenates and generates a VASP POTCAR file from a flattened directory structure, based on the list of elements in the POSCAR file.

    Args:

    base_path (str): The root directory for the material calculation. The function assumes the POSCAR is located at base_path/1_scf/POSCAR.

    Returns:

    None. Prints status messages to the console.
    """
    # --- Path Definitions ---
    # The root directory of the POTCAR library, now points directly to your potcars folder.

    pot_path = Path('/home/student/zylu/Lu_potpaw_copy')
    base_path = Path(base_path)
    poscar_path = base_path / '1_scf' / 'POSCAR'
    output_path = base_path / '1_scf' / 'POTCAR'
    
    print(f"Generating POTCAR from {poscar_path} to {output_path}...")

    # --- Reading element symbols from POSCAR file ---
    try:
        with open(poscar_path, 'r') as f:
            lines = f.readlines()
        # In VASP 5+ format, element symbols are on the 6th line (index 5)
        elements = lines[5].strip().split()
        if not elements:
            print(f"Error: No element symbols found on line 6 of POSCAR file '{poscar_path}'.")
            return
        print(f"Elements found in POSCAR: {elements}")
    except FileNotFoundError:
        print(f"Error: File '{poscar_path}' not found.")
        return
    except IndexError:
        print(f"Error: Unable to read element line from '{poscar_path}'. It may not be a valid VASP 5+ format.")
        return

    # --- Concatenating POTCAR files ---
    try:
        with open(output_path, 'wb') as potcar_file:
            for element in elements:
                # *** Main modification part ***
                # Construct the path according to your 'POTCAR' filename format
                element_pot_path = pot_path / element / "POTCAR"

                if not os.path.exists(element_pot_path):
                    print(f"Error: POTCAR file for element '{element}' not found at '{element_pot_path}'")
                    # Clean up partially created file
                    potcar_file.close()
                    os.remove(output_path)
                    return

                with open(element_pot_path, 'rb') as element_file:
                    potcar_file.write(element_file.read())
        
        print(f"Successfully generated POTCAR for elements: {' '.join(elements)} at '{output_path}'")

    except Exception as e:
        (f"An unexpected error occurred during file writing: {e}")

def create_vasp_project(material_id, material_formula, base_path, POSCAR_content):

    base_dir = os.path.join(base_path, f"{material_id}_{material_formula}")
    scf_dir = os.path.join(base_dir, "1_scf")
    dos_dir = os.path.join(base_dir, "2_dos")
    bs_dir = os.path.join(base_dir, "3_bs") 

    # Create directories
    for d in [base_dir, scf_dir, dos_dir, bs_dir]:
        os.makedirs(d, exist_ok=True)

    # Create POSCAR, INCAR, job.sh in 1_scf directory
    # with open(os.path.join(scf_dir, "POSCAR"), "w") as f:
    #     f.write(f"{POSCAR_content}")
    target_poscar = os.path.join(scf_dir, "POSCAR")
    
    if os.path.exists(POSCAR_content) and os.path.isfile(POSCAR_content):
        # It's a file path, copy it
        shutil.copy(POSCAR_content, target_poscar)
    else:
        # It's a string content, write it
        with open(target_poscar, "w") as f:
            f.write(POSCAR_content)

    with open(os.path.join(scf_dir, "INCAR"), "w") as f:
        f.write(INCAR_SCF)
    with open(os.path.join(scf_dir, "job.sh"), "w") as f:
        f.write(job_sh_template)

    # Create INCAR in 2_dos 
    with open(os.path.join(dos_dir, "INCAR"), "w") as f:
        f.write(INCAR_DOS)

    # Create INCAR in 3_bs 
    with open(os.path.join(bs_dir, "INCAR"), "w") as f:
        f.write(INCAR_BS)
    with open(os.path.join(bs_dir, "KPOINTS"),"w") as f:
        f.write(KPOINTS)
    
    generate_potcar_from_poscar_flat(base_path = base_dir)
    print(f"Project structure has been created for {material_id}_{material_formula}.")

def fetch_cubic_docs():
    with MPRester(API_KEY) as mpr:
        results = mpr.materials.summary.search(
            formula=["ABC3"],
            crystal_system=CrystalSystem.cubic,
        )
    docs = list(results)
    print(f"found {len(docs)} ABX3 cubic materials.")
    return docs

for idx, entry in sorted(_REMOVED_SYS_PATH_ENTRIES):
    sys.path.insert(idx, entry)


API_KEY = "5ypNyStQsaDUA2neNLD17P8v2KNB241U"
MATERIALS_BASE_PATH = "/home/student/zylu/materials_5"

def Final_Calculation(launchpad, material_id, formula, material_path: str = ''):
    if not material_path: raise ValueError("material_path must be provided.")
    

    _submit_SCF_Task_ = PyTask(func = "Auto_baby.scf_ready",
                               args=[material_path])


    

    _transfer_POSCAR_of_SCF_to_DOS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/POSCAR', 
                                                                   'dest': f'{material_path}/2_dos'}], 'mode': 'copy'})
    _transfer_POTCAR_of_SCF_to_DOS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/POTCAR', 
                                                       'dest': f'{material_path}/2_dos'}], 'mode': 'copy'})
    _transfer_JOB_of_SCF_to_DOS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/job.sh', 
                                                  'dest': f'{material_path}/2_dos'}], 'mode': 'copy'})
    _transfer_WAVECAR_of_SCF_to_DOS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/WAVECAR', 
                                                       'dest': f'{material_path}/2_dos'}], 'mode': 'copy'})
    _transfer_CONTCAR_of_SCF_to_DOS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/CHGCAR', 
                                                                   'dest': f'{material_path}/2_dos/CHGCAR'}], 'mode': 'copy'})

    _submit_DOS_Task_ = PyTask(func= "Auto_baby.dos_ready",
                               args=[material_path])



    _transfer_CONTCAR_of_SCF_to_BS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/CHGCAR', 
                                                    'dest': f'{material_path}/3_bs/POSCAR'}], 'mode': 'copy'})
    _transfer_POTCAR_of_SCF_to_BS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/POTCAR', 
                                                    'dest': f'{material_path}/3_bs'}], 'mode': 'copy'})
    _transfer_JOB_of_SCF_to_BS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/job.sh', 
                                                    'dest': f'{material_path}/3_bs'}], 'mode': 'copy'})
    _transfer_CHGCAR_of_SCF_to_BS = FileTransferTask({'files': [{'src': f'{material_path}/1_scf/CHGCAR', 
                                                    'dest': f'{material_path}/3_bs'}], 'mode': 'copy'})
    

    _submit_BS_Task_ = PyTask(func = "Auto_baby.bs_ready",
                             args=[material_path])

    scf_calculation = Firework(
        tasks=[_submit_SCF_Task_],
        name="SCF Calculation" 
    )

    DOS_calculation = Firework(
        tasks=[_transfer_POSCAR_of_SCF_to_DOS, _transfer_CONTCAR_of_SCF_to_DOS, _transfer_POTCAR_of_SCF_to_DOS, _transfer_JOB_of_SCF_to_DOS, _transfer_WAVECAR_of_SCF_to_DOS,_submit_DOS_Task_],
        name="DOS Calculation",
        parents=[scf_calculation]
    )

    BAND_structure_calculation = Firework(
        tasks=[_transfer_CONTCAR_of_SCF_to_BS, _transfer_POTCAR_of_SCF_to_BS, _transfer_JOB_of_SCF_to_BS, _transfer_CHGCAR_of_SCF_to_BS, _submit_BS_Task_],
        name="Band Structure Calculation",
        parents=[scf_calculation]
    )



    wf = Workflow([scf_calculation, DOS_calculation, BAND_structure_calculation], name = f"{material_id}_{formula}")
    launchpad.add_wf(wf)

def queue_workflows(docs, launchpad, base_path: str = MATERIALS_BASE_PATH):
    for doc in docs:
        structure = doc.structure
        formula = doc.formula_pretty
        material_id = doc.material_id

        print(f"Processing material: {material_id}, formula: {formula}, material_id: {material_id}")

        create_vasp_project(
            material_id,
            formula,
            base_path=base_path,
            POSCAR_content=structure.to(fmt="poscar")
        )

        material_path = os.path.join(base_path, f"{material_id}_{formula}")
        Final_Calculation(
            launchpad=launchpad,
            material_path=material_path,
            material_id=material_id,
            formula=formula
        )

def _run_single_worker(worker_idx: int):
    """Spawn a dedicated FireWorks rocket that keeps pulling ready jobs."""
    worker_launchpad = LaunchPad()
    fworker = FWorker(name=f"parallel-worker-{worker_idx}")
    rapidfire(worker_launchpad, fworker=fworker, nlaunches=0)

def run_workflows_parallel(num_workers: int):
    if num_workers < 1:
        raise ValueError("num_workers must be >= 1")

    processes = []
    try:
        for idx in range(num_workers):
            process = mp.Process(target=_run_single_worker, args=(idx,))
            process.start()
            processes.append(process)

        for process in processes:
            process.join()
    finally:
        for process in processes:
            if process.is_alive():
                process.terminate()

def parse_args():
    parser = argparse.ArgumentParser(description="Submit VASP workflows and run them in parallel via FireWorks.")
    default_workers = max(1, (os.cpu_count() or 1))
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=default_workers,
        help=f"Number of parallel rockets to spawn (default: {default_workers})."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)

    docs = fetch_cubic_docs()
    if not docs:
        print("No matching materials were found. Exiting without launching workflows.")
        return

    queue_workflows(docs, launchpad)
    run_workflows_parallel(args.workers)

if __name__ == "__main__":
    main()
