import numpy as np

import qe.runner as runner

DEGAUSS = 0.01
WINDOW_DOS = (-10, 10)
DELTA_E = 0.02

# RY_TO_EV = 

def _write_input(input_path, outdir, output_path, prefix, fermi_energy, window=WINDOW_DOS):
    
    with open(input_path, "w") as input_file:
        input_file.write(f"""&DOS
                         prefix = "{prefix}"
                         outdir = "{outdir}"  ! directory containing the input data, i.e. the pw.x metadata
                         bz_sum = "smearing"  ! integration using gaussian smearing
                         ngauss = 0  ! type of gaussian broadening - 0: Simple Gaussian (default)
                         degauss = {DEGAUSS}  ! gaussian broadening, Ry (not eV!)
                         emin = {fermi_energy + window[0]}
                         emax = {fermi_energy + window[1]}
                         deltaE = {DELTA_E}  ! energy grid step (eV)
                         fildos = "{output_path}"  ! output file containing DOS(E)
                         /
                         """)


def _read_output(path):
    
    # idos - integrated dos
    energy, dos, idos = np.loadtxt(path, unpack=True)
    
    return energy, dos, idos


def calculate(path, fermi_energy):
    
    path.mkdir(parents=True, exist_ok=True)
    
    outdir = path / "data"
    input_path = path / "dos.in"
    log_path = path / "dos.log"  # log file
    output_path = path / "dos.out"  # data file
    
    print(f"\nCREATING {input_path.name}")
    _write_input(input_path, outdir, output_path, path.name, fermi_energy)

    print(f"RUNNING dos.x WITH {input_path.name}")
    runner.run("dos.x", input_path, log_path)
    
    print(f"READING {output_path.name}")
    return _read_output(output_path)
