import numpy as np

import qe.runner as runner

DEGAUSS = 0.01
WINDOW_DOS = (-10, 10)

def write_input(input_path, outdir, data_path, prefix, fermi_energy, window=WINDOW_DOS):
    
    with open(input_path, "w") as input_file:
        input_file.write(f"""&DOS  ! QE input begins
                         prefix = "{prefix}"
                         outdir = "{outdir}"  ! directory containing the input data, i.e. the pw.x metadata
                         bz_sum = "smearing"  ! integration using gaussian smearing
                         ngauss = 0  ! type of gaussian broadening - 0: Simple Gaussian (default)
                         degauss = {DEGAUSS}  ! gaussian broadening, Ry (not eV!)
                         emin = {fermi_energy + window[0]}
                         emax = {fermi_energy + window[1]}
                         deltaE = 0.02  ! energy grid step (eV)
                         fildos = "{data_path}"  ! output file containing DOS(E)
                         /
                         """)


def read_output(path):
    
    # idos - integrated dos
    energy, dos, idos = np.loadtxt(path, unpack=True)
    
    return energy, dos, idos


def calculate(path, fermi_energy):
    
    path.mkdir(parents=True, exist_ok=True)
    
    input_path = path / "dos.in"
    output_path = path / "dos.out"  # log file
    data_path = path / "dos.data"  # data file
    
    print(f"\nCREATING {input_path.name}")
    write_input(input_path, path / "data", data_path, path.name, fermi_energy)

    print(f"RUNNING dos.x WITH {input_path.name}")
    runner.run("dos.x", input_path, output_path)
    
    print(f"READING {data_path.name}")
    
    return read_output(data_path)
