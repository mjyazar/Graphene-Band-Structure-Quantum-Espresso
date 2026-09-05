import numpy as np

import qe.runner as runner


def write_input(input_path, outdir, data_path, prefix):
    
    with open(input_path, "w") as input_file:
        input_file.write(f"""&DOS  ! QE input begins
                         prefix = '{prefix}'
                         outdir = '{outdir}'  ! directory containing the input data, i.e. the pw.x metadata
                         bz_sum = 'smearing'  ! integration using gaussian smearing
                         ngauss = 0  ! type of gaussian broadening - 0: Simple Gaussian (default)
                         degauss = 0.005  ! gaussian broadening, Ry (not eV!)
                         deltae = 0.01  ! energy grid step (eV)
                         fildos = '{data_path}'  ! output file containing DOS(E)
                         /
                         """)


def read_output(path):
    
    # idos - integrated dos
    energy, dos, idos = np.loadtxt(path, unpack=True)
    
    return energy, dos, idos


def calculate(path):
    
    path.mkdir(parents=True, exist_ok=True)
    
    dos_path = path / "dos"
    dos_path.mkdir(parents=True, exist_ok=True)
    
    input_path = dos_path / "dos.in"
    output_path = dos_path / "dos.out"  # log file
    data_path = dos_path / "dos.data"  # data file
    
    print(f"\nCREATING {input_path.name}")
    write_input(input_path, path / "data", data_path, path.name)

    print(f"RUNNING dos.x WITH {input_path.name}")
    runner.run("dos.x", input_path, output_path)
    
    print(f"\nREADING {data_path.name}")
    
    return read_output(data_path)
