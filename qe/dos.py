from pathlib import Path
import numpy as np

import qe.runner as runner


ROOT = Path(__file__).resolve().parent


def write_input(path, outdir, fildos, prefix="graphene"):
    
    with open(path, "w") as input_file:
        input_file.write(f"""&DOS  # QE input begins
                         prefix = '{prefix}'
                         outdir = '{outdir}'  # directory containing the input data, i.e. the same as in pw.x
                         bz_sum = 'smearing'  # integration using gaussian smearing
                         ngauss = 0  # type of gaussian broadening - 0: Simple Gaussian (default)
                         degauss = 0.01  # gaussian broadening, Ry (not eV!)
                         deltae = 0.01  # energy grid step (eV)
                         fildos = '{fildos}'  # output file containing DOS(E)
                         /
                         """)


def read_output(path):
    
    data = np.loadtxt(path)
    
    energy = data[:, 0]
    dos = data[:, 1]
    
    return energy, dos


def calculate(path, input_data_path):

    path.mkdir(parents=True, exist_ok=True)
    
    input_path = path / f"dos.in"
    output_path = path / f"dos.out"
    fildos = path / "dos.dat"
    
    print(f"\nCREATING {input_path.name}")
    write_input(input_path, input_data_path, fildos)

    print(f"RUNNING dos.x WITH {input_path.name}")
    runner.run(input_path, output_path, "dos.x")
    
    print(f"\nCREATING {input_path.name}")
    energy, dos = read_output(fildos)
    
    return energy, dos
