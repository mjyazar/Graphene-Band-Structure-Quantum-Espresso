import numpy as np

from results import *
import qe.runner as runner
from config import *


def _write_input(input_path, input_data_path):
    
    with open(input_path, "w") as input_file:
        input_file.write(f"""
        1
        "{input_data_path}"
        1.0
        1000
        3
        5
        """)


def _read_output(path):
    
    coordinates, planar_average, macroscopic_average = np.loadtxt(path, unpack=True)
    
    return coordinates, planar_average, macroscopic_average


def _calculate(path, input_data_path):
    
    path.mkdir(parents=True, exist_ok=True)

    process = input_data_path.name.removesuffix(".pp.dat")
    
    input_path = path / f"{process}.avg.in"
    log_path = path / f"{process}.avg.log"  # log file
    output_path = path / f"{process}.avg.dat"
    
    print(f"\nCREATING {input_path.name}")
    _write_input(input_path, input_data_path)

    print(f"RUNNING average.x WITH {input_path.name}")
    runner.run("average.x", input_path, log_path, cwd=path, nproc=1)

    qe_output_path = path / "avg.dat"
    qe_output_path.replace(output_path)
    
    print(f"READING {output_path.name}")
    return _read_output(output_path)


def potential(path, input_data_path):
    
    coordinates, planar_average, macroscopic_average = _calculate(path, input_data_path)
    
    return PotentialAveraged(coordinates=coordinates, 
                             planar=planar_average, 
                             macroscopic=macroscopic_average)
