import numpy as np

from results import *
import qe.runner as runner
from config import *


def _write_input(input_path, input_data_path):
    
    with open(input_path, "w") as input_file:
        input_file.write(f"1\n"
                         f"{input_data_path}\n"
                         f"1.0\n"
                         f"1000\n"
                         f"3\n"
                         f"5.0\n")


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
    
    z, planar, macroscopic = _calculate(path, input_data_path)
    
    return PotentialAveraged(coordinates=z, 
                             planar=planar, 
                             macroscopic=macroscopic)


def charge_density(path, input_data_path):
    
    z, planar, macroscopic = _calculate(path, input_data_path)

    return ChargeDensityAveraged(coordinates=z, 
                                 planar=planar,
                                 macroscopic=macroscopic)


def ldos(path, input_data_path):
    
    averaged_dos = []
    coordinates = None
        
    for file in sorted(input_data_path.parent.glob(f"{input_data_path.name}*")):
        print(file)
        
        z, planar, _ = _calculate(path, file)
        
        if coordinates is None:
            coordinates = z
        
        else:
            np.testing.assert_allclose(coordinates, z)
        
        averaged_dos.append(planar)
    
    averaged_dos = np.asarray(averaged_dos)
        
    energies = np.arange(WINDOW_LDOS[0], WINDOW_LDOS[1] + DELTA_E, DELTA_E)
    
    return LDOSAveraged(energies=energies, coordinates=coordinates, planar=averaged_dos)
