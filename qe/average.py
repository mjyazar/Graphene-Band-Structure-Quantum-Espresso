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


def _calculate(path, input_data_path, output_path = None, verbose=True):
    
    path.mkdir(parents=True, exist_ok=True)

    process = input_data_path.name.split(".pp.dat")[0]
    
    input_path = path / f"{process}.avg.in"
    log_path = path / f"{process}.avg.log"  # log file
    
    if output_path is None:
        output_path = path / f"{process}.avg.dat"
    
    if verbose:
        print(f"\nCREATING {input_path.name}")
    _write_input(input_path, input_data_path)

    if verbose:
        print(f"RUNNING average.x WITH {input_path.name}")
    runner.run("average.x", input_path, log_path, cwd=path, nproc=1)

    qe_output_path = path / "avg.dat"
    qe_output_path.replace(output_path)
    
    if verbose:
        print(f"READING {output_path.name}")
    return _read_output(output_path)


def potential(path, input_data_path):
    
    z, planar, _ = _calculate(path, input_data_path)
    
    return PotentialAveraged(coordinates=z, 
                             planar=planar)


def charge_density(path, input_data_path):
    
    z, planar, _ = _calculate(path, input_data_path)

    return ChargeDensityAveraged(coordinates=z, 
                                 planar=planar)


def ldos(path, input_data_path, energies):
    
    averaged_ldos = []
    coordinates = None
        
    files = sorted(input_data_path.parent.glob(f"{input_data_path.name}[0-9]*"), 
                   key=lambda f: int(f.name.split("dat")[-1]))
    file_count = len(files)
    
    print("\nCREATING ldos.avg.in")

    for i, file in enumerate(files, start=1):
        
        file_number = file.name.split("dat")[-1]
        averaged_path = path / "data" / f"ldos.avg.dat{file_number}"
        
        # Read averaged files if already averaged
        if averaged_path.exists():
            z, planar, _ = _read_output(averaged_path)
        
        # otherwise run average.x
        else:
            print(f"\rRUNNING average.x WITH ldos.avg.in [{i}/{file_count}]", flush=True)
        
            z, planar, _ = _calculate(path, file, output_path=averaged_path,  verbose=False)

            print("\033[2A", end="")  # move back up to the LDOS line
        
        if coordinates is None:
            coordinates = z
    
        else:
            np.testing.assert_allclose(coordinates, z)
        
        averaged_ldos.append(planar)

    print("\nREADING ldos.avg.dat")
    averaged_dos = np.asarray(averaged_ldos)
    
    return LDOSAveraged(energies=energies, coordinates=coordinates, planar=averaged_dos)
