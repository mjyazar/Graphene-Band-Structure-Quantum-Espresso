import numpy as np
from ase.io.cube import read_cube_data

from results import *
import qe.runner as runner
from config import *


def _write_input(computation, input_path, outdir, intermediate_path, fileout, prefix, fermi_energy, iflag, output_format, window=WINDOW_LDOS):
    
    with open(input_path, "w") as input_file:
        input_file.write(f"""&INPUTPP
                         prefix = "{prefix}"
                         outdir = "{outdir}"  ! directory containing the input data, i.e. the pw.x metadata
                         plot_num = {computation}
                         filplot = "{intermediate_path}"  ! intermediate metadata path
                         """)
        
        # 0  = electron (pseudo-)charge density
        # total charge
        if computation == 0:
            input_file.write(f"""
                             spin_component = 0  ! total charge
                             """)
        
        # 1 = total potential V_bare + V_H + V_xc
        # total potential including xc
        elif computation == 1:
            input_file.write(f"""
                             spin_component = 0  ! spin averaged potential
                             """)
          
        # 3 = local density of states at specific energy or grid of energies
        # (number of states per volume, in bohr^3, per energy unit, in Ry)
        # LDOS is plotted on grid [emin, emax] with spacing delta_e. 
        elif computation == 3:
            input_file.write(f"""
                             emin = {fermi_energy + window[0]}
                             emax = {fermi_energy + window[1]}
                             delta_e = {DELTA_E}
                             degauss_ldos = {DEGAUSS_LDOS}
                             """)       
        
        # 10 = integrated local density of states (ILDOS) from emin to emax 
        # (emin, emax in eV) if emax is not specified, emax=E_fermi
        elif computation == 10:
            input_file.write(f"""
                             emin = {fermi_energy + window[0]}
                             emax = {fermi_energy + window[1]}
                             spin_component = 0  ! spin-up + spin-down
                             """)
        
        input_file.write("/\n")
        
        input_file.write(f"""&PLOT
                         filepp(1) = "{intermediate_path}"
                         iflag = {iflag}
                         output_format = {output_format}
                         fileout = "{fileout}"
                         """)
        
        if iflag == 2:
            input_file.write("""
                             e1(1) = 1.0
                             e1(2) = 0.0
                             e1(3) = 0.0

                             e2(1) = 0.0
                             e2(2) = 0.0
                             e2(3) = 9.491870

                             x0(1) = 0.0
                             x0(2) = 0.2886752
                             x0(3) = 0.0

                             nx = 30
                             ny = 288
            """)
        
        input_file.write("/\n")


def _read_output(path, output_format):

    if output_format == 6:
        # data has shape (nx, ny, nz)
        data, atoms = read_cube_data(path)
        
        return data, atoms
    
    elif output_format == 7:
        x, y, z = np.loadtxt(path, unpack=True)
        
        return x, y, z
    
    else:
        raise ValueError("Invalid output_format.")
    

def _calculate(computation, path, iflag, fileout, fermi_energy=None):
    
    path.mkdir(parents=True, exist_ok=True)
    
    outdir = path / "data"
    input_path = path / f"{fileout}.pp.in"
    log_path = path / f"{fileout}.pp.log"  # log file
    intermediate_path = path / "data" / f"{fileout}.pp.dat"  # intermediate metadata
    
    # 2D plot
    if iflag == 2:
        output_path = path / f"{fileout}.dat"  # output data file
        output_format = 7  # 7 = format suitable for gnuplot (2D): x, y, f(x,y)
    
    # 3D plot
    elif iflag == 3:
        output_path = path / f"{fileout}.cube"
        output_format = 6  # 6  = format as gaussian cube file (3D)
    
    print(f"\nCREATING {input_path.name}")
    _write_input(computation, input_path, outdir, intermediate_path, output_path, path.name, fermi_energy, iflag, output_format)

    print(f"RUNNING pp.x WITH {input_path.name}")
    runner.run("pp.x", input_path, log_path)
    
    print(f"READING {output_path.name}")
    return _read_output(output_path, output_format), intermediate_path


def potential(path):
    
    (data, atoms), intermediate_path = _calculate(11, path, 3, "potential")
    
    results = Potential(data=data, atoms=atoms)
    
    return results, intermediate_path


def charge_density(path):
    (data, atoms), intermediate_path = _calculate(0, path, 3, "charge",)
    
    results = ChargeDensity(data=data, atoms=atoms)
    
    return results, intermediate_path


def ldos(path, fermi_energy):

    # (data, atoms), intermediate_path = _calculate(3, path, 3, "ldos", fermi_energy)
    # results = LDOS(data=data, atoms=atoms)
    # return results, intermediate_path

    path.mkdir(parents=True, exist_ok=True)
    
    outdir = path / "data"
    input_path = path / "ldos.pp.in"
    log_path = path / "ldos.pp.log"  # log file
    intermediate_path = path / "data" / "ldos.pp.dat"  # intermediate metadata
    
    # output_path = path / f"ldos.cube"
    # output_format = 6  # 6  = format as gaussian cube file (3D)

    print(f"\nCREATING {input_path.name}")
    with open(input_path, "w") as input_file:
        input_file.write(f"""&INPUTPP
                         prefix = "{path.name}"
                         outdir = "{outdir}"  ! directory containing the input data, i.e. the pw.x metadata
                         plot_num = 3
                         filplot = "{intermediate_path}"  ! intermediate metadata path
                         
                         emin = {fermi_energy + WINDOW_LDOS[0]}
                         emax = {fermi_energy + WINDOW_LDOS[1]}
                         delta_e = {DELTA_E}
                         degauss_ldos = {DEGAUSS_LDOS}
                         /
                         """)
    
    print(f"RUNNING pp.x WITH {input_path.name}")
    runner.run("pp.x", input_path, log_path)
    
    return intermediate_path
