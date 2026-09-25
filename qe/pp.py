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
filplot = "{intermediate_path}"  ! QE intermediate 3D FFT i.e. real-space grid data
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
        # LDOS is plotted on grid [emin, emax] with spacing delta_e_ldos
        elif computation == 3:
            input_file.write(f"""
emin = {fermi_energy + window[0]}
emax = {fermi_energy + window[1]}
delta_e = {DELTA_E_LDOS}
degauss_ldos = {DEGAUSS_LDOS}
use_gauss_ldos = .true.
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
        
        xy_averaged = np.mean(data, axis=(0, 1))
        z = _coordinates(data, atoms)
        
        return data, atoms, xy_averaged, z
    
    elif output_format == 7:
        x, y, z = np.loadtxt(path, unpack=True)
        
        return x, y, z
    
    else:
        raise ValueError("Invalid output_format.")
    

def _read_ldos(intermediate_path):

    averaged_ldos = []
    z = None
    atoms_ = None
    
    files = sorted(intermediate_path.parent.glob(f"{intermediate_path.name}[0-9]*.cube"))
    file_count = len(files)
    
    for i, file in enumerate(files, start=1):
        print(f"\rREADING {file.name} [{i}/{file_count}]", flush=True, end="")
        
        data, atoms = read_cube_data(file)

        if i == 1:
            z = _coordinates(data, atoms)
            atoms_ = atoms
            
        averaged_ldos.append(np.mean(data, axis=(0, 1)))
    
    averaged_ldos = np.asarray(averaged_ldos)
    
    # shape[0] gets the number of rows i.e. the energy count
    energies = WINDOW_LDOS[0] + np.arange(averaged_ldos.shape[0]) * DELTA_E_LDOS
    
    return LDOS(averaged=averaged_ldos, z=z, energies=energies, atoms=atoms_)


def _calculate(computation, path, fileout, iflag=3, fermi_energy=None):
    
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
        if fileout == "ldos":
            output_path = ".cube"
        
        else:
            output_path = path / f"{fileout}.cube"
        
        # output_format = 6  # 6 = gaussian cube file (3D)
        output_format = 5  # 5 = XCRYSDEN  (3D, using entire FFT grid)
    
    else:
        raise ValueError("iflag MUST BE A VALID VALUE")
    
    # delete old files if running new ldos calculations - prevents issues when changing parameters
    if fileout == "ldos":
        for old_file in intermediate_path.parent.glob(f"{intermediate_path.name}[0-9]*"):
            old_file.unlink()
    
    print(f"\nCREATING {input_path.name}")
    _write_input(computation, input_path, outdir, intermediate_path, output_path, path.name, fermi_energy, iflag, output_format)

    print(f"RUNNING pp.x WITH {input_path.name}")
    runner.run("pp.x", input_path, log_path)        

    if fileout == "ldos":
        return _read_ldos(intermediate_path), intermediate_path
    
    print(f"READING {output_path.name}")
    return _read_output(output_path, output_format), intermediate_path


def _coordinates(data, atoms):
    number_of_z_points = data.shape[2]  # number of z grid points
    total_length = atoms.cell.lengths()[2]

    return np.linspace(0, total_length, number_of_z_points, endpoint=False)


def potential(path):
    
    (data, atoms, xy_averaged, z), intermediate_path = _calculate(11, path, "potential")

    results = Potential(data=data, atoms=atoms, averaged=xy_averaged, z=z)
    
    return results, intermediate_path


def charge_density(path):
    (data, atoms, xy_averaged, z), intermediate_path = _calculate(0, path, "charge",)
    
    results = ChargeDensity(data=data, atoms=atoms, averaged=xy_averaged, z=z)
    
    return results, intermediate_path


def ldos(path, fermi_energy):
    ldos, intermediate_path = _calculate(3, path, "ldos", fermi_energy=fermi_energy)
    
    return ldos, intermediate_path, ldos.energies
