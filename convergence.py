import qe.pw as pw
import numpy as np
from config import *

from pathlib import Path
import matplotlib.pyplot as plt
from ase.dft.dos import DOS

ROOT = Path(__file__).resolve().parent

FIG_DIR = ROOT / "convergence" / "figures"
DATA_DIR = ROOT / "convergence" / "data"

FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)


def _graph_convergence(values, energies, field, header, parameter):
    
    fig, ax = plt.subplots()
    
    ax.plot(values, energies, "o-", color="black")
    ax.set_title(f"{parameter} Parameter Convergence at E-field={field}au")
    
    if parameter.startswith("kgrid"):
        ax.set_xlabel(f"{parameter} Parameter x -> (x, x, 1)")
    
    else:
        ax.set_xlabel(f"{parameter} Parameter")

    ax.set_ylabel(header)
    
    fig.savefig(FIG_DIR / f"{parameter}_Convergence_{field}au", bbox_inches="tight")
    plt.close(fig)


def test_kgrid(structure, field, computation):
    """
    Convergence testing of kgrid by iterating through from values of x and y 
    ranging from 3 to upper bound, with z being kept constant at 1.
    Involves running the scf process and plotting the total energy against kgrid values.
    return: scf.pwo file read using ASE
    """
    
    # list to store values for plotting
    kgrid_values = []
    energy_values = []
    doses = []
    
    if computation == "scf":
        calc_path = DATA_DIR / "scf"
        l_bound = 3
        u_bound = SCF_BOUND
        label = "scf"

    elif computation == "nscf":
        calc_path = DATA_DIR / "nscf"
        l_bound = NSCF_CONVERGENCE_KGRID
        u_bound = NSCF_BOUND
        label = "nscf"
        
        scf = pw._calculate("scf", structure, calc_path, (NSCF_CONVERGENCE_KGRID, NSCF_CONVERGENCE_KGRID, 1), efield=field)
    
    for i in range(l_bound, u_bound+1, 3):
        kgrid = (i, i, 1)
                
        if computation == "scf":
            scf = pw._calculate("scf", structure, calc_path, kgrid, efield=field)
            total_energy = scf.get_potential_energy()  # extract total energy
            energy_values.append(total_energy)
            print(f"KGRID = ({i}, {i}, 1) --> TOTAL ENERGY: {total_energy}")
            
        elif computation == "nscf":
            nscf = pw._calculate("nscf", structure, calc_path, kgrid, efield=field)
            fermi_energy = nscf.calc.get_fermi_level()
            energy_values.append(fermi_energy)
            print(f"KGRID = ({i}, {i}, 1) --> FERMI ENERGY: {fermi_energy}")
            doses.append(DOS(nscf.calc, width=0.1, window=(-10, 10), npts=1000).get_dos())
        
        kgrid_values.append(i)
    
    kgrid_values = np.array(kgrid_values)
    energy_values = np.array(energy_values)
    
    # save the results as a .txt file for reference
    header=f"kgrid {'total_energy_eV' if computation == 'scf' else 'fermi_energy_eV'}"
    np.savetxt(DATA_DIR / f"{label} kgrid Convergence.txt", np.column_stack((kgrid_values, energy_values)), header=header)
    
    _graph_convergence(kgrid_values, energy_values, field, header, f"kgrid_{label}")
    
    if computation == "nscf":
        error = [np.abs(d - doses[-1]).sum() / doses[-1].sum() for d in doses[:-1]]
        _graph_convergence(kgrid_values[:-1], error, field, "Relative DOS error", "kgrid_nscf_dos")
        
    return (kgrid_values, energy_values)


def test_ecutwfc(structure, field):
    """
    Convergence testing of ecutwfc by iterating through values 30 to upper bound in increments of 10.
    Involves running the scf process and plotting the total energy against ecutwfc values.
    return: scf.pwo file read using ASE
    """
    
    # list to store values for plotting
    ecutwfc_values = []
    energy_values = []
    
    for ecutwfc in range(30, ECUTWFC_BOUND+1, 10):
        
        calculation_path = DATA_DIR / "ecutwfc" / str(ecutwfc)
        
        scf = pw._calculate("scf", structure, calculation_path, KGRID, ecutwfc=ecutwfc, ecutrho=8*ecutwfc, efield=field)
        
        total_energy = scf.get_potential_energy()  # extract total energy
    
        ecutwfc_values.append(ecutwfc)
        energy_values.append(total_energy) 
        print(f"ECUTWFC = {ecutwfc} --> TOTAL ENERGY: {total_energy}")

    ecutwfc_values = np.array(ecutwfc_values)
    energy_values = np.array(energy_values)
    
    # save the results as a .txt file for reference
    header="ecutwfc total_energy_eV"
    np.savetxt(DATA_DIR / "ecutwfc Convergence.txt", np.column_stack((ecutwfc_values, energy_values)), header=header)
    
    _graph_convergence(ecutwfc_values, energy_values, field, header, "ecutwfc")
    
    return (ecutwfc_values, energy_values)
