import qe.pw as pw
import qe.dos as dos
import numpy as np
from config import *

from pathlib import Path
import matplotlib.pyplot as plt
import shutil

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
        ax.set_xlabel(f"{parameter} Parameter x --> (x, x, 1)")
    
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
        step = 3
        label = "scf"

    elif computation == "nscf":
        calc_path = DATA_DIR / "nscf"
        l_bound = NSCF_LOWER_BOUND
        u_bound = NSCF_UPPER_BOUND
        step = 6
        label = "nscf"
        
        pw._calculate("scf", structure, calc_path, (NSCF_LOWER_BOUND, NSCF_LOWER_BOUND, 1), efield=field)
    
    for i in range(l_bound, u_bound+1, step):
        kgrid = (i, i, 1)
                
        if computation == "scf":
            scf = pw._calculate("scf", structure, calc_path, kgrid, efield=field)
            
            total_energy = scf.get_potential_energy()  # extract total energy
            energy_values.append(total_energy)
            
            print(f"KGRID = ({i}, {i}, 1) --> TOTAL ENERGY: {total_energy}")
            
        elif computation == "nscf":
            nscf = pw._calculate("nscf", structure, calc_path, kgrid, efield=field)
            
            fermi_energy = nscf.calc.get_fermi_level()
            dos_ = dos.calculate(calc_path, fermi_energy)
            doses.append((dos_.energy, dos_.dos))
            
            print(f"KGRID = ({i}, {i}, 1) --> DOS: {dos_}")
        
        kgrid_values.append(i)
    
    kgrid_values = np.array(kgrid_values)
    energy_values = np.array(energy_values)

    if computation == "scf":
        data = np.column_stack((kgrid_values, energy_values))
        np.savetxt(DATA_DIR / "scf kgrid Convergence.txt", data, header="kgrid total_energy_eV")
        _graph_convergence(kgrid_values, energy_values, field, "Total Energy (eV)", "kgrid_scf")

    # plot dos for each kgrid density compared against last density i.e. [-1]
    if computation == "nscf":
        dos_ref = doses[-1][1]  # take densest grid's dos as reference
        error = [np.abs(dos - dos_ref).sum() / dos_ref.sum() for _, dos in doses] # calc relative error for all grids
        data = np.column_stack((kgrid_values, error))
        np.savetxt(DATA_DIR / "nscf kgrid Convergence.txt", data, header="kgrid relative_dos_error")
        _graph_convergence(kgrid_values, error, field, "Relative DOS error", "kgrid_nscf_dos")
        
    # delete large metadata
    shutil.rmtree(calc_path, ignore_errors=True)  # don't crash if folder not found
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
        
        calc_path = DATA_DIR / "ecutwfc" / str(ecutwfc)
        
        scf = pw._calculate("scf", structure, calc_path, KGRID, ecutwfc=ecutwfc, ecutrho=4*ecutwfc, efield=field)
        
        total_energy = scf.get_potential_energy()  # extract total energy
        
        ecutwfc_values.append(ecutwfc)
        energy_values.append(total_energy)
        
        # delete large metadata
        shutil.rmtree(calc_path, ignore_errors=True)  # don't crash if folder not found
    
        print(f"ECUTWFC = {ecutwfc} --> TOTAL ENERGY: {total_energy}")

    ecutwfc_values = np.array(ecutwfc_values)
    energy_values = np.array(energy_values)
    
    # save the results as a .txt file for reference
    header="ecutwfc total_energy_eV"
    np.savetxt(DATA_DIR / "ecutwfc Convergence.txt", np.column_stack((ecutwfc_values, energy_values)), header=header)
    
    _graph_convergence(ecutwfc_values, energy_values, field, header, "ecutwfc")

    return (ecutwfc_values, energy_values)
