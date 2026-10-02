import qe.pw as pw
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent

FIG_DIR = ROOT / "convergence" / "figures"
DATA_DIR = ROOT / "convergence" / "data"

FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

ECUTWFC = 100.0
KGRID = (12, 12, 1)


def _graph_convergence(values, energies, field, parameter):
    
    fig, ax = plt.subplots()
    
    ax.plot(values, energies, "o-")
    ax.set_title(f"{parameter} Parameter Convergence for {parameter} at E-field={field}au")
    
    if parameter == "kgrid":
        ax.set_xlabel(f"{parameter} Parameter x -> (x, x, 1)")
    
    else:
        ax.set_xlabel(f"{parameter} Parameter")

    ax.set_ylabel("Total Energy (eV)")
    
    fig.savefig(FIG_DIR / f"{parameter} Convergence, E-field = {field}au", bbox_inches="tight")
    plt.close(fig)


def test_kgrid(structure, field, upper):
    """
    Convergence testing of kgrid by iterating through from values of x and y 
    ranging from 2 to 16, with z being kept constant at 1.
    Involves running the scf process and plotting the total energy against kgrid values.
    return: scf.pwo file read using ASE
    """
    
    # list to store values for plotting
    kgrid_values = []
    energy_values = []
    
    for i in range(2, upper+1, 2):
        
        kgrid = (i, i, 1)
        
        calculation_path = DATA_DIR / "kgrid" / str(i)
        
        scf = pw._calculate("scf", structure, calculation_path, kgrid, efield=field)
                
        total_energy = scf.get_potential_energy()  # extract total energy
        
        kgrid_values.append(i)
        energy_values.append(total_energy)
        print(f"KGRID = ({i}, {i}, 1) --> TOTAL ENERGY: {total_energy}")
    
    kgrid_values = np.array(kgrid_values)
    energy_values = np.array(energy_values)
    
    # save the results as a .txt file for reference
    np.savetxt(DATA_DIR / "kgrid Convergence.txt", np.column_stack((kgrid_values, energy_values)), header="kgrid total_energy_eV")
    
    _graph_convergence(kgrid_values, energy_values, field, "kgrid")
    
    return (kgrid_values, energy_values)


def test_ecutwfc(structure, field, upper):
    """
    Convergence testing of ecutwfc by iterating through values 10 to 100 in increments of 10.
    Involves running the scf process and plotting the total energy against ecutwfc values.
    return: scf.pwo file read using ASE
    """
    
    # list to store values for plotting
    ecutwfc_values = []
    energy_values = []
    
    for ecutwfc in range(10, upper+1, 10):
        
        calculation_path = DATA_DIR / "ecutwfc" / str(ecutwfc)
        
        scf = pw._calculate("scf", structure, calculation_path, kgrid, ecutwfc)
        
        total_energy = scf.get_potential_energy()  # extract total energy
    
        ecutwfc_values.append(ecutwfc)
        energy_values.append(total_energy) 
        print(f"ECUTWFC = {ecutwfc} --> TOTAL ENERGY: {total_energy}")

    ecutwfc_values = np.array(ecutwfc_values)
    energy_values = np.array(energy_values)
    
    # save the results as a .txt file for reference
    np.savetxt(DATA_DIR / "ecutwfc Convergence.txt", np.column_stack((ecutwfc_values, energy_values)), header="kgrid total_energy_eV")
    
    _graph_convergence(ecutwfc_values, energy_values, field, "ecutwfc")
    
    return (ecutwfc_values, energy_values)
