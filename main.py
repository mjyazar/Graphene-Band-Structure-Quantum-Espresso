import numpy as np
from ase import Atoms
from ase.io import write, read
from ase.spectrum.band_structure import get_band_structure, BandStructure
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent

LATTICE_CONSTANT = 2.46
VACUUM = 10.0  # QE requires 2D Coulomb truncation of the cell to have min z-length ~10.58 A
BANDPATH = 'GMKG'

# Paths
PSEUDO_DIR = ROOT / "pseudo"
OUT_DIR = ROOT / "outputs"
TEMP_DIR = OUT_DIR / "temp"
FIG_DIR = OUT_DIR / "figs"
CONV_DIR = OUT_DIR / "conv"

# Ensure folders exist
OUT_DIR.mkdir(exist_ok=True)
TEMP_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)
CONV_DIR.mkdir(exist_ok=True)

# Carbon pseudopotential (from https://sssp.materialscloud.org/pseudopotentials/PBE/efficiency)
PSEUDO = "C.pbe-n-kjpaw_psl.1.0.0.UPF"

# parameters
ecutwfc = 100.0
ECUTRHO = 400.0
NBND = 10  # number of bands
CONV_THRESHOLD = 1.0e-8
DEGAUSS = 0.01
SMEARING = "mv"
kgrid = (12, 12, 1)


# May make False if script already ran and want to work with existing files
# True if running for the first time or want to create new files with new parameters
RUN_QE = True



def graphene(a=LATTICE_CONSTANT, vacuum=VACUUM):
    """
    Creating graphene as an Atoms object
    Monolayer graphene: honeycomb, 2 atoms/unit cell
    """
    
    cell = np.array([[a, 0, 0], 
                    [-a/2, a*np.sqrt(3)/2, 0], 
                    [0, 0, 2*vacuum]])
    
    scaled_positions = np.array([[0, 0, 0.5], [2/3, 1/3, 0.5]])
    
    return Atoms(symbols="C2", cell=cell, scaled_positions=scaled_positions, pbc=True)
    

def band_path(graphene, path=BANDPATH):
    """
    Determining the band- (k-) path
    
    pbc: Whether cell is periodic in each direction
         If cell has three nonzero cell vectors, use pbc=[1, 1, 0] to request a 2D bandpath
    """
    
    return graphene.cell.bandpath(path=path, pbc=[True, True, False], npoints=100)


def input_data(calculation, ecutwfc=ecutwfc, ecutrho=ECUTRHO):
    """
    Functioncalled by write_input to create input file
    """

    control = {"calculation": calculation, 
               "prefix": "graphene",
               "verbosity": "high",  # amount of information written in QE output - high -> slower
               "outdir": str(TEMP_DIR),  # path for temporary/intermediate calculation files
               "pseudo_dir": str(PSEUDO_DIR),  # directory containing pseudopotentials.
               "tprnfor": True}  # print atomic forces on each atom
    
    system = {"assume_isolated": "2D",
              "ecutwfc": ecutwfc,  # kinetic energy (1Ry ~13.6eV) upto which plane waves are included in compuation
              "ecutrho": ecutrho,  # kinetic energy upto which electron charge desnity is computed
              "nbnd": NBND}  # number of energy/eigenbands calculated at every k point - C: 4 valence e -> 2 atoms
                             # = 8 electrons = 4 filled bands (spin degeneracy). Anything above is empty states
    
    electrons = {"conv_thr": CONV_THRESHOLD,  # scf convergence thmaximum number of electronic SCF iterations allowed treshold (Ry)
                 "mixing_beta": 0.7,  # How much electron density updates with scf
                 "electron_maxstep": 200}  # maximum number of electronic SCF iterations allowed
    
    # additional parameters if doing relaxation or s
    if calculation in ["relax", "scf"]:
        system["occupations"] = "smearing"  # smoothing out Fermi level 0-1 jump - prevent oscillation of SCF
        system["smearing"] = SMEARING  # 'mv'
        system["degauss"] = DEGAUSS
        
    namelists = {"control": control, "system": system, "electrons": electrons}

    if calculation == "relax":
        namelists["ions"] = {"ion_dynamics": "bfgs"}

    return namelists


def write_input(path, object, calculation, kpts, ecutwfc, ecutrho=ECUTRHO):
    """
    Write QE input file using ASE
    """
    
    write(path, object, format="espresso-in", input_data=input_data(calculation, ecutwfc, ecutrho), pseudopotentials={"C": PSEUDO}, kpts=kpts)


def run_qe(input, output):
    """
    Atomate the process of manually running pw.x << input_file.pwi >> output_file.pwo 
    in the terminal after creating an input file.
    """
    
    with open(input, mode="r") as input_file, open(output, mode="w") as output_file:
        subprocess.run(["pw.x"], stdin=input_file, stdout=output_file, check=True)


def convergence_testing_kgrids(object, upper):
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

        write_input(OUT_DIR / "scf.pwi", object, "scf", kgrid, ecutwfc)
        print("\nCREATED scf.pwi")
        run_qe(OUT_DIR / "scf.pwi", OUT_DIR / "scf.pwo")   
        print("RAN QE WITH scf.pwi")
        
        scf = read(OUT_DIR / "scf.pwo", format="espresso-out")
        print("READ scf.pwo")
        
        total_energy = scf.get_potential_energy()  # extract total energy
        
        kgrid_values.append(i)
        energy_values.append(total_energy)
        print(f"KGRID = ({i}, {i}, 1) --> TOTAL ENERGY: {total_energy}")
        
        
    kgrid_values = np.array(kgrid_values)
    energy_values = np.array(energy_values)
    
    # save the results as a .txt file for reference
    np.savetxt(CONV_DIR / "KGRID Convergence.txt", np.column_stack((kgrid_values, energy_values)), header="kgrid total_energy_eV")
    
    fig, ax = plt.subplots()
    
    ax.plot(kgrid_values, energy_values, "o-")
    ax.set_title("kgrid Parameter Convergence")
    ax.set_xlabel("kgrid Parameter x -> (x, x, 1)")
    ax.set_ylabel("Total Energy (eV)")
    
    fig.savefig(FIG_DIR / "KGRIDS Convergence", bbox_inches="tight")
    plt.close(fig)
    
    return scf


def convergence_testing_ecutwfc(object, upper):
    """
    Convergence testing of ecutwfc by iterating through values 10 to 100 in increments of 10.
    Involves running the scf process and plotting the total energy against ecutwfc values.
    return: scf.pwo file read using ASE
    """
    
    # list to store values for plotting
    ecutwfc_values = []
    energy_values = []
    
    for i in range(10, 100, 10):
        ecutwfc = i
        write_input(OUT_DIR / "scf.pwi", object, "scf", kgrid, ecutwfc)
        print("\nCREATED scf.pwi")
        run_qe(OUT_DIR / "scf.pwi", OUT_DIR / "scf.pwo")   
        print("RAN QE WITH scf.pwi")
        
        scf = read(OUT_DIR / "scf.pwo", format="espresso-out")
        print("READ scf.pwo")

        total_energy = scf.get_potential_energy()  # extract total energy
    
        ecutwfc_values.append(i)
        energy_values.append(total_energy) 
        print(f"ECUTWFC = {i} --> TOTAL ENERGY: {total_energy}")

    energy_values = np.array(energy_values)
    
    # save the results as a .txt file for reference
    np.savetxt(CONV_DIR / "ECUTWFC Convergence.txt", np.column_stack((ecutwfc_values, energy_values)), header="kgrid total_energy_eV")
    
    fig, ax = plt.subplots()
    
    ax.plot(ecutwfc_values, energy_values, "o-")
    ax.set_title("ECUTWFC Parameter Convergence")
    ax.set_xlabel("ECUTWFC Parameter")
    ax.set_ylabel("Total Energy (eV)")
    
    fig.savefig(FIG_DIR / "ECUTWFC Convergence", bbox_inches="tight")
    plt.close(fig)
    
    return scf

    

if __name__ == "__main__":
    print("\nCREATING GRAPHENE...")
    graphene = graphene()
    print("GRAPHENE CREATED...")
    
    # relaxation process
    if RUN_QE:
        write_input(OUT_DIR / "relax.pwi", graphene, "relax", kgrid, ecutwfc)
        print("\nCREATED relax.pwi")
        run_qe(OUT_DIR / "relax.pwi", OUT_DIR / "relax.pwo")
        print("RAN QE WITH relax.pwi")
    
    relax = read(OUT_DIR / "relax.pwo", format="espresso-out", index=-1)  # index=-1 gets last structure
    print("READ relax.pwo")
    
    # Before implementing convergence testing:
    # if RUN_QE:
    #     write_input(OUT_DIR / "scf.pwi", relax, "scf", kgrid)
    #     run_qe(OUT_DIR / "scf.pwi", OUT_DIR / "scf.pwo")
        
    # scf = read(OUT_DIR / "scf.pwo", format="espresso-out")
    
    # Self Consistent Field analysis
    if RUN_QE:
        #scf = convergence_testing_kgrids(relax, 16)
        #print("\nRAN KGRIDS CONVERGENCE")
        scf = convergence_testing_ecutwfc(relax, 100)
        print("\nRAN ECUTWFC CONVERGENCE")
    
    bandpath = band_path(relax)
    
    # creating graphene bands
    if RUN_QE:
        write_input(OUT_DIR / "bands.pwi", relax, "bands", bandpath, ecutwfc)
        print("\nCREATED bands.pwi")
        run_qe(OUT_DIR / "bands.pwi", OUT_DIR / "bands.pwo")
        print("RAN QE WITH bands.pwi")

    bands = read(OUT_DIR / "bands.pwo", format="espresso-out")
    print("READ bands.pwo")

    band_structure = get_band_structure(atoms=bands, calc=bands.calc)
    print("\nBAND STRUCTURE COMPUTED")
    
    # print some data:
    print("\nAtoms:")
    print(graphene)

    print("\nCell:")
    print(graphene.cell)

    print("\nPositions:")
    print(graphene.positions)

    print("\nDistances:")
    print(graphene.get_all_distances(mic=True))

    print(f"\nNumber of Atoms: {len(graphene)}")

    print("\nScaled positions:")
    print(graphene.get_scaled_positions())

    # use mic=True to use the Minimum Image Convention
    # vector=True gives the distance vector (from a0 to a1)
    print(f"\nC-C Distance (A): {graphene.get_distance(0, 1, mic=True)}")
    print(f"\nC-C Distance (A): {graphene.get_distance(0, 1, vector=True)}")
    
    print("\nFinal Relaxed Atomic Coordinates:")
    print(relax.positions)
    
    print("\nBand energy array shape:")
    print(band_structure.energies.shape)

    band_energies = band_structure.energies
    print("\nBand energies:")
    print(band_energies)

    total_energy = scf.get_potential_energy()
    fermi_energy = scf.calc.get_fermi_level()
    print(f"\nTotal Energy: {total_energy} eV")
    print(f"\nFermi Energy: {fermi_energy} eV")

    energies = band_structure.energies - fermi_energy
    
    # plot and save band structure
    band_structure = BandStructure(path=bandpath, energies=energies, reference=0.0)  # reference is now zero after shifting
    ax = band_structure.plot()
    ax.set_ylim(-10, 10)
    ax.set_ylabel("Energy - $E_F$ (eV)")
    ax.set_title("Graphene Band Structure")
    plt.savefig(FIG_DIR / "Graphene Band Structure.png", dpi=300)
