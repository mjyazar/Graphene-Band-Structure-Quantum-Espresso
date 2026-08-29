from ase.io import write, read
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent

# Carbon pseudopotential (from https://sssp.materialscloud.org/pseudopotentials/PBE/efficiency)
PSEUDO_DIR = ROOT / "pseudo"
PSEUDO = "C.pbe-n-kjpaw_psl.1.0.0.UPF"

# parameters
ECUTRHO = 400.0
CONV_THRESHOLD = 1.0e-8
DEGAUSS = 0.01
SMEARING = "mv"
UNOCCUPIED_BANDS = 4  # number of unoccupied bands to run the calculations for


def input_data(calculation, data_path, ecutwfc, nbnd):
    """
    Functioncalled by write_input to create input file
    """

    control = {"calculation": calculation, 
               "prefix": "graphene",
               "verbosity": "high",  # amount of information written in QE output - high -> slower
               "outdir": str(data_path),  # path for temporary/intermediate calculation files
               "pseudo_dir": str(PSEUDO_DIR),  # directory containing pseudopotentials.
               "tprnfor": True}  # print atomic forces on each atom
    
    system = {"assume_isolated": "2D",
              "ecutwfc": ecutwfc,  # kinetic energy (1Ry ~13.6eV) upto which plane waves are included in compuation
              "ecutrho": ECUTRHO,  # kinetic energy upto which electron charge desnity is computed
              "nbnd": nbnd}  # number of energy/eigenbands calculated at every k point - C: 4 valence e -> 2 atoms
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


def write_input(path, structure, calculation, data_path, kpts, ecutwfc, nbnd):
    """
    Write QE input file using ASE
    """
    
    write(path, structure, format="espresso-in", input_data=input_data(calculation, data_path, ecutwfc, nbnd), pseudopotentials={"C": PSEUDO}, kpts=kpts)


def run_qe(input_path, output_path):
    """
    Atomate the process of manually running pw.x << input_file.pwi >> output_file.pwo 
    in the terminal after creating an input file.
    """
    
    with open(input_path, mode="r") as input_file, open(output_path, mode="w") as output_file:
        subprocess.run(["pw.x"], stdin=input_file, stdout=output_file, check=True)


def read_output(path, index=-1):
    structure = read(path, format="espresso-out", index=index)  # index=-1 gets last structure
        
    return structure


def calculate(structure, calculation, path, kpts, ecutwfc):
    
    path.mkdir(parents=True, exist_ok=True)
    
    data_path = path / "data"
    input_path = path / f"{calculation}.pwi"
    output_path = path / f"{calculation}.pwo"

    # 2 * (number of atoms in Atoms object) -> number of occupied bands
    nbnd = 2 * len(structure) + UNOCCUPIED_BANDS
    
    write_input(input_path, structure, calculation, data_path, kpts, ecutwfc, nbnd)
    print(f"\nCREATED {input_path.name}")

    run_qe(input_path, output_path)
    print(f"RAN QE WITH {input_path.name}")
    
    output = read_output(output_path)
    print(f"READ {output_path.name}")
    
    return output
