from ase.io import write, read
from pathlib import Path

import qe.runner as runner


ROOT = Path(__file__).resolve().parents[1]

# Carbon pseudopotential (from https://sssp.materialscloud.org/pseudopotentials/PBE/efficiency)
PSEUDO_DIR = ROOT / "pseudo"
PSEUDO = "C.pbe-n-kjpaw_psl.1.0.0.UPF"

# parameters
ECUTRHO = 400.0
CONV_THRESHOLD = 1.0e-8
DEGAUSS = 0.005
SMEARING = "gauss"
UNOCCUPIED_BANDS = 8  # number of unoccupied bands to run the calculations for


def input_data(calculation, data_path, ecutwfc, nbnd, prefix, efield):
    """
    https://www.quantum-espresso.org/Doc/INPUT_PW.html#id3
    Function called by write_input to create input file
    """

    control = {"calculation": calculation, 
               "prefix": prefix,
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
    if calculation in ["relax", "scf", "nscf"]:
        system["occupations"] = "smearing"  # smoothing out Fermi level 0-1 jump - prevent oscillation of SCF
        system["smearing"] = SMEARING  # 'mv'
        system["degauss"] = DEGAUSS
        
    if efield != 0:
        """
        Perpendicular saw-like potential 
        
        The saw-like potential increases with slope eamp in the region from (emaxpos+eopreg-1) to (emaxpos), 
        then decreases to 0 until (emaxpos+eopreg), in units of the crystal vector edir. 
        Important: the change of slope of this potential must be located in the empty region, 
        or else unphysical forces will result.
        """
        
        control["tefield"] = True  # saw-like potential simulating an E-field is added to the bare ionic potential
        system["edir"] = 3  # z-direction
        system["emaxpos"] = 0.9  # Position of the maximum of the saw-like potential along crystal axis edir
        system["eopreg"] = 0.1  # Zone in the unit cell where the saw-like potential decreases
        system["eamp"] = efield
        
    namelists = {"control": control, "system": system, "electrons": electrons}

    if calculation == "relax":
        namelists["ions"] = {"ion_dynamics": "bfgs"}        

    return namelists


def write_input(path, structure, calculation, data_path, kpts, ecutwfc, nbnd, prefix, efield):
    """
    Write QE input file using ASE
    """
    
    write(path, structure, format="espresso-in", input_data=input_data(calculation, data_path, ecutwfc, nbnd, prefix, efield), pseudopotentials={"C": PSEUDO}, kpts=kpts)


def read_output(path, index=-1):
    structure = read(path, format="espresso-out", index=index)  # index=-1 gets last structure
        
    return structure


def calculate(structure, calculation, path, kpts, ecutwfc, efield=0):
    
    path.mkdir(parents=True, exist_ok=True)
    
    data_path = path / "data"
    input_path = path / f"{calculation}.pwi"
    output_path = path / f"{calculation}.pwo"

    # 2 * (number of atoms in Atoms object) -> number of occupied bands
    nbnd = 2 * len(structure) + UNOCCUPIED_BANDS
    
    print(f"\nCREATING {input_path.name}")
    write_input(input_path, structure, calculation, data_path, kpts, ecutwfc, nbnd, path.name, efield)

    print(f"RUNNING pw.x WITH {input_path.name}")
    runner.run("pw.x", input_path, output_path)

    print(f"READING {output_path.name}")
    output = read_output(output_path)
    
    return output
