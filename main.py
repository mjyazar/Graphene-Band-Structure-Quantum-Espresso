from graphene import GrapheneStructure
import qe
import convergence

from ase.spectrum.band_structure import get_band_structure, BandStructure
import matplotlib.pyplot as plt
from pathlib import Path


ROOT = Path(__file__).resolve().parent

# Paths
OUT_DIR = ROOT / "outputs"
FIG_DIR = OUT_DIR / "figures"
MONOLAYER = OUT_DIR / "monolayer"
BILAYER = OUT_DIR / "bilayer"

# Ensure folders exist
OUT_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)
MONOLAYER.mkdir(exist_ok=True)
BILAYER.mkdir(exist_ok=True)


# RUN_QE may be False if script already ran and want to work with existing files
# True if running for the first time or want to create new files with new parameters
RUN_QE = True
RUN_CONVERGENCE = True

# parameters
BANDPATH = 'GMKG'
ecutwfc = 100.0
kgrid = (12, 12, 1)


def band_path(structure, path=BANDPATH):
    """
    Determining the band- (k-) path
    
    pbc: Whether cell is periodic in each direction
         If cell has three nonzero cell vectors, use pbc=[1, 1, 0] to request a 2D bandpath
    """
    
    return structure.cell.bandpath(path=path, pbc=[True, True, False], npoints=100)


def print_structure_data(name, structure, relaxed, band_structure, total_energy, fermi_energy):
    print(f"\n{name.upper()}:")
    print(structure)

    print("\nCell:")
    print(structure.cell)

    print("\nPositions:")
    print(structure.positions)

    print("\nDistances:")
    print(structure.get_all_distances(mic=True))

    print(f"\nNumber of Atoms: {len(structure)}")

    print("\nScaled positions:")
    print(structure.get_scaled_positions())

    # use mic=True to use the Minimum Image Convention
    # vector=True gives the distance vector (from a0 to a1)
    print(f"\nC-C Distance (A): {structure.get_distance(0, 1, mic=True)}")
    print(f"\nC-C Distance (A): {structure.get_distance(0, 1, vector=True)}")
    
    print("\nFinal Relaxed Atomic Coordinates:")
    print(relaxed.positions)
    
    print("\nBand energy array shape:")
    print(band_structure.energies.shape)

    band_energies = band_structure.energies
    print("\nBand energies:")
    print(band_energies)

    print(f"\nTotal Energy: {total_energy} eV")
    print(f"\nFermi Energy: {fermi_energy} eV")


def plot_band_structure(name, bandpath, energies):
    """
    Plot and save band structure
    """
    
    band_structure = BandStructure(path=bandpath, energies=energies, reference=0.0)  # reference is now zero after shifting
    
    ax = band_structure.plot()
    ax.set_ylim(-10, 10)
    ax.set_ylabel("Energy - $E_F$ (eV)")
    ax.set_title(f"Graphene {name.capitalize()} Band Structure")
    
    plt.savefig(FIG_DIR / f"{name.capitalize()} Graphene Band Structure.png", dpi=300)
    plt.close()
    

def main():
    graphene = GrapheneStructure()

    print("\nCREATING GRAPHENE MONOLAYER")
    monolayer = graphene.monolayer()
    print("CREATED GRAPHENE MONOLAYER")

    print("\nCREATING GRAPHENE BILAYER")
    bilayer = graphene.bilayer()
    print("CREATED GRAPHENE BILAYER\n")

    atoms = {"monolayer": (monolayer, MONOLAYER), "bilayer": (bilayer, BILAYER)}

    for name, (structure, path) in atoms.items():
        print(f"\nGRAPHENE {name.upper()}")

        if RUN_QE:
                                
            relaxed = qe.calculate(structure, "relax", path, kgrid, ecutwfc)
            scf = qe.calculate(relaxed, "scf", path, kgrid, ecutwfc)
            bandpath = band_path(relaxed)
            bands = qe.calculate(relaxed, "bands", path, bandpath, ecutwfc)
            
        else:
            
            relaxed = qe.read_output(path / "relax.pwo")
            scf = qe.read_output(path / "scf.pwo")
            bandpath = band_path(relaxed)
            bands = qe.read_output(path / "bands.pwo")

        print("GETTING BAND STRUCTURE")
        band_structure = get_band_structure(atoms=bands, calc=bands.calc)
        
        print("CALCULATING ENERGIES")
        total_energy = scf.get_potential_energy()
        fermi_energy = scf.calc.get_fermi_level()
        energies = band_structure.energies - fermi_energy

        print("PLOTTING BAND STRUCTURE")
        plot_band_structure(name, bandpath, energies)

        print_structure_data(name, structure, relaxed, band_structure, total_energy, fermi_energy)

    
    if RUN_CONVERGENCE:
        
        for name, (structure, path) in atoms.items():
            print(f"\nCONVERGENCE TESTING {name.upper()} GRAPHENE")
            convergence_path = path / "convergence"
            
            print("\nRUNNING KGRID CONVERGENCE TEST")
            kgrid_values, kgrid_energies = convergence.test_kgrid(structure, name, convergence_path, 16, ecutwfc)
            print("\nRUNNING ECUTWFC CONVERGENCE TEST")
            ecutwfc_values, ecutwfc_energies = convergence.test_ecutwfc(structure, name, convergence_path, 100, kgrid)


if __name__ == "__main__":
    main()
