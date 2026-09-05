from graphene import GrapheneStructure
import qe.pw as pw
import qe.dos as dos
import convergence
import plotting

from ase.spectrum.band_structure import get_band_structure, BandStructure
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Paths
OUT_DIR = ROOT / "outputs"
#MONOLAYER = OUT_DIR / "monolayer"
BILAYER = OUT_DIR / "bilayer"

# Ensure folders exist
OUT_DIR.mkdir(exist_ok=True)
#MONOLAYER.mkdir(exist_ok=True)
BILAYER.mkdir(exist_ok=True)

# RUN_QE may be False if script already ran and want to work with existing files
# True if running for the first time or want to create new files with new parameters
RUN_QE = True
RUN_CONVERGENCE = True

# parameters
BANDPATH = 'GMKG'
ecutwfc = 100.0
KGRID = (12, 12, 1)
KGRID_DENSE = (15, 15, 1)


def band_path(structure, path=BANDPATH):
    """
    Determining the band- (k-) path
    
    pbc: Whether cell is periodic in each direction
         If cell has three nonzero cell vectors, use pbc=[1, 1, 0] to request a 2D bandpath
    """
    
    return structure.cell.bandpath(path=path, pbc=[True, True, False], npoints=100)


def print_structure_data(name, structure, relaxed, band_structure, total_eamp, fermi_eamp):
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
    
    print("\nBand eamp array shape:")
    print(band_structure.energies.shape)

    # band_energies = band_structure.energies
    # print("\nBand energies:")
    # print(band_energies)

    print(f"\nTotal eamp: {total_eamp} eV")
    print(f"\nFermi eamp: {fermi_eamp} eV")


def main():
    graphene = GrapheneStructure()
    
    energies = [0, 0.005, 0.01]
    
    for eamp in energies:
        
        path = BILAYER /  f"field_{str(eamp)}Ry"
        
        PATH_COUPLED = path / "coupled"
        PATH_BOTTOM = path / "bottom"
        PATH_TOP = path / "top"
        
        
        print(f"\nE-field = {str(eamp)}Ry")
        print("-" * 30)

        print("CREATING GRAPHENE BILAYERS")
        bilayer = graphene.bilayer()
        
        print("RELAXING COUPLED BILAYER")
        relaxed_coupled = pw.calculate(bilayer, "relax", PATH_COUPLED, KGRID, ecutwfc, eamp)
        
        print("EXTRACTING FROZEN LAYERS")
        bilayer_bottom, bilayer_top = graphene.isolate_bilayer(relaxed_coupled)
        
        results = {}
        
        print("\nCOUPLED LAYERS COMPUTATIONS")
        scf_coupled = pw.calculate(relaxed_coupled, "scf", PATH_COUPLED, KGRID, ecutwfc, eamp)
        nscf_coupled = pw.calculate(relaxed_coupled, "nscf", PATH_COUPLED, KGRID_DENSE, ecutwfc, eamp)
        dos_coupled = dos.calculate(PATH_COUPLED)
        fermi_e_coupled = nscf_coupled.calc.get_fermi_level()
        
        results["coupled"] = [dos_coupled[0], dos_coupled[1], fermi_e_coupled]
        
        print("\nBOTTOM LAYER COMPUTATIONS")
        scf_bottom = pw.calculate(bilayer_bottom, "scf", PATH_BOTTOM, KGRID, ecutwfc, eamp)
        nscf_bottom = pw.calculate(bilayer_bottom, "nscf", PATH_BOTTOM, KGRID_DENSE, ecutwfc, eamp)
        dos_bottom = dos.calculate(PATH_BOTTOM)
        fermi_e_bottom = nscf_bottom.calc.get_fermi_level()

        results["bottom"] = [dos_bottom[0], dos_bottom[1], fermi_e_bottom]
        
        print("\nTOP LAYER COMPUTATIONS")
        scf_top = pw.calculate(bilayer_top, "scf", PATH_TOP, KGRID, ecutwfc, eamp)
        nscf_top = pw.calculate(bilayer_top, "nscf", PATH_TOP, KGRID_DENSE, ecutwfc, eamp)
        dos_top = dos.calculate(PATH_TOP)
        fermi_e_top = nscf_top.calc.get_fermi_level()

        results["top"] = [dos_top[0], dos_top[1], fermi_e_top]
        
        plotting.plot_dos(results, eamp)


if __name__ == "__main__":
    main()
