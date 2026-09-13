from graphene import GrapheneStructure
import qe.pw as pw
import qe.dos as dos
import convergence
import plotter

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
    
    energies = [0]
    
    for eamp in energies:
        
        path = BILAYER /  f"field_{str(eamp)}"
        
        PATH_COUPLED = path / "coupled"
        PATH_BOTTOM = path / "bottom"
        PATH_TOP = path / "top"
        
        
        print(f"\nE-field = {str(eamp)}au")
        print("-" * 30)

        print("CREATING GRAPHENE BILAYERS")
        bilayer = graphene.bilayer()
        
        print("RELAXING COUPLED BILAYER")
        relaxed_coupled = pw.relax(bilayer, PATH_COUPLED, eamp)
        
        print("\nEXTRACTING FROZEN LAYERS")
        bilayer_bottom, bilayer_top = graphene.isolate_bilayer(relaxed_coupled)
        
        results = {}
        
        print("\nCOUPLED LAYERS COMPUTATIONS")
        scf_coupled = pw.scf(relaxed_coupled, PATH_COUPLED, eamp)
        nscf_coupled = pw.nscf(relaxed_coupled, PATH_COUPLED, eamp)
        fermi_e_coupled = nscf_coupled.calc.get_fermi_level()
        dos_coupled = dos.calculate(PATH_COUPLED, fermi_e_coupled)
        
        results["coupled"] = [dos_coupled[0], dos_coupled[1], fermi_e_coupled]
        
        
        print("\nBOTTOM LAYER COMPUTATIONS")
        scf_bottom = pw.scf(bilayer_bottom, PATH_BOTTOM, eamp)
        nscf_bottom = pw.nscf(bilayer_bottom, PATH_BOTTOM, eamp)
        fermi_e_bottom = nscf_bottom.calc.get_fermi_level()
        dos_bottom = dos.calculate(PATH_BOTTOM, fermi_e_bottom)

        results["bottom"] = [dos_bottom[0], dos_bottom[1], fermi_e_bottom]
        
        
        print("\nTOP LAYER COMPUTATIONS")
        scf_top = pw.scf(bilayer_top, PATH_TOP, eamp)
        nscf_top = pw.nscf(bilayer_top, PATH_TOP, eamp)
        fermi_e_top = nscf_top.calc.get_fermi_level()
        dos_top = dos.calculate(PATH_TOP, fermi_e_top)

        results["top"] = [dos_top[0], dos_top[1], fermi_e_top]
        
        
        plotter.plot_dos(results, eamp)


if __name__ == "__main__":
    main()
