import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ase.spectrum.band_structure import BandStructure

ROOT = Path(__file__).resolve().parent

OUT_DIR = ROOT / "outputs"
FIG_DIR = OUT_DIR / "figures"

FIG_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)

WINDOW = (-10, 10)


def plot_band_structure(bandpath, energies, name):
    """
    Plot and save band structure
    """
    print("\nPLOTTING BAND STRUCTURE")
    
    band_structure = BandStructure(path=bandpath, energies=energies, reference=0.0)  # reference is now zero after shifting
    
    fig, ax = plt.subplots()

    ax = band_structure.plot()
    assert ax is not None  # ensures BandStructure doesn't return None to prevent an error
    ax.set_ylim(-10, 10)
    ax.set_ylabel("Energy - $E_F$ (eV)")
    ax.set_title(f"Graphene {name.capitalize()} Band Structure")
    
    plt.savefig(FIG_DIR / f"{name.capitalize()} Graphene Band Structure.png", dpi=300)
    plt.close(fig)


def plot_individual_dos(results, field, window=WINDOW):

    for name, (energy, dos, fermi_energy) in results.items():
        print(f"\nPLOTTING {name} DOS")

        fig, ax = plt.subplots()

        ax.plot(energy - fermi_energy, dos, linewidth=0.75, color='red')
        ax.axvline(0, linestyle="--")
        ax.set_title(f"1D {name} Layer DOS {field}au")
        ax.set_xlabel(r"$E - E_{Fermi}$ (eV)")
        ax.set_ylabel("DOS (states/eV/cell)")
        ax.set_xlim(window)
        # ax.set_ylim()
        
        plt.savefig(FIG_DIR / f"{name.capitalize()} Layer DOS {field}au.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_dos_comparison(results, field, window=WINDOW):
    print("\nPLOTTING DOS COMPARISON")
    
    fig, ax = plt.subplots()

    for name, (energy, dos, fermi_energy) in results.items():
        ax.plot(energy - fermi_energy, dos, label=name)    
    
    ax.axvline(0, linewidth=0.8, linestyle="--")
    ax.set_title(f"Graphene 1D DOS Comparison {field}au")
    ax.set_xlabel(r"$E - E_{Fermi}$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.set_xlim(window)
    ax.legend()
 
    plt.savefig(FIG_DIR / f"DOS Comparison {field}au.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_dos_added(results, field, delta_e=0.01):
    print("\nPLOTTING DOS ADDED")

    energy_coupled, dos_coupled, fermi_coupled = results["coupled"]
    energy_bottom, dos_bottom, fermi_bottom = results["bottom"]
    energy_top, dos_top, fermi_top = results["top"]
    
    min_energy = max(energy_coupled.min(), energy_bottom.min(), energy_top.min())
    max_energy = min(energy_coupled.max(), energy_top.max(), energy_bottom.max())
    
    grid = np.arange(min_energy, max_energy, delta_e)
    
    dos_coupled_grid = np.interp(grid, energy_coupled, dos_coupled)
    dos_bottom_grid = np.interp(grid, energy_bottom, dos_bottom)
    dos_top_grid = np.interp(grid, energy_top, dos_top)

    dos_sum = dos_bottom_grid + dos_top_grid
    
    fig, ax = plt.subplots()

    ax.plot(grid, dos_coupled_grid, label="$DOS_{coupled}$")
    ax.plot(grid, dos_sum, label="$DOS_{top} + DOS_{bottom}$")

    ax.set_title(f"DOS Comparison {field}au")
    ax.set_xlabel(r"$Energy$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.legend()
        
    plt.savefig(FIG_DIR / f"DOS Added {field}au.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_2d_dos():
    pass


def plot_dos(results, field):
    
    plot_individual_dos(results, field)
    plot_dos_comparison(results, field)
    plot_dos_added(results, field)
    plot_2d_dos()
