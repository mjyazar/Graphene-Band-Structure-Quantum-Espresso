import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ase.spectrum.band_structure import BandStructure

ROOT = Path(__file__).resolve().parent

OUT_DIR = ROOT / "outputs"
FIG_DIR = OUT_DIR / "figures"

FIG_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)


def plot_band_structure(bandpath, energies, name):
    """
    Plot and save band structure
    """
    
    band_structure = BandStructure(path=bandpath, energies=energies, reference=0.0)  # reference is now zero after shifting
    
    fig, ax = plt.subplots()

    ax = band_structure.plot()
    assert ax is not None  # ensures BandStructure doesn't return None to prevent an error
    ax.set_ylim(-10, 10)
    ax.set_ylabel("Energy - $E_F$ (eV)")
    ax.set_title(f"Graphene {name.capitalize()} Band Structure")
    
    plt.savefig(FIG_DIR / f"{name.capitalize()} Graphene Band Structure.png", dpi=300)
    plt.close(fig)


def plot_dos_comparison(results, field, window=(-5, 5)):
    print("\nPLOTTING DOS")
    
    fig, ax = plt.subplots()

    for name, (energy, dos, fermi_energy) in results.items():
        ax.plot(energy - fermi_energy, dos, label=name)    
        ax.fill_between(energy - fermi_energy, 0, dos, where=(energy - fermi_energy < 0), facecolor='red', alpha=0.15)

    
    ax.axvline(0, linewidth=0.8, linestyle="--")
    ax.set_title(f"Graphene 1D DOS Comparison")
    ax.set_xlabel(r"$E - E_{Fermi}$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.set_xlim(window)
    # ax.set_ylim()
    
    ax.text(0, 2.5, "Fermi energy", fontsize=16, rotation=90)

    fig.tight_layout()
    plt.savefig(FIG_DIR / f"DOS Comparison {field}au", dpi=300, )
    plt.close(fig)


def plot_dos_added(results, field, delta_e=0.01):
    
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

    ax.set_title(f"DOS Comparison")
    ax.set_xlabel(r"$Energy$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.legend()
    
    fig.tight_layout()
    
    plt.savefig(FIG_DIR / f"DOS Added {field}au", dpi=300, )
    plt.close(fig)


def plot_2d_dos():
    pass


def plot_dos(results, field):
    
    plot_dos_comparison(results, field)
    plot_dos_added(results, field)
    plot_2d_dos()
