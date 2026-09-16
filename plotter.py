import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ase.spectrum.band_structure import get_band_structure, BandStructure

from results import * 
from config import *

ROOT = Path(__file__).resolve().parent

OUT_DIR = ROOT / "outputs"
FIG_DIR = OUT_DIR / "figures"
POTENTIAL_DIR = FIG_DIR / "potential"
CHARGE_DENSITY_DIR = FIG_DIR / "charge_density"
DOS_DENSITY_DIR = FIG_DIR / "dos"

FIG_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)
POTENTIAL_DIR.mkdir(exist_ok=True)
CHARGE_DENSITY_DIR.mkdir(exist_ok=True)
DOS_DENSITY_DIR.mkdir(exist_ok=True)


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


def individual_dos(results, field, window=WINDOW):

    for name, (energy, dos, fermi_energy) in results["dos"].items():
        print(f"\nPLOTTING {name} DOS")

        fig, ax = plt.subplots()

        ax.plot(energy - fermi_energy, dos, linewidth=0.75, color='red')
        ax.axvline(0, linestyle="--")
        ax.set_title(f"1D {name} Layer DOS {field}au")
        ax.set_xlabel(r"$E - E_{Fermi}$ (eV)")
        ax.set_ylabel("DOS (states/eV/cell)")
        ax.set_xlim(window)
        # ax.set_ylim()
        
        plt.savefig(FIG_DIR / f"DOS_{name.capitalize()}_Layer_{field}au.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def fermi_aligned(results:Results, field, delta_e=0.01):
    print(f"PLOTTING FERMI-ALIGNED SUBTRACTED DOS")

    coupled = results.coupled
    bottom = results.bottom
    top = results.top
    
    energy_coupled = coupled.dos.energy - coupled.fermi_energy
    energy_bottom = bottom.dos.energy - bottom.fermi_energy
    energy_top = top.dos.energy - top.fermi_energy
    
    dos_coupled = coupled.dos.dos
    dos_bottom = bottom.dos.dos
    dos_top = top.dos.dos
    
    min_energy = max(energy_coupled.min(), energy_bottom.min(), energy_top.min())
    max_energy = min(energy_coupled.max(), energy_top.max(), energy_bottom.max())
    
    grid = np.arange(min_energy, max_energy, delta_e)
    
    dos_coupled_grid = np.interp(grid, energy_coupled, dos_coupled)
    dos_bottom_grid = np.interp(grid, energy_bottom, dos_bottom)
    dos_top_grid = np.interp(grid, energy_top, dos_top)
    
    dos_subtracted = dos_coupled_grid - dos_bottom_grid - dos_top_grid
    
    fig, ax = plt.subplots()

    # ax.plot(grid, dos_coupled_grid, label="$DOS_{coupled}$")
    ax.plot(grid, dos_subtracted, label="$DOS_{coupled} - DOS_{top} - DOS_{bottom}$")

    ax.set_title(f"Fermi-Aligned DOS Subtracted, E-field={field}au")
    ax.set_xlabel(r"$Energy$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.legend()
        
    plt.savefig(FIG_DIR / f"DOS_Subtracted_Fermi-aligned_{field}au.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def vacuum_aligned(results, field):
    pass



def dos_comparison(results, field, window=WINDOW):
    print("\nPLOTTING DOS COMPARISON")
    
    fig, ax = plt.subplots()

    for name, (energy, dos, fermi_energy) in results["dos"].items():
        ax.plot(energy - fermi_energy, dos, label=name)    
    
    ax.axvline(0, linewidth=0.8, linestyle="--")
    ax.set_title(f"Graphene 1D DOS Comparison, E-field={field}au")
    ax.set_xlabel(r"$E - E_{Fermi}$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.set_xlim(window)
    ax.legend()
 
    plt.savefig(FIG_DIR / f"DOS_Comparison_{field}au.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def potential_subtracted(results:Results, field):
    print("\nPLOTTING SUBTRACTED POTENTIAL")
    
    coupled = results.coupled.potential_averaged
    bottom = results.bottom.potential_averaged
    top = results.top.potential_averaged

    np.testing.assert_allclose(coupled.coordinates, bottom.coordinates)
    np.testing.assert_allclose(coupled.coordinates, top.coordinates)
    
    z = coupled.coordinates * BOHR_TO_ANGSTROM
    
    potential = (coupled.planar - bottom.planar - top.planar) * RY_TO_EV

    # coordinates = coupled.potential_averaged.coordinates
    # potential = coupled.potential_averaged.planar - bottom.potential_averaged.planar - top.potential_averaged.planar
    
    fig, ax = plt.subplots()
    
    ax.plot(z, potential, label=r"U_{coupled layers} - U_{bottom layer} - U_{top layer}")
    
    ax.set_xlabel(r"z (($\AA$))")
    ax.set_ylabel("Potential Energy (eV)")
    ax.set_title(f"Subtracted Potential, E-field={field}au")
    ax.legend(loc="upper right")


    plt.savefig(POTENTIAL_DIR / f"Subtracted_Potential_{field}au.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    

def potential_individual(results:Results, field):
    print("\nPLOTTING INDIVIDUAL POTENTIALS")

    for system in results.__dict__.values():
        system:System

        coordinates = system.potential_averaged.coordinates * BOHR_TO_ANGSTROM
        potential = system.potential_averaged.planar * RY_TO_EV
        
        print(f"\nPLOTTING {system.name} POTENTIAL")

        fig, ax = plt.subplots()

        ax.plot(coordinates, potential, linewidth=0.75, color='red')
        ax.set_title(f"1D {system.name.capitalize()} Layer Potential V(z), E-field={field}au")
        ax.set_xlabel(r"z ($\AA$)")
        ax.set_ylabel("Potential Energy (eV)")
        
        plt.savefig(POTENTIAL_DIR / f"Potential_{system.name.capitalize()}_Layer_{field}au.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def charge_density_subtracted(results:Results, field):
    print("\nPLOTTING SUBTRACTED CHARGE DENSITY")

    coupled = results.coupled.charge_density_averaged
    bottom = results.bottom.charge_density_averaged
    top = results.top.charge_density_averaged
    
    z = coupled.coordinates * BOHR_TO_ANGSTROM
    
    charge_density = (coupled.planar - bottom.planar - top.planar) / (BOHR_TO_ANGSTROM)**3
    
    fig, ax = plt.subplots()
    
    ax.plot(z, charge_density)
    
    ax.set_xlabel(r"z (($\AA$))")
    ax.set_ylabel(r"Charge Density ($e/{\AA}^3)$")
    ax.set_title(f"Subtracted Charge Density, E-field={field}au")

    plt.savefig(CHARGE_DENSITY_DIR / f"Subtracted_Charge_Density_{field}au.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_potential(results, field):
    
    potential_subtracted(results, field)
    potential_individual(results, field)


def plot_charge_density(result, field):
    
    charge_density_subtracted(result, field)
    

def plot_dos(results, field):
    
    #individual_dos(results, field)
    #dos_comparison(results, field)
    fermi_aligned(results, field)
