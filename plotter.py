import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import TwoSlopeNorm
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
DOS_DIR = FIG_DIR / "dos"
LDOS_DIR = FIG_DIR / "ldos"

FIG_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)
POTENTIAL_DIR.mkdir(exist_ok=True)
CHARGE_DENSITY_DIR.mkdir(exist_ok=True)
DOS_DIR.mkdir(exist_ok=True)
LDOS_DIR.mkdir(exist_ok=True)


def plot_band_structure(bandpath, energies, name):
    """
    Plot and save band structure
    """
    print("PLOTTING BAND STRUCTURE")
    
    band_structure = BandStructure(path=bandpath, energies=energies, reference=0.0)  # reference is now zero after shifting
    
    fig, ax = plt.subplots()

    ax = band_structure.plot()
    assert ax is not None  # ensures BandStructure doesn't return None to prevent an error
    ax.set_ylim(-10, 10)
    ax.set_ylabel("Energy - $E_F$ (eV)")
    ax.set_title(f"Graphene {name.capitalize()} Band Structure")

    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{name.capitalize()} Graphene Band Structure.png", dpi=300)
    plt.close(fig)


def _get_vacuum_level(V, z):
        
    assert V.shape == z.shape
    
    # old implementation
    # atomic_positions = structure.get_positions()[:, 2]
    # z_bottom = atomic_positions.min()
    # z_top = atomic_positions.max()
    # vacuum_mask = (z < z_bottom) | (z > z_top)
    # positions = V[vacuum_mask]
    
    dV_dz = np.gradient(V, z)
    tolerance = VACUUM_LEVEL_TOLERANCE
    vacuum_level = V[np.abs(dV_dz) < tolerance]

    # if tolerance is too low prevent code crash by increasing magnitude
    while vacuum_level.size == 0:
        tolerance *= 10
        vacuum_level = V[np.abs(dV_dz) < tolerance]
    
    vacuum_level = np.min(vacuum_level)
    
    return vacuum_level


def potential_subtracted(results: Results, field, vdw_scheme):
    print("PLOTTING SUBTRACTED POTENTIAL")
    
    coupled = results.coupled.potential
    bottom = results.bottom.potential
    top = results.top.potential

    np.testing.assert_allclose(coupled.z, bottom.z)
    np.testing.assert_allclose(coupled.z, top.z)
    
    z = coupled.z
    
    potential = (coupled.averaged - bottom.averaged - top.averaged) * RY_TO_EV

    # coordinates = coupled.potential_averaged.coordinates
    # potential = coupled.potential_averaged.planar - bottom.potential_averaged.planar - top.potential_averaged.planar
    
    fig, ax = plt.subplots()
    
    ax.plot(z, potential, color='black', label=r"$U_{coupled layers} - U_{bottom layer} - U_{top layer}$")
    
    ax.set_xlabel(r"z ($\AA$)")
    ax.set_ylabel("Potential Energy (eV)")
    ax.set_title(f"Subtracted Potential, E-field={field}au, Correction: {vdw_scheme}")
    ax.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(POTENTIAL_DIR / f"Subtracted_Potential_{field}au_{vdw_scheme}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    

def potential_individual(results: Results, field, vdw_scheme):
    print("PLOTTING INDIVIDUAL POTENTIALS")
    
    # iterate through coupled, top, and bottom
    for system in results.__dict__.values():
        system:System

        z = system.potential.z
        potential = system.potential.averaged * RY_TO_EV
        
        print(f"PLOTTING {system.name} POTENTIAL")

        fig, ax = plt.subplots()

        ax.plot(z, potential, linewidth=0.75, color='black')
        ax.set_title(f"{system.name.capitalize()} Layer Potential V(z), E-field={field}au, Correction: {vdw_scheme}")
        ax.set_xlabel(r"z ($\AA$)")
        ax.set_ylabel("Potential Energy (eV)")
        
        if field == 0:
            vacuum_level = _get_vacuum_level(potential, z)
            ax.axhline(vacuum_level, color="red", linestyle="--", linewidth=0.8, label="Vacuum Level")

        plt.tight_layout()
        plt.savefig(POTENTIAL_DIR / f"Potential_{system.name.capitalize()}_Layer_{field}au_{vdw_scheme}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def charge_density_subtracted(results: Results, field, vdw_scheme):
    print("PLOTTING SUBTRACTED CHARGE DENSITY")

    coupled = results.coupled.charge_density
    bottom = results.bottom.charge_density
    top = results.top.charge_density
    
    z = coupled.z
    
    charge_density = (coupled.averaged - bottom.averaged - top.averaged) / (BOHR_TO_ANGSTROM)**3
    
    fig, ax = plt.subplots()
    
    ax.plot(z, charge_density, color="black")
    
    ax.set_xlabel(r"z ($\AA$)")
    ax.set_ylabel(r"Charge Density ($e/{\AA}^3)$")  # averaged but not multiplied by area
    ax.set_title(f"Subtracted Charge Density, E-field={field}au, Correction: {vdw_scheme}")

    plt.tight_layout()
    plt.savefig(CHARGE_DENSITY_DIR / f"Subtracted_Charge_Density_{field}au_{vdw_scheme}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_vdw_charge_difference(results, field):
    """
    Plot the plane averaged charge density difference between dn_c09 and dn_pbe, 
    where each is the subtracted charge density i.e. coupled - top - bottom
    """
    
    
    print(f"PLOTTING VDW Charge Density Difference")
    pbe: Results = results["d3"]
    c09: Results = results["c09"]
    
    np.testing.assert_allclose(pbe.coupled.charge_density.z, pbe.bottom.charge_density.z)
    np.testing.assert_allclose(pbe.coupled.charge_density.z, pbe.top.charge_density.z)
    np.testing.assert_allclose(c09.coupled.charge_density.z, c09.bottom.charge_density.z)
    np.testing.assert_allclose(c09.coupled.charge_density.z, c09.top.charge_density.z)
    np.testing.assert_allclose(pbe.coupled.charge_density.z, c09.coupled.charge_density.z)
    
    z = pbe.coupled.charge_density.z
    
    dn_pbe = pbe.coupled.charge_density.averaged - pbe.bottom.charge_density.averaged - pbe.top.charge_density.averaged
    dn_c09 = c09.coupled.charge_density.averaged - c09.bottom.charge_density.averaged - c09.top.charge_density.averaged
    
    cell = c09.coupled.atoms.cell.array
    area = np.linalg.norm(np.cross(cell[0], cell[1]))  # cross product of a1 and a2 (i.e. in-plane vecors)
    
    dn_vdw = ((dn_c09 - dn_pbe) / BOHR_TO_ANGSTROM**3) * area

    fig, ax = plt.subplots(figsize=(4, 7))
    
    ax.plot(dn_vdw, z, color="black")
    ax.axvline(0, color="blue", linestyle="--", linewidth=0.8)

    # print("\nc09 Coupled Atoms Positions: \n", c09.coupled.atoms.positions)
    positions_z = c09.coupled.atoms.positions[:, 2]
    
    bottom_z = positions_z.min()
    top_z = positions_z.max()
    
    ax.axhline(bottom_z, color="blue", linestyle=":", linewidth=0.8)
    ax.axhline(top_z, color="red", linestyle=":", linewidth=0.8)

    ax.set_ylim(0, 20)
    
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    
    ax.tick_params(axis="both", direction="in")
    
    ax.set_xlabel(r"$\Delta n_{vdW}(z)$ ($e/{\AA})$")
    ax.set_ylabel(r"z ($\AA$)")
    ax.set_title(rf"vdW Charge Density $\Delta n_{{c09}}(z)-\Delta n_{{PBE + d3}}(z)$, E-field={field}au")

    plt.tight_layout()
    plt.savefig(CHARGE_DENSITY_DIR / f"vdW_Charge_Density_{field}au.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def dos_subtracted(results: Results, field, vdw_scheme, delta_e=DELTA_E_DOS):
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
    ax.plot(grid, dos_subtracted, color="black", label=r"$DOS_{coupled} - DOS_{top} - DOS_{bottom}$")

    ax.set_title(f"Fermi-Aligned DOS Subtracted, E-field={field}au, Correction: {vdw_scheme}")
    ax.set_xlabel(r"$Energy$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.legend()

    plt.tight_layout()
    plt.savefig(DOS_DIR / f"DOS_Subtracted_{field}au_{vdw_scheme}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def dos_individual(results: Results, field, vdw_scheme, delta_e=LDOS_GRID_DELTA_E):
    print(f"PLOTTING FERMI-ALIGNED INDIVIDUAL DOS")

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
    
    plots = [(dos_coupled_grid, "Coupled Layers", "black", "-"), (dos_bottom_grid, "Bottom Layer", "red", "--"), (dos_top_grid, "Top Layer", "green", ":")]
    
    fig, ax = plt.subplots()

    for plot, label, color, linestyle in plots:
        ax.plot(grid, plot, color=color, label=f"{label}", linestyle=linestyle)

    ax.set_title(f"Fermi-Aligned DOS Subtracted, E-field={field}au, Correction: {vdw_scheme}")
    ax.set_xlabel(r"$Energy$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(DOS_DIR / f"DOS_Individual_{field}au_{vdw_scheme}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)



def ldos_individual(results: Results, field, vdw_scheme):
    print(f"PLOTTING FERMI-ALIGNED INDIDIVDUAL LDOS")
        
    coupled = results.coupled
    bottom = results.bottom
    top = results.top
    
    np.testing.assert_allclose(coupled.ldos.z, bottom.ldos.z)
    np.testing.assert_allclose(coupled.ldos.z, top.ldos.z)
    
    z = coupled.ldos.z
    V = coupled.potential.averaged * RY_TO_EV - coupled.fermi_energy
    
    plots = [(coupled, "coupled"), (bottom, "bottom"), (top, "top")]

    for layer, label in plots:
        fig, ax = plt.subplots()
        
        E = layer.ldos.energies

        limit = np.max(np.abs(layer.ldos.averaged))
        symmetric_norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
        mesh = ax.pcolormesh(z, E, layer.ldos.averaged, cmap="seismic", shading="auto", norm=symmetric_norm)
        
        fig.colorbar(mesh, ax=ax, label=r"$\Delta LDOS$")

        ax.plot(layer.potential.z, V, linewidth=0.75, color='black')
        ax.set_title(f"Fermi-Aligned LDOS, E-field={field}au, Correction: {vdw_scheme}")
        ax.set_xlabel(r"$z$ ($\AA$)")
        ax.set_ylabel(r"$E-E_F$ (eV)")
        
        plt.tight_layout()
        plt.savefig(LDOS_DIR / f"LDOS_{label}_{field}au_{vdw_scheme}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def ldos_subtracted(results: Results, field, vdw_scheme):
    print(f"PLOTTING FERMI-ALIGNED SUBTRACTED LDOS")

    coupled = results.coupled.ldos
    bottom = results.bottom.ldos
    top = results.top.ldos
    
    z = coupled.z
    
    np.testing.assert_allclose(coupled.z, bottom.z)
    np.testing.assert_allclose(coupled.z, top.z)
    
    E = coupled.energies
    
    ldos_subtracted = coupled.averaged - bottom.averaged - top.averaged
    
    fig, ax = plt.subplots()
    
    limit = np.max(np.abs(ldos_subtracted))
    symmetric_norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    mesh = ax.pcolormesh(z, E, ldos_subtracted, cmap="seismic", shading="auto", norm=symmetric_norm)
    
    fig.colorbar(mesh, ax=ax, label=r"$\Delta LDOS$")
    
    ax.set_title(f"Fermi-Aligned Subtracted LDOS, E-field={field}au, Correction: {vdw_scheme}")
    ax.set_xlabel(r"$z$ ($\AA$)")
    ax.set_ylabel(r"$E-E_F$ (eV)")
    
    plt.tight_layout()
    plt.savefig(LDOS_DIR / f"LDOS_Subtracted_{field}au_{vdw_scheme}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_potential(results, field, vdw_scheme):
    print("")
    potential_subtracted(results, field, vdw_scheme)
    potential_individual(results, field, vdw_scheme)


def plot_charge_density(results, field, vdw_scheme):
    charge_density_subtracted(results, field, vdw_scheme)


def plot_dos(results, field, vdw_scheme):
    #individual_dos(results, field)
    #dos_comparison(results, field)
    dos_subtracted(results, field, vdw_scheme)


def plot_ldos(results, field, vdw_scheme):
    ldos_subtracted(results, field, vdw_scheme)
    ldos_individual(results, field, vdw_scheme)
