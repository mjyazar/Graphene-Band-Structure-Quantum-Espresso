import matplotlib.pyplot as plt
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


def plot_dos(energy, dos, fermi_energy, name, window=(-15, 15)):

    fig, ax = plt.subplots()
    
    ax.plot(energy - fermi_energy, dos, linewidth=0.75, color='red')
    ax.axvline(0, linestyle="--")
    ax.set_xlabel(r"$E - E_{Fermi}$ (eV)")
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.set_xlim(window)
    # ax.set_ylim()
    
    ax.fill_between(energy - fermi_energy, 0, dos, where=(energy - fermi_energy < 0), facecolor='red', alpha=0.25)
    ax.text(0, 2.5, "Fermi energy", fontsize=16, rotation=90)

    fig.tight_layout()
    plt.savefig(FIG_DIR / f"{name.capitalize()} Layer DOS.png", dpi=300)
    plt.close(fig)
