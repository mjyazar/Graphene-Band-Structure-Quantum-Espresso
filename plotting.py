import matplotlib.pyplot as plt
from pathlib import Path
from ase.spectrum.band_structure import get_band_structure, BandStructure

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
    
    ax = band_structure.plot()
    assert ax is not None  # ensures BandStructure doesn't return None to prevent an error
    ax.set_ylim(-10, 10)
    ax.set_ylabel("Energy - $E_F$ (eV)")
    ax.set_title(f"Graphene {name.capitalize()} Band Structure")
    
    plt.savefig(FIG_DIR / f"{name.capitalize()} Graphene Band Structure.png", dpi=300)
    plt.close()


def plot_dos(energy, dos, fermi_energy, name):
    
    plt.figure(figsize = (12, 6))
    plt.plot(energy, dos, linewidth=0.75, color='red')
    plt.yticks([])
    plt.xlabel('Energy (eV)')
    plt.ylabel('DOS')
    plt.axvline(x=fermi_energy, linewidth=0.5, color='k', linestyle=(0, (8, 10)))
    plt.xlim(-15, 15)
    plt.ylim(0, )
    plt.fill_between(energy, 0, dos, where=(energy < fermi_energy), facecolor='red', alpha=0.25)
    plt.text(0, 2.5, 'Fermi energy', fontsize= 16, rotation=90)
    plt.show()
