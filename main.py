from bgw_qe_tools.read_qe_xml import traj_from_qe_xml
from ase.visualize import view
from matplotlib import pyplot as plt
from graphene import GrapheneStructure
from results import *
import qe.pw as pw
import qe.dos as dos
import qe.pp as pp
import qe.average as average
import convergence
import plotter

from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Paths
OUT_DIR = ROOT / "outputs"
# MONOLAYER = OUT_DIR / "monolayer"
# BILAYER = OUT_DIR / "bilayer"

# Ensure folders exist
OUT_DIR.mkdir(exist_ok=True)
# MONOLAYER.mkdir(exist_ok=True)
# BILAYER.mkdir(exist_ok=True)

# RUN_QE may be False if script already ran and want to work with existing files
# True if running for the first time or want to create new files with new parameters
RUN_QE = True
RUN_CONVERGENCE = True


def run_calculations(structure, path, eamp):
    print(f"\n{path.name.upper()} LAYERS COMPUTATIONS")

    scf = pw.scf(structure, path, eamp)
    nscf = pw.nscf(structure, path, eamp)

    fermi_energy = nscf.calc.get_fermi_level()

    potential_raw, potential_pp_path = pp.potential(path)
    potential_averaged = average.potential(path, potential_pp_path)

    plt.contourf(potential_raw.data[15, :, :])
    plt.show()

    charge_raw, charge_pp_path = pp.charge_density(path)
    charge_averaged = average.charge_density(path, charge_pp_path)

    dos_ = dos.calculate(path, fermi_energy)

    ldos_pp_path = pp.ldos(path, fermi_energy)
    ldos_averaged = average.ldos(path, ldos_pp_path)

    results = System(name=f"{path.name}",
                     atoms=structure,
                     path=path,
                     fermi_energy=fermi_energy,

                     potential_raw=potential_raw,
                     potential_averaged=potential_averaged,
                     charge_density_raw=charge_raw,
                     charge_density_averaged=charge_averaged,
                     dos=dos_,
                     ldos=ldos_averaged)

    return results


def main():
    graphene = GrapheneStructure()

    energies = [0]

    for eamp in energies:

        path = OUT_DIR / f"field_{str(eamp)}"

        PATH_COUPLED = path / "coupled"
        PATH_BOTTOM = path / "bottom"
        PATH_TOP = path / "top"

        print(f"\nE-field = {str(eamp)}au")
        print("-" * 30)

        print("CREATING GRAPHENE BILAYERS")
        bilayer = graphene.bilayer()

        print("RELAXING COUPLED BILAYER")
        reuse_data = True
        if reuse_data:
            relaxed_coupled = pw.relax(bilayer, PATH_COUPLED, eamp)
            atoms, ecut, atoms_list = traj_from_qe_xml(PATH_COUPLED / "data" / "coupled.xml")
            if len(atoms_list) > 0:
                view(atoms_list)
            relaxed_coupled = atoms
        else:
            relaxed_coupled = pw.relax(bilayer, PATH_COUPLED, eamp)

        print("\nEXTRACTING FROZEN LAYERS")
        isolated_bottom, isolated_top = graphene.isolate_bilayer(relaxed_coupled)

        coupled = run_calculations(relaxed_coupled, PATH_COUPLED, eamp)
        bottom = run_calculations(isolated_bottom, PATH_BOTTOM, eamp)
        top = run_calculations(isolated_top, PATH_TOP, eamp)

        results = Results(coupled, bottom, top)

        plotter.plot_potential(results, eamp)
        plotter.plot_charge_density(results, eamp)
        plotter.plot_dos(results, eamp)
        plotter.plot_ldos(results, eamp)


if __name__ == "__main__":
    main()
