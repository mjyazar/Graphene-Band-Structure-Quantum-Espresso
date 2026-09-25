from config import *
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
#MONOLAYER = OUT_DIR / "monolayer"
#BILAYER = OUT_DIR / "bilayer"

# Ensure folders exist
OUT_DIR.mkdir(exist_ok=True)
#MONOLAYER.mkdir(exist_ok=True)
#BILAYER.mkdir(exist_ok=True)


def run_calculations(structure, path, eamp, xc):
    
    print(f"\n\n{path.name.upper()} LAYER(S) COMPUTATIONS")
    
    if RUN_QE:
        scf = pw.scf(structure, path, eamp, xc)
        nscf = pw.nscf(structure, path, eamp, xc)
    else: 
        scf = pw._read_output(path / "scf.pwo")
        nscf = pw._read_output(path / "nscf.pwo")
    
    fermi_energy = nscf.calc.get_fermi_level()
    
    if RUN_POTENTIAL:
        potential, potential_pp_path = pp.potential(path)
        # potential_averaged = average.potential(path, potential_pp_path)
    else:
        # potential, potential_pp_path = pp._read_output(path / "potential.cube", 6), path / "data" / f"potential.pp.dat"        
        data, atoms, xy_averaged, z = pp._read_output(path / "potential.cube", 5)
        potential = Potential(data=data, atoms=atoms, averaged=xy_averaged, z=z)
        
    if RUN_CHARGE_DENSITY:
        charge, charge_pp_path = pp.charge_density(path)
    else:
        # charge, charge_pp_path = pp._read_output(path / "charge.cube", 6), path / "data" / f"charge.pp.dat"
        data, atoms, xy_averaged, z = pp._read_output(path / "charge.cube", 5)
        charge = ChargeDensity(data=data, atoms=atoms, averaged=xy_averaged, z=z)
    
    if RUN_DOS:
        dos_ = dos.calculate(path, fermi_energy)
    else:
        dos_ = dos._read_output(path / "dos.out")
    
    if RUN_LDOS:
        ldos, ldos_pp_path, energies = pp.ldos(path, fermi_energy)
        # ldos_averaged = average.ldos(path, ldos_pp_path, energies)
    else:
        ldos = pp._read_ldos(path / "data" / f"ldos.pp.dat")
        
    # ldos_averaged = average.ldos(path, ldos_pp_path, energies)
        

    results = System(name=f"{path.name}",
                            atoms=structure,
                            path=path,
                            fermi_energy=fermi_energy,
                            potential=potential,
                            charge_density=charge,
                            dos = dos_,
                            ldos = ldos)
    
    return results


def main():
    graphene = GrapheneStructure()
    
    energies = [0]
    
    for eamp in energies:
        
        path = OUT_DIR /  f"field_{str(eamp)}"
        
        print(f"\nE-field = {str(eamp)}au")
        print("-" * 30)

        print("CREATING GRAPHENE BILAYERS")
        bilayer = graphene.bilayer()

        if RUN_QE:
            print("RELAXING COUPLED BILAYER")
            relaxed_coupled = pw.relax(bilayer, path, eamp)
        
        else:
            relaxed_coupled = pw._read_output(path / "relax.pwo")
        
        print("\nEXTRACTING FROZEN LAYERS")
        isolated_bottom, isolated_top = graphene.isolate_bilayer(relaxed_coupled)
        
        results = {}
        for xc in XC:
            if xc is None:
                label = "d3"
            else:
                label = xc.split("-")[0]
            
            PATH_COUPLED = path / (f"PBE-{label}" if xc is None else xc) / "coupled" 
            PATH_BOTTOM = path / (f"PBE-{label}" if xc is None else xc) / "bottom"
            PATH_TOP = path / (f"PBE-{label}" if xc is None else xc) / "top"
            
            coupled = run_calculations(relaxed_coupled, PATH_COUPLED, eamp, xc)
            bottom = run_calculations(isolated_bottom, PATH_BOTTOM, eamp, xc)
            top = run_calculations(isolated_top, PATH_TOP, eamp, xc)
            
            result = Results(coupled, bottom, top)
            results[label] = result
            
            plotter.plot_potential(result, eamp, label)
            plotter.plot_charge_density(result, eamp, label)
            plotter.plot_dos(result, eamp, label)
            plotter.plot_ldos(result, eamp, label)

        plotter.plot_vdw_charge_difference(results, eamp)


if __name__ == "__main__":
    main()
