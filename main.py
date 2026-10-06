from config import *
from graphene import GrapheneStructure
from results import *
import qe.pw as pw
import qe.dos as dos
import qe.pp as pp
import qe.average as average
import convergence
import plotter as plot

# Paths
OUT_DIR = ROOT / "outputs"
#MONOLAYER = OUT_DIR / "monolayer"
#BILAYER = OUT_DIR / "bilayer"

# Ensure folders exist
OUT_DIR.mkdir(exist_ok=True)
#MONOLAYER.mkdir(exist_ok=True)
#BILAYER.mkdir(exist_ok=True)


def run_calculations(structure, path, eamp, xc, label):
    
    print(f"\n{label} | {path.name.upper()} LAYER(S) COMPUTATIONS")
    print("-" * 35)

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
        data, atoms, xy_averaged, z = pp._read_output(path / "potential.xsf", 5)
        potential = Potential(data=data, atoms=atoms, averaged=xy_averaged, z=z)


    if RUN_CHARGE_DENSITY:
        charge, charge_pp_path = pp.charge_density(path)
    else:
        # charge, charge_pp_path = pp._read_output(path / "charge.cube", 6), path / "data" / f"charge.pp.dat"
        data, atoms, xy_averaged, z = pp._read_output(path / "charge.xsf", 5)
        charge = ChargeDensity(data=data, atoms=atoms, averaged=xy_averaged, z=z)
    
    if RUN_DOS:
        dos_ = dos.calculate(path, fermi_energy)
    else:
        dos_ = dos._read_output(path / "dos.out")
    
    if RUN_LDOS:
        ldos, ldos_pp_path, energies = pp.ldos(path, fermi_energy)
        for f in (path / "data").glob("*.save/wfc*.dat"):
            f.unlink()
        # ldos_averaged = average.ldos(path, ldos_pp_path, energies)
    else:
        ldos = pp._load_ldos(path)
        
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
    
    #energies = [0, 0.001, 0.005]
    energies = [0, 0.001, 0.005]
    
    for eamp in energies:
        
        path = OUT_DIR /  f"field_{str(eamp)}"

        print()
        print("*" * 35)
        print(f"E-field = {str(eamp)}au")
        print("*" * 35)

        print("\nCREATING GRAPHENE BILAYERS")
        bilayer = graphene.bilayer()

        if RUN_QE:
            print("RELAXING COUPLED BILAYER")
            relaxed_coupled = pw.relax(bilayer, path, eamp)
        
        else:
            relaxed_coupled = pw._read_output(path / "relax.pwo")
        
        print("\nEXTRACTING FROZEN LAYERS")
        isolated_bottom, isolated_top = graphene.isolate_bilayer(relaxed_coupled)
        
        if RUN_CONVERGENCE:
            convergence.test_kgrid(relaxed_coupled, eamp, "scf")
            convergence.test_kgrid(relaxed_coupled, eamp, "nscf")
            convergence.test_ecutwfc(relaxed_coupled, eamp)
        
        # store pbe+d3 results as d3 and vdw results as c09
        results = {}
        for xc in XC:
            if xc is None:
                label = "d3"
            else:
                label = xc.split("-")[-1]
            
            PATH_COUPLED = path / (f"PBE-{label}" if xc is None else xc) / "coupled" 
            PATH_BOTTOM = path / (f"PBE-{label}" if xc is None else xc) / "bottom"
            PATH_TOP = path / (f"PBE-{label}" if xc is None else xc) / "top"
            
            coupled = run_calculations(relaxed_coupled, PATH_COUPLED, eamp, xc, label)
            bottom = run_calculations(isolated_bottom, PATH_BOTTOM, eamp, xc, label)
            top = run_calculations(isolated_top, PATH_TOP, eamp, xc, label)
            
            result = Results(coupled, bottom, top)
            results[label] = result

            plot.potential_individual(result, eamp, label)
            plot.potential_subtracted(result, eamp, label)         
            plot.charge_density_subtracted(result, eamp, label)
            plot.dos_subtracted(result, eamp, label)
            plot.dos_individual(result, eamp, label)
            plot.ldos_subtracted(result, eamp, label)
            plot.ldos_individual(result, eamp, label)
        
        pp.vdw_charge_difference(path, eamp)  # creates corresponding xsf file
        plot.plot_vdw_charge_difference(results, eamp)


if __name__ == "__main__":
    main()
