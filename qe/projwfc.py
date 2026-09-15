import qe.runner as runner


WINDOW_PROJWFC = (-10, 10)


def write_input(input_path, outdir, data_path, prefix, fermi_energy, window=WINDOW_PROJWFC):
    """
    emin / emax -> E_F ± 8-10 eV
    deltae - resolution
    degauss - broadening of Gaussian wavefunctions
    
    """
    
    with open(input_path, "w") as input_file:
        input_file.write(f"""&PROJWFC
                         prefix = pw.x
                         outdir = "{outdir}"
                         ngauss = 0
                         degauss = 0.05
                         emin = -10
                         emax = 10
                         deltaE = 0.01  ! energy grid step (eV)
                         kresolveddos = .true
                         filproj = "{data_path}"
                         """)


def read_output(path):

    # = np.loadtxt(path, unpack=True)
    
    pass


def calculate(path, fermi_energy):
    
    path.mkdir(parents=True, exist_ok=True)
    
    input_path = path / "projwfc.in"
    output_path = path / "projwfc.out"
    data_path = path / "projwfc.dat"

    print(f"\nCREATING {input_path.name}")
    write_input(input_path, path / "data", data_path, path.name, fermi_energy)
    
    print(f"RUNNING projwfc.x WITH {input_path.name}")
    runner.run("projwfc.x", input_path, output_path)
    
    print(f"READING {data_path.name}")

    return read_output(data_path)
