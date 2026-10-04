from pathlib import Path
ROOT = Path(__file__).resolve().parent


"""
RUN
False for efficient run if data from previous run present, True if running for the first time or want to create new files with new parameters
"""
RUN_QE = True
RUN_POTENTIAL = RUN_QE
RUN_CHARGE_DENSITY = RUN_QE
RUN_DOS = RUN_QE
RUN_LDOS = RUN_QE
RUN_CONVERGENCE = False


"""
RUNNER
"""
NPROC = 16
NK = 16  # split kgrid computations into n pools with


"""
vdW SCHEME
"""
VDW_CORR = "grimme-d3"
XC = [None, "vdw-df2-c09"]


"""
CONSTANTS AND CONVERSIONS
"""
BOHR_TO_ANGSTROM = 0.529177210903  # Bohr radius
RY_TO_EV = 13.605693122994
EV_TO_RY = 1 / RY_TO_EV


"""
GRAPHENE
"""
LATTICE_CONSTANT = 2.46
VACUUM = 10.0  # QE requires 2D Coulomb truncation of the cell to have min z-length ~10.58 A
INTERLAYER_DISTANCE = 3.35  # for bilayer (angstrom)
STACKING = "AB"


"""
Carbon pseudopotential
from https://sssp.materialscloud.org/pseudopotentials/PBE/efficiency
"""
PSEUDO = "C.upf"
# PSEUDO = "C.pbe-n-kjpaw_psl.1.0.0.UPF"


"""
pw.x
"""
ECUTWFC = 90.0
ECUTRHO = 4 * ECUTWFC
CONV_THRESHOLD = 1.0e-8
DEGAUSS_PW = 0.01
SMEARING = "gauss"
SCF_EXTRA_BANDS = 6  # number of unoccupied bands to run the calculations for nscf
NSCF_NBND = 30
NSCF_EXTRA_BANDS_PER_ATOM = 4  # number of unoccupied bands to run the calculations for nscf
FORCE_CONVERGENCE_THRESHOLD = 0.001
BANDPATH = 'GMKG'
KGRID = (15, 15, 1)  # scf
KGRID_DENSE = (75, 75, 1)  # nscf


"""
pp.x
"""
DEGAUSS_LDOS = 0.2
WINDOW_LDOS = (-15, 15)
DELTA_E_LDOS = 0.02
VACUUM_LEVEL_TOLERANCE = 0.0000001  # do not go lower -> error


"""
dos.x
"""
DEGAUSS_DOS = DEGAUSS_LDOS * EV_TO_RY
WINDOW_DOS = (-15, 15)
DELTA_E_DOS = 0.01


"""
PLOTTER
"""
WINDOW = (-15, 15)
LDOS_GRID_DELTA_E = 0.01
FIGSIZE_CHARGE_DENSITY = (5, 8)


"""
CONVERGENCE
"""
ECUTWFC_BOUND = 110
SCF_BOUND = 30
NSCF_UPPER_BOUND = 100
NSCF_LOWER_BOUND = 24  # lower bound
