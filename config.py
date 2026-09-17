
BOHR_TO_ANGSTROM = 0.529177210903
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
pw.x
"""
ECUTRHO = 400.0
CONV_THRESHOLD = 1.0e-8
DEGAUSS = 0.01
SMEARING = "gauss"
SCF_EXTRA_BANDS = 6  # number of unoccupied bands to run the calculations for nzcf
NSCF_EXTRA_BANDS_PER_ATOM = 4  # number of unoccupied bands to run the calculations for nscf
FORCE_CONVERGENCE_THRESHOLD = 0.001
BANDPATH = 'GMKG'
ecutwfc = 60.0
KGRID = (15, 15, 1)
KGRID_DENSE = (30, 30, 1)


"""
pp.x
"""
DEGAUSS_LDOS = 0.05 * EV_TO_RY  # eV
WINDOW_LDOS = (-10, 10)
DELTA_E_LDOS = 0.1


"""
dos.x
"""
DEGAUSS = 0.01
WINDOW_DOS = (-10, 10)
DELTA_E_DOS = 0.02

# RY_TO_EV = 


"""
PLOTTER
"""

WINDOW = (-10, 10)
