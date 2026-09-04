import numpy as np
from ase import Atoms
# from ase.constraints import FixCartesian

LATTICE_CONSTANT = 2.46
VACUUM = 10.0  # QE requires 2D Coulomb truncation of the cell to have min z-length ~10.58 A

INTERLAYER_DISTANCE = 3.35  # for bilayer (angstrom)
STACKING = "AB"


class GrapheneStructure:
    
    def __init__(self, lattice_const=LATTICE_CONSTANT, vacuum=VACUUM):
        self.lattice_const = lattice_const
        self.vacuum = vacuum


    def cell(self, dz=0.):
        """
        Construct unit cell lattice vectors
        https://mc2d.materialscloud.org/details/mc2d-71
        
        dz: thickness of the layer
        """
        
        a = self.lattice_const
                
        return np.array([[a, 0, 0], 
                         [-a/2, a*np.sqrt(3)/2, 0],
                         [0, 0, 2*self.vacuum + dz]])
    
    
    def monolayer(self):
        """
        Create Monolayer Graphene as Atoms object: honeycomb, 2 atoms/unit cell
        """
        
        # atomic positions in units of the unit cell 
        atomic_positions = np.array([[0, 0, 0.5],      # atom 1
                                     [2/3, 1/3, 0.5]]) # atom 2
        
        return Atoms(symbols="C2", cell=self.cell(), scaled_positions=atomic_positions, pbc=True)
    
    
    def bilayer(self, interlayer_dist=INTERLAYER_DISTANCE, stacking=STACKING):
        """
        Create Bilayer Graphene as an Atoms object
        """
        
        # calculating z coordinates in units of unit cell (i.e. 0 to 1)
        dz = interlayer_dist / (self.vacuum * 2 + interlayer_dist)
        
        # centre the bilayer around around z=0.5
        z1 = 0.5 - dz/2
        z2 = 0.5 + dz/2
            
        if stacking == "AB":
            atomic_positions = np.array([[0, 0, z1], # atom 1
                                   [2/3, 1/3, z1],   # atom 2
                                   [2/3, 1/3, z2],    # atom 3
                                   [1/3, 2/3, z2]])   # atom 4

        elif stacking == "BA":
            atomic_positions = np.array([[0, 0, z1],
                                         [2/3, 1/3, z1],
                                         [1/3, 2/3, z2],
                                         [0, 0, z2]])

        # must be AA otherwise
        else:
            atomic_positions = np.array([[0, 0, z1],
                                         [2/3, 1/3, z1],
                                         [0, 0, z2],
                                         [2/3, 1/3, z2]])
            
        return Atoms(symbols="C4", cell=self.cell(interlayer_dist), scaled_positions=atomic_positions, pbc=True)

    
    def isolate_bilayer(self, bilayer, layer):
        
        # returns atoms 1 and 2 - in bottom layer
        z = bilayer.positions[:, 2]
        
        mid = (z.min() + z.max()) / 2
        
        bottom = bilayer[np.where(z < mid)[0]]
        top = bilayer[np.where(z >= mid)[0]]

        return bottom, top
