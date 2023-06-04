from math import pi, cos, sin
import numpy as np
from base import Molecule
import sys

"""
I'm not sure lattice dimensions are properly handled here if I want several
flakes. Maybe one needs to make several flakes separately and then put them
together into the Lammps data files.

This will require some tinkering though, and would be nicer if one cold do
in this code itself.

One should alse be able to rotate flakes, but this is probably done in the
crystal class.
"""

class Rectangle_Graphene(Molecule):
    def __init__(self, x_length, y_length, cc_bond_length=1.4148, 
                 ch_bond_length = 1.077, layer_gap=3.4827):
        # zigzag edges run along the x axis
        # armchair edges run along the y axis
        # config = self.crystal_params()

        self.CC = cc_bond_length          #Bond length C-C
        self.CH = ch_bond_length         #Nond length C-H edge
        self.layer_gap = layer_gap #Interlayer distance

        self.x_length = x_length
        self.y_length = y_length
        x_unit_length = 2.0 * self.CC * cos(pi / 6.0)
        y_unit_length = 1.5 * self.CC



        #Units along x-axis
        # 
        #         
        #         |       These are our base units with three 
        #         C       Carbons per row. The x-unit-length of each unit 
        #       /   \     and y-unit-length, seems not to be in the class.
        #      /     \    only the number of units in x and y directions
        #     C       C    
        #     |       |
        #   <-- x_unit_length -->    



        # number of rows of hexagons along each axis
        self.x = int((x_length - self.CC * cos(pi / 6.0)) / x_unit_length)
        self.y = int((y_length - self.CC * sin(pi / 6.0)) / y_unit_length)
        # mandate odd number of rows of hexagons for symmetry
        if not self.y % 2:
            self.y -= 1

        if (self.x <= 0) or (self.y <= 0):
            raise Exception("Dimensions too small")

        carbons_per_row = 1 + self.x * 2
        self.n_Cs = carbons_per_row * (self.y + 1)
        
        

        hydrogens_x = 2 * self.x
        hydrogens_y = 2 * (self.y + 1)
        self.n_Hs = hydrogens_x + hydrogens_y

        self.natoms = self.n_Cs + self.n_Hs

    #This one ought to have an additional argument.....
    def cell_shape(self):
        a = 0 + self.x_length    #Why the hard-coded 100 Å?
        b = 0 + self.y_length    #These are the unit vectors if 
        c = self.layer_gap         #if we create a lattice of flakes
        return [a, b, c]

    #This one is ugly. But it works
    def cell_coords(self):
        CC = self.CC
        CH = self.CH
        cos_CC = cos(pi / 6.0) * CC
        sin_CC = 0.5 * CC

        coords = np.zeros((self.natoms, 3))

        global atom
        atom = 0

        def add_coord(coord):
            global atom
            coords[atom] = coord
            atom += 1

        # first row of carbons
        add_coord([0, 0, 0])
        for col in range(self.x):
            add_coord([cos_CC + col * (2 * cos_CC), -sin_CC, 0])
            add_coord([2 * cos_CC + col * (2 * cos_CC), 0, 0])

        # other rows of carbons
        for row in range(self.y):
            if not row % 2:  # even (0 is even)
                y_offset = CC + (row * 3 / 2) * CC
                add_coord([0, y_offset, 0])
                for col in range(self.x):
                    add_coord([cos_CC + col * (2 * cos_CC), y_offset + sin_CC, 0])
                    add_coord([2 * cos_CC + col * (2 * cos_CC), y_offset, 0])
            else:  # odd
                y_offset = (row + 1) * CC * 3 / 2
                add_coord([0, y_offset, 0])
                for col in range(self.x):
                    add_coord([cos_CC + col * (2 * cos_CC), y_offset - sin_CC, 0])
                    add_coord([2 * cos_CC + col * (2 * cos_CC), y_offset, 0])

        # zig-zag Hs
        for col in range(self.x):
            add_coord([cos_CC + 2 * col * cos_CC, -sin_CC - CH, 0])
        for col in range(self.x):
            add_coord(
                [cos_CC + 2 * col * cos_CC, CC + (row * 3 / 2) * CC + sin_CC + CH, 0]
            )

        # armchair Hs
        cos_CH = cos(pi / 6.0) * CH
        sin_CH = 0.5 * CH
        for row in range(self.y):
            if row % 2:
                continue  # odd rows of hexagons have no Hs on ends
            y_offset = row * 3 / 2 * CC
            add_coord([-cos_CH, y_offset - sin_CH, 0])
            add_coord([-cos_CH, y_offset + CC + sin_CH, 0])
            add_coord([2 * self.x * cos_CC + cos_CH, y_offset - sin_CH, 0])
            add_coord([2 * self.x * cos_CC + cos_CH, y_offset + CC + sin_CH, 0])

        return coords
    #Why does this guy only go along the z-direction???
    def assign_molecules(self, lattice_dimensions):
        molecule_labels = []
        for z in range(lattice_dimensions[2]):
            labels = np.ones(self.natoms, dtype=int)
            molecule_labels.extend(list(labels + z))
        return molecule_labels

    #This one I understand, but does not invoke lattice dimensions and the
    #fact that we have copies.
    def assign_atom_labels(self, lattice_dimensions):
        atom_labels = np.zeros(self.natoms, dtype=int)
        atom_labels[: self.n_Cs] = 1  #This makes use of the fact that C
        atom_labels[self.n_Cs :] = 2  #were added first, and then the H were
        return list(atom_labels)      #added

    def assign_atom_states(self, lattice_dimensions):
        return list(np.zeros(self.natoms, dtype=int))
        

    #This should be done in the parametrization. For now, just set q=0.115
    #It seems it would be a lot easier to simply find all the edge H and use 
    #the bond graph to find the connecting CA or CT
    def assign_atom_charges(self, lattice_dimensions):
        q=0.115
        atom_charges = np.zeros(self.natoms)
        # zigzag carbons
        for col in range(self.x):
            # bottom row
            atom = 1 + 2 * col
            atom_charges[atom] = -q
            # top row
            atom = self.n_Cs - 2 - 2 * col
            atom_charges[atom] = -q

        # armchair carbons
        carbons_per_row = 1 + self.x * 2
        for row in range(self.y + 1):
            atom = carbons_per_row * row   #Left carbon
            atom_charges[atom] = -q
            atom = carbons_per_row * row + carbons_per_row - 1
            atom_charges[atom] = -q        #Right carbon

        # hydrogens
        atom_charges[self.n_Cs :] = q      #All the remaining are H.

        return list(atom_charges)

    #Where does lattice dimensions enter? Can I only make one flake?
    #Or is it a dummy as other molecules need it?
    #For a single flake I get it. Again, it's not too nice of a code but 
    #it does the job.
    def assign_bonds(self, lattice_dimensions):
        #This first section is just to figure out number of bonds and allocate
        #space for them in bonds=np.array
        firstrow = 3 + 4 * self.x
        a = 3 + 3 * self.x  # even rows
        b = 2 + 3 * self.x  # odd rows
        nbonds = firstrow + a * ((self.y - 1) / 2) + b * ((self.y + 1) / 2)
        bonds = np.empty((int(nbonds), 2), dtype=int) 
    
        global bond
        bond = 0

        def add_bond(a, b):
            global bond
            bonds[bond] = [a + 1, b + 1]
            bond += 1

        # along each row
        carbons_per_row = 1 + self.x * 2
        for row in range(self.y + 1):
            for atom in range(self.x * 2):
                atom1 = carbons_per_row * row + atom
                atom2 = carbons_per_row * row + atom + 1
                add_bond(atom1, atom2)

        # between rows of carbons
        for row in range(self.y):
            if row % 2:
                for col in range(self.x):
                    atom1 = carbons_per_row * row + col * 2 + 1
                    atom2 = carbons_per_row * (row + 1) + col * 2 + 1
                    add_bond(atom1, atom2)
            else:
                for col in range(self.x + 1):
                    atom1 = carbons_per_row * row + col * 2
                    atom2 = carbons_per_row * (row + 1) + col * 2
                    add_bond(atom1, atom2)

        # hydrogen bonds
        # bottom row
        for col in range(self.x):
            atom_c = col * 2 + 1
            atom_h = self.n_Cs + col
            add_bond(atom_c, atom_h)
        # top row
        for col in range(self.x):
            atom_c = col * 2 + 1 + self.y * carbons_per_row
            atom_h = self.n_Cs + col + self.x
            add_bond(atom_c, atom_h)
        # armchair hydrogens
        h_offset = self.n_Cs + 2 * self.x
        for row in range(self.y + 1):
            if row % 2:
                continue  # odd rows of hexagons have no Hs on ends
            atom_c = carbons_per_row * row
            atom_h = h_offset + row * 2
            add_bond(atom_c, atom_h)

            atom_c = carbons_per_row * (row + 1)
            atom_h = h_offset + 1 + row * 2
            add_bond(atom_c, atom_h)

            atom_c = carbons_per_row * (row + 1) - 1
            atom_h = h_offset + 2 + row * 2
            add_bond(atom_c, atom_h)

            atom_c = carbons_per_row * (row + 2) - 1
            atom_h = h_offset + 3 + row * 2
            add_bond(atom_c, atom_h)
        return bonds

    #This routine is never called. It is automatically fixed in connector.
    #This i because the oxidiser rearranges stuff, and we need to do everything
    #all over again. 
    def connection_types(self):
        sys.exit("who called me?")
        bond_types = [[1, 1], [1, 2]]
        angle_types = [[1, 1, 1], [1, 1, 2]]
        dihedral_types = [[1, 1, 1, 1], [1, 1, 1, 2], [2, 1, 1, 2]]
        improper_types = [[1, 1, 1, 1], [1, 1, 1, 2]]
        return bond_types, angle_types, dihedral_types, improper_types
