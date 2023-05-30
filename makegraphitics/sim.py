from connector import Connector
import numpy as np
import noyaml

class Sim(object):
    def generate_connections(self):
       #If we don't have atom labels or bonds between them it won't work
        try:
            self.atom_labels, self.bonds
        except AttributeError:
            raise Exception("Simulation has not been assigned atom types and bonds yet")

        connect = Connector()   #Does not have an __init__

        #Start with the bonds. 
        self.bond_types = connect.find_bond_types(self.atom_labels, self.bonds)
        
        self.bond_labels = connect.bond_labels(self.atom_labels, 
                                               self.bonds, 
                                               self.bond_types)


        #The bond graph is used in subsequent functions.
        self.bond_graph = self.generate_bond_graph(self.bonds)


        #Next, we sort out the angles.
        self.angles = connect.angles(self.bonds, self.bond_graph)
        
        self.angle_types = connect.find_angle_types(self.atom_labels, 
                                                    self.angles)
        
        self.angle_labels = connect.angle_labels(self.atom_labels, 
                                                 self.angles, 
                                                 self.angle_types)

        #Dihedrals,..... 
        self.dihedrals = connect.dihedrals(self.bonds, self.bond_graph)
        
        self.dihedral_types = connect.find_dihedral_types(self.atom_labels, 
                                                          self.dihedrals)
        
        self.dihedral_labels = connect.dihedral_labels(self.atom_labels, 
                                                       self.dihedrals, 
                                                       self.dihedral_types)

        #and finally the impropers. 
        self.impropers = connect.impropers(self.bonds, self.bond_graph)

        self.improper_types = connect.find_improper_types(self.atom_labels, 
                                                          self.impropers)

        self.improper_labels = connect.improper_labels(self.atom_labels, 
                                                       self.impropers, 
                                                       self.improper_types)


    #The bond graph is a dictionary which tells which bonds connect to an 
    #atom. In contrast to other lists, there is no +1 to atoms or bonds.
    #Hence, in the bond-graph, an entry N:{n1,n2,n3} refers to atom N+1
    #and bonds n1+1, n2+1 and n3+1, etc.
    def generate_bond_graph(self, bonds):
        N = int(np.amax(bonds))  # Number of atoms (clever girl....)
        bond_graph = dict()
        for i in range(N):          #Initialize for each atom- 
            bond_graph[i] = set()   #looks clumsy. Python must have one-liner
                                    #for this.
                                    
        #I need to read up on dictionaries and sets. I don't immediately see how this 
        #loop works. But if I try to add an item to a set that is already
        #in it, it will not create a duplicate.
        for bond in bonds:        
            bond_graph[bond[0] - 1].add(bond[1] - 1)
            bond_graph[bond[1] - 1].add(bond[0] - 1)
        return bond_graph

    def bonded_to(self, centre):
        return list(self.bond_graph[centre])

    def validate(self):
        n_atoms = len(self.coords)
        for attr in ["molecule_labels", "atom_charges", "atom_labels"]:
            if hasattr(self, attr):
                assert len(getattr(self, attr)) == n_atoms, attr
        if abs(np.sum(self.atom_charges)) > 0.01:
            print("WARNING: charges do not sum to zero",
                  np.sum(self.atom_charges))

    def crystal_params(self):
        return noyaml.get_params()
