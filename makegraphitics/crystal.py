# ---------------------------------------------------------------------------
#
# Version 2022-04-01:      Removing all references to the force field...
#                          Removing all stuff that is not needed for single
#                          sheet PBC
#
# ---------------------------------------------------------------------------
from lattice import Lattice
from sim import Sim


class Crystal(Sim):
    def __init__(self, molecule, lattice_dimensions):

        self.molecule = molecule
        self.lattice = self.init_lattice()
        self.lattice_dimensions = lattice_dimensions
        self.generate_structure()
        self.generate_bonds()
        self.generate_connections()

    def init_lattice(self):
        self.cell_coords = self.molecule.cell_coords()
        return Lattice(self.molecule.cell_shape())

    def generate_structure(self):
        self.lattice_points = self.lattice.create_lattice_points(
            self.lattice_dimensions)

        self.coords = self.lattice.cell_onto_lattice(self.cell_coords,
                                                     self.lattice_points)

        self.box_dimensions = self.lattice.system_size(self.lattice_dimensions)

        self.molecule_labels = self.molecule.assign_molecules(
            self.lattice_dimensions)

        self.atom_labels = self.molecule.assign_atom_labels(
            self.lattice_dimensions)

        self.atom_charges = self.molecule.assign_atom_charges(
            self.lattice_dimensions)

        self.atom_states = self.molecule.assign_atom_states(
            self.lattice_dimensions)  # This is just an array of zeros.

    def generate_bonds(self):
        self.bonds = self.molecule.assign_bonds(self.lattice_dimensions)


# def determine_lattice(self, lattice_dimensions):
#     if not lattice_dimensions:
#         return [1, 1, 1]
#     elif lattice_dimensions == "vdw":
#         return self.lattice.lattice_size_vdw(
#             self.config["system"]["vdw_cutoff"])
#     elif lattice_dimensions == "layers":
#         return self.lattice.lattice_size_layers(
#             self.config["system"]["vdw_cutoff"],
#             self.config["system"]["N_layers"])
#     else:
#         return lattice_dimensions
