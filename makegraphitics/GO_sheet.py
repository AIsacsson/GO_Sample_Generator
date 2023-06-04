"""
Version 2022-04-01.
Works with makegraphitics_2022_04_01.zip

Main revision from version 2022-03-31: Removed the need to specify the force
field when creating the pristine sheet. This is now replaced by specifying
constants

CC_BOND_LENGTH = 1.4148 Å
LAYER_GAP = 3.4827 Å

The corresponding change have been made in the file graphene_cell.py



Creates a graphene oxide sheet with dimensions Lx, Ly
specified C:O ratio (co) and OH/epoxy ratio (ae):

There are some parts from a previous code that was general purpose and
allowed non periodic BC. Right now it only does PBC. The flakes are soon
coming back.

"""
import numpy as np
import sys
from math import pi, cos
from graphene_cell import Graphene
from oxidiser import Oxidiser
from crystal import Crystal
from params import Parameterise
from write_coords import Writer
from rectangle_graphene import Rectangle_Graphene

# ---------------------------------------------------------------------------
#
# Output file path management.
#
# mountpoint:    This I just use as I run it quite often from
#                the portable drive which mounts in different
#                locations depending which computer it is
#                currently plugged in to.
#
# destination_path: The output file eventually ends up.
#
# ---------------------------------------------------------------------------
mountpoint = "G:/"
destination_path = mountpoint + "Samples/"

# ---------------------------------------------------------------------------
#
# Geometry settings.
#
# ---------------------------------------------------------------------------
PBC = True               # Current version only runs with PBC

Nx = 10                # Nx, Ny are really Lx and Ly in Å
Ny = 10

# ---------------------------------------------------------------------------
#
# Oxidiser settings.
#
# The node_target variable (node_targ) is the biggest change here.
#
# node_targ = -1        The oxidizer runs as in the original version of
#                       makegraphitics
#
# node_targ = 0         No oxidation is made. The result is pristine graphene
#
#
# node_targ = integer   The number of nodes (islands) that the oxidizer
#                       targets. Once it has generated this number of
#                       islands, it stops making new ones. Note that one
#                       may still get fewer clusters, as islands can merge
#                       as they oxidize.
#
#                       The oxidizer will abort the oxidation if it converges
#                       too slowly towards the desired number of nodes.
#                       It will also abort if it has not reached the target
#                       because it ended up stuck.
#
#                       The oxidizer will return a flag that it
#                       aborted. If this flag is True, then the oxidizer
#                       aborted.
#
#                       Note that if the C:O ratio is small the targeted number
#                       of nodes may be impossible to reach.
#
# meth = "rf"           Oxidizer method. The random forrest "rf" is slow.
# or                    However, for few clusters, it gives more "natural"
# meth = "empirical"    shapes than the "empirical".
#                       The empirical method is the better choice for nearly
#                       homogenous samples.
#
#
# nif                   New island formation frequency. A low number will
#                       always create few nodes. Larger node targets will
#                       then never be reached. Too high values may cause
#                       an error in numpy, throwing an exception that we can't
#                       catch (see comments in oxidiser.py).
#
#                       How high values one can have depends on ae.
#                       if ae=1.0, then typically values of the order 1e21
#                       are needed. However, such large values may throw
#                       exceptions for lower values of ae. Thus, this
#                       must be tuned a bit carefully.
#
# force_OH_neighbor     Boolean. If set to True (default) nearest neighbor
#                       -OH groups cannot be on the same side of the basal
#                       plane. If set to False, neighboring -OH groups can
#                       appear on the same side.
#
# ---------------------------------------------------------------------------
co = 5.60                # C:O ratio,   between 1 and infinity
ae = 0.50                # -OH : -O- ratio

force_OH_side = True     # nn. -OH are on opposite sides of basal plane

node_targ = 2            # The number of nodes that the oxidizer should make.

Pristine = False         # If Pristine = False, we will oxidize.
# Pristine = True

meth = "rf"              # Oxidizer method: Random forrest regressor
# meth = "empirical"     # Oxidizer method: No clue really.

Box_size_z = [-40, 40]   # Box size, zlo and zhi in Å. Sheet will be at z=0.

nif = 1.0e17             # New island formation frequency

if Pristine:
    node_targ = 0

# ----------------------------------------------------------------------------
#
# What comes below is a bit silly (it's from the original programming)
# Essentially it just extracts the bond lengths between the carbons in
# graphene. The line config[focefield]["CC"] can just be replaced by the
# number 1.4293 which is what it reads from the OPLS parameters.
# Something is a bit weird here though, if one opens the noyaml.py.
# The C-C distance is different between the different forcefields.
# But once the structure is relaxed in LAMMPS, the distance should become
# 1.4148 (which is the value usually quoted, rather than 1.43 as above).
# This doesn't matter much though, if we get a bond-length that is correct
# after relaxation.
#
# layout = [x_cells, y_cells, 1] makes an array of unit cells with this
# dimension. Note that the last number must be 1 at the moment.
#
# ----------------------------------------------------------------------------
CC_BOND_LENGTH = 1.4148
LAYER_GAP = 3.4827

if PBC:
    unit_cell_x = 2.0 * CC_BOND_LENGTH * cos(pi / 6.0)
    unit_cell_y = 3.0 * CC_BOND_LENGTH
    x_cells = int(Nx / unit_cell_x)
    y_cells = int(Ny / unit_cell_y)
    layout = [x_cells, y_cells, 1]


# ----------------------------------------------------------------------------
#
# The two lines below genereates the PBC-graphene flake. However, if one wants
# to generate a pristine sheet, the coding is a bit stupid as the vdW
# parameters are not included in the flake, but must be set specifically.
# This is temporarily fixed by setting node_target=0 and calling the
# oxidiser anyway.
#
# The force-field parameters used to be a parameter here. The only things
# that were used were CC_BOND_LENGTH and LAYER_GAP
# which should not really depend on the force field.
# They can now be entered as parameters by calling the routine as
#
# motif = mg.molecules.Graphene(cc_bond_length=CC_BOND_LENGTH,
#                               layer_gap=LAYER_GAP)
#
# If mg.molecules are called without these parameters, they default to
# CC_BOND_LENGTH = 1.4148 and LAYER_GAP = 3.4827
#
# ----------------------------------------------------------------------------
# motif = mg.molecules.Graphene(cc_bond_length=CC_BOND_LENGTH,
#                               layer_gap=LAYER_GAP)
if PBC:
    motif = Graphene()
    flake = Crystal(motif, layout)
else:
    motif = Rectangle_Graphene(Nx, Ny)
    flake = Crystal(motif, [1, 1, 1])


# ----------------------------------------------------------------------------
#
# Oxidiser init -  Here there are some default parameters hidden. I disabled
#                  a few as I wen't along and changed things. Above all,
#                  it cannot write xyz-files, it doesn't do edge oxidation at
#                  the moment, and I need to fix some flags....
#
# ----------------------------------------------------------------------------
oxidiser = Oxidiser(ratio=co, surface_OHratio=ae,
                    new_island_freq=nif, method=meth,
                    node_target=node_targ,
                    force_OH_side=True)

# ----------------------------------------------------------------------------
#
# Oxidiser step   -  In the original version, only the flake was
# (reactor)          returned. I have added two more outputs from it
#
#                    nodes: Integer that tells how many nodes where genereated
#
#                    stuck: Boolean that returns true if the oxidizer
#                           got stuck and/or deemed it could not reach the
#                           number of desired nodes. If it returns False,
#                           it achieved the node_target. It's an ugly
#                           and experimental implementation that works at the
#                           moment, but it must be fixed.
#
# ----------------------------------------------------------------------------
flake, nodes, stuck = oxidiser.react(flake)

# ----------------------------------------------------------------------------
#
# PARAMETRIZATION - I've started on this bit. At the moment I still focus on
#                   PBC. The idea is to put all force-field info in a .json
#                   file, and that we then later can chose different force
#                   fields and parameters.
#
# ----------------------------------------------------------------------------
ffield = "OPLS_Jorgensen2009"
if ((not Pristine) and (not stuck)) or (Pristine):
    print("Parametrizing sample with forcefield " + ffield)
    Parameterise(flake, forcefield=ffield)

# ----------------------------------------------------------------------------
#
# SOME STUFF - This part is just stuff I used for my self when building
#              libriaries etc.
#
# Here I make use of the fact that I know which atom types are which.
# Although the oxidizer targets a certaing co and ae ratio, it does not always
# succeed. For a 200x200 Å2 flake it misses by a few percent in either
# direction when it comes to the ae. The co is usually more spot on.
# I use these when I make large sample batches and continously loop to make
# libraries to decide in which library to put them in.
#
# I also set up my file-name convetion here. Depending on versions I either
# put them in a target library if they met a tolerance criterion, or I put
# them in folders that matched the actual co and ae ratio. I go with the
# former one at the moment as I otherwise get silly co-ratios like 49.9 etc.
#
# ----------------------------------------------------------------------------
if (not Pristine) and (not stuck):
    N_OH = np.sum(np.array(flake.atom_labels) == 4)
    N_epoxy = np.sum(np.array(flake.atom_labels) == 6)
    N_O = N_OH + N_epoxy
    N_C = np.sum(np.array(flake.atom_labels) == 1)
    N_C = N_C + np.sum(np.array(flake.atom_labels) == 3)

    actual_ae = N_OH / N_O
    actual_co = N_C / N_O
    Onoderatio = N_O / nodes

    ae_ratio_str = "{:.2f}".format(actual_ae)
    co_ratio_str = "{:.2f}".format(actual_co)
    target_ae_str = "{:.2f}".format(ae)
    target_co_str = "{:.2f}".format(co)
    Onoderatio_str = "{:.2f}".format(Onoderatio)

    print("Actual ae: " + ae_ratio_str + "  (Target:" + target_ae_str + ")")
    print("Actual co: " + co_ratio_str + "  (Target:" + target_co_str + ")")
    print("Current nif:" + str(nif))
elif Pristine:
    print("No oxidation done. Pristine graphene/flake made")
else:
    print("Oxidizer did not meet node target. Stopping")

# ----------------------------------------------------------------------------
#
# LAMMPS data-file writer
#
# Here I added the possibility set the z-dimensions of the box, and
# to automatically add the OPLS TIP3P/ew parameters for water
# to the .data file so that we don't have to do this.
#
# For a future update, this should also include different forcefields and
# water models, each to go with each other.
#
# There is a parameter verbose=True/False. Setting it to False will supress
# a bit of output to the console.
#
# ----------------------------------------------------------------------------
if (not Pristine) and (not stuck):
    fname = "go_"+str(Nx)+"x"+str(Ny)+"A2_co"+target_co_str+"_ae"+target_ae_str
    fname = destination_path + fname
elif Pristine:
    fname = destination_path + "graphene_"+str(Nx)+"x"+str(Ny)+"A2"

if not stuck:
    print("Writing Lammps file....")
    flake.box_dimensions[2][:] = Box_size_z
    output = Writer(flake, fname, verbose=False, fix_the_water=True)
    output.write_lammps(fname + ".data")

# ----------------------------------------------------------------------------
#
# PDB - WRITER
#
# I haven't used the pdb-write since I played around with ATB server, LigParGen
# and the weird Russian server. Not sure it works anymore. I need to check this
# This one doesn't include water. It just says where the atoms (C,O,H) are, and
# where the bonds between them are, and that's basically it.
#
# ----------------------------------------------------------------------------
print("Writing pdb-file")
output.write_pdb(fname+".pdb")

# ----------------------------------------------------------------------------
#
# XYZ - Writer
#
# As I recall this one I broke when I included the water into the LAMMPS
# .data file. It's an easy fix to get it back, but I never really saw the
# need for it. As it doesn't work anymore, it's not possible to make movies
# of the oxidation process. But I don't care too much about that at the moment.
# It may come handy later on if we wish to make other oxidezer options other
# than "rf" and "empirical". Then it might be good if one can track the
# evolution for debugging purposes
#
# ----------------------------------------------------------------------------
# output.write_xyz(fname+".xyz")

print(flake.angle_types)
print(flake.angle_coeffs)
