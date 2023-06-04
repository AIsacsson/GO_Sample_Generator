# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 20:09:34 2022

@author: aisac
"""

# ----------------------------------------------------------------------------
#
# This is actually misleading. The numbers below are the UA numbers to be
# used to find the vdW parameters when reading from the Jorgensen .prm
# file. So a more appropriate name would be get UA_defs
#
# ----------------------------------------------------------------------------
def get_vdw_defs(force_field):
    vdw_defs = {"OPLS_Jorgensen2009":
                {
                    1: 90,  # Cg, graphitic (aromatic)
                    2: 91,  # Hg, graphitic edge
                    3: 101,  # Ct, tertiary C-OH
                    4: 96,  # Oa, C-OH
                    5: 97,  # Ha, C-OH
                    6: 122,  # Oe, epoxy
                    7: 109,  # Oa, C-OH
                    8: 209,  # Cc, Carboxylic carbon
                    9: 210,  # Oc, Ketone oxygen
                    10: 211,  # Oa, alcohol
                    11: 108,  # Cb, Benzyl
                    }
                }
    return vdw_defs[force_field]


# ----------------------------------------------------------------------------
#
# Partial charges for OPLS-AA (Jorgensen).
#
#
# atom_label    atom_states charge   Description
#     1             0          0         CA
#     2             0        0.115      Edge hydrogen
#     1             3       -0.115      Edge carbon connected to hydrogen
#     3             1        0.265      CT connected to -HO, above plane
#     3            -1        0.265      CT connected to -OH, below plane
#     3             2        0.200      CT connected to -O-, above plane
#     3            -2        0.200      CT connected to -O-, below plane
#     4             0       -0.683      OH connected to CT
#     5             0        0.418      HO connected to CT
#     6             0       -0.400      -O- connected to CTs
#     7             0
#     8             0
#     9             0
#    10             0
#    11             0                   TIP3P O
#    12             0                   TIP3P H
# ----------------------------------------------------------------------------
def get_charges(force_field):
    charges = {"OPLS_Jorgensen2009":
               {
                   1: {0:  0.000, 3: -0.115},
                   2: {0:  0.115},
                   3: {1:  0.265, 2:  0.200},
                   4: {0: -0.683},
                   5: {0:  0.418},
                   6: {0: -0.400}
                   }
               }

    return charges[force_field]


def get_masses(force_field):
    masses = {"OPLS_Jorgensen2009":
              {
                   1:  {1: 12.011, 2: "# C, sp2 basal plane"},
                   2:  {1:  1.008, 2: "# H, edge termination"},
                   3:  {1: 12.011, 2: "# C, sp3 basal plane"},
                   4:  {1: 15.999, 2: "# O, basal plane -OH"},
                   5:  {1:  1.008, 2: "# H, basal plane -OH"},
                   6:  {1: 15.999, 2: "# O, basal plane epoxy"},
                   7:  {1: 15.999, 2: "# Edge Oxygen"},
                   8:  {1: 12.011, 2: "# Edge Carbon"},
                   9:  {1: 15.999, 2: "# Edge Oxygen"},
                   10: {1: 15.999, 2: "# Edge Oxygen"},
                   11: {1: 12.011, 2: "# Edge carbon"},
                   12: {1: 15.999, 2: "# O, H2O"},
                   13: {1:  1.008, 2: "# H, H20"}
                   }
              }

    return masses[force_field]


# ----------------------------------------------------------------------------
#
# Pair coeffs for OPLS-AA (Jorgensen2009).
# Mixing geometric
#
#
# atom_label     UA        Pair coeff.          Description
#     1          90       0.070   3.55          CA    "Aromatic C"
#     2          91       0.030   2.42          HA    "Aromatic H-C"
#     3         101       0.066   3.50          CT    "Alcohol R3COH"
#     4          96       0.170   3.12          OH    "Alcohol -OH"
#     5          97       0.000   0.00          HO    "Alcohol -OH"
#     6         122       0.140   2.90          OS    "Dialkyl Ether -O-"
#     7         109       0.170   3.07          OH    "Phenol -OH"
#     8         209       0.105   3.75          C     "Carboxylic Acid -COOH"
#     9         210       0.210   2.96          O     "Carboxylic Acid C=O"
#    10         211       0.170   3.00          OH    "Carboxylic Acid -OH"
#    11         108       0.070   3.55          CA    "Phenol C-OH"
#
# ----------------------------------------------------------------------------
def get_pair_coeffs(force_field):
    pair_coeffs = {"OPLS_Jorgensen2009":
                   {
                       1:  {1: 0.070, 2: 3.550, 3: "# CA Aromatic C,           OPLS-UA  90 (Jorgensen2009)"},
                       2:  {1: 0.030, 2: 2.420, 3: "# HA Aromatic H-C          OPLS-UA  91 (Jorgensen2009)"},
                       3:  {1: 0.066, 2: 3.500, 3: "# CT Alcohol R3COH         OPLS-UA 101 (Jorgensen2009)"},
                       4:  {1: 0.170, 2: 3.120, 3: "# OH Alcohol -OH           OPLS-UA  96 (Jorgensen2009)"},
                       5:  {1: 0.000, 2: 0.000, 3: "# HO Alcohol -OH           OPLS-UA  97 (Jorgensen2009)"},
                       6:  {1: 0.140, 2: 2.900, 3: "# OS Dialkyl Ether -O-     OPLS-UA 122 (Jorgensen2009)"},
                       7:  {1: 0.170, 2: 3.070, 3: "# OH Phenol -OH            OPLS-UA 109 (Jorgensen2009)"},
                       8:  {1: 0.105, 2: 3.750, 3: "# C  Carboxylic Acid -COOH OPLS-UA 209 (Jorgensen2009)"},
                       9:  {1: 0.210, 2: 2.960, 3: "# O  Carboxylic Acid C=O   OPLS-UA 210 (Jorgensen2009)"},
                       10: {1: 0.170, 2: 3.000, 3: "# OH Carboxylic Acid -OH   OPLS-UA 211 (Jorgensen2009)"},
                       11: {1: 0.070, 2: 3.550, 3: "# CA Phenol C-OH           OPLS-UA 108 (Jorgensen2009)"},
                       12: {1: 0.102, 2: 3.188, 3: "# O TIP3P/x "},
                       13: {1: 0.000, 2: 0.000, 3: "# H TIP3P/x "}
                       }
                   }
    return pair_coeffs[force_field]


# ----------------------------------------------------------------------------
#
# Bond  coeffs for OPLS-AA (Jorgensen2009).
# Mixing geometric
#
#
#   bond_type   AA        Bond coeff.          Description
#     [1 1]    48 48     469.00 1.4000          CA - CA
#     [1 2]    48 49     367.00 1.0800          CA - HA
#     [1 3]    48 13     317.00 1.5100          CA - CT
#     [3 1]    13  ?                            CT - HX     
#     [3 3]    13 13     268.00 1.5290          CT - CT
#     [3 4]    13  5     320.00 1.4100          CT - OH
#     [4 5]     5  7     553.00 0.9450          OH - HO
#     [3 6]    13 20     320.00 1.4100          CT - OS
#     [x y]           Edges
#     [x y]           Edges
#     [x y]           Edges
#     [x y]           Edges
#     [x y]           Edges
#
# ----------------------------------------------------------------------------
def get_bond_coeffs(ffield, bond_type_list):

    bond_pot = {"OPLS_Jorgensen2009":
                {
                    1: {1: [469.00, 1.4000, "# CA-CA  OPLS-AA [48 48] (Jorgensen2009)"],
                        2: [367.00, 1.08,   "# CA-CH  OPLS-AA [48 49] (Jorgensen2009)"],
                        3: [317.00, 1.5100, "# CA-CT  OPLS-AA [13 48] (Jorgensen2009)"]},
                    3: {3: [268.00, 1.5290, "# CT-CT  OPLS-AA [13 13] (Jorgensen2009)"],
                        4: [320.00, 1.4100, "# CT-OH  OPLS-AA [ 5 13] (Jorgensen2009)"],
                        6: [320.00, 1.4100, "# CT-OS  OPLS-AA [13 20] (Jorgensen2009)"]},
                    4: {5: [553.00, 0.9450, "# OH-HO  OPLS-AA [ 5  7] (Jorgensen2009)"]}
                    }
                }

    btl = bond_type_list     # make a copy in order not to change anything
    bond_coeffs = dict()

    for i in range(len(btl)):
        btl[i].sort()
        c1 = bond_pot[ffield][btl[i][0]][btl[i][1]][0]
        c2 = bond_pot[ffield][btl[i][0]][btl[i][1]][1]
        tag = bond_pot[ffield][btl[i][0]][btl[i][1]][2]
        bond_coeffs[i+1] = {1: c1, 2: c2, 3: tag}

    return bond_coeffs


def get_atom_labels():
    atom_label = {
        1: "C",  # Cg, graphitic (aromatic)
        2: "H",  # Hg, graphitic edge
        3: "C",  # Ct, tertiary C-OH
        4: "O",  # Oa, C-OH
        5: "H",  # Ha, C-OH
        6: "O",  # Oe, epoxy
        7: "O",  # Oa, C-OH
        8: "C",  # Cc, Carboxylic carbon
        9: "O",  # Oc, Ketone oxygen
        10: "O",  # Oa, alcohol
        11: "C",  # Cb, Benzyl
    }  # OPLS definitions
    return atom_label

"""
  # all OPLS atom types that are introduced by oxidation
  sim.vdw_defs = {
      1: 90,  # Cg, graphitic (aromatic)
      2: 91,  # Hg, graphitic edge
      3: 101,  # Ct, tertiary C-OH 
      4: 96,  # Oa, C-OH
      5: 97,  # Ha, C-OH
      6: 122,  # Oe, epoxy
      7: 109,  # Oa, C-OH
      8: 209,  # Cc, Carboxylic carbon
      9: 210,  # Oc, Ketone oxygen
      10: 211,  # Oa, alcohol
      11: 108,  # Cb, Benzyl
      12: 213, # C, carboxylate -COO
      13: 214, # O, carboxylate -COO
      14: 349, # Na+
      15: 354, # Ca 2+ 
  }  # OPLS definitions
  """