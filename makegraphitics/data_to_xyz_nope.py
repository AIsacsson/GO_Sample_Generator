'''
This is a remnant of the original make-graphitics. It is not used.
'''

import yaml
from scripts.molecules import *
from scripts import *
import numpy as np
import sys

assert 1 == 0

sim = ReadLammpsData(sys.argv[1])


name = "out"
output = Writer(sim, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
