
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('..')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

def longitudinal_position_landing_gear(aircraft):
    pass

def lateral_position_landing_gear(aircraft):
    pass

def height_landing_gear(aircraft):
    pass

def run(aircraft):
    longitudinal_position_landing_gear(aircraft)
    lateral_position_landing_gear(aircraft)
    height_landing_gear(aircraft)

# check most fwd position cg
# check most aft position cg
