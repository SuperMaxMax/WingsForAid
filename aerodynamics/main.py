import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
from parameters import UAV
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate

import wing_planform as wp
import horizontal_tail_design as htd

aircraft = UAV('aircraft')
#
airfoil = aircraft.airfoil 

#print(wp.plot_lift_distr(aircraft))
wing_data = wp.plot_lift_distr(aircraft)
AR = wing_data[0]
taper = wing_data[1]
twist = wing_data[2]
span_eff = wing_data[3]
CL_wing = wing_data[4]
CD_induced_wing = wing_data[5]
incidence_wing = wing_data[6]

 


#AR, Lambda, alpha_twist, span_eff, C_L_wing, CD_induced, span_eff, i_w


