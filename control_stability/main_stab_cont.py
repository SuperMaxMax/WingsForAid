import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('.')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

from Class_II_cg_estimation import cg_calc
from horizontal_tailplane import hor_run
from vertical_tailplane import ver_run
#from landing_gear import run

cg_calc(aircraft)
hor_run(aircraft)
#ver_run(aircraft)
# run()

print("\n-----------------Summary-------------------\n")
print(f"CG_range is from {aircraft.X_cg_fwd} to {aircraft.X_cg_aft}\n")
print(f"Sh_S is {aircraft.CS_Sh_S}")
print(f"Sv_V is {aircraft.CS_Sv_S}")
print(f"\nC_m_alpha is {aircraft.C_m_alpha} [1/rad]\n")

