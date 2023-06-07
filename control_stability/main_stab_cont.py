import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('..')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

from Class_II_cg_estimation import cg_calc
from horizontal_tailplane import hor_run
from vertical_tailplane import ver_run
from control_surfaces import main_control_surface

cg_calc(aircraft)
hor_run(aircraft)
ver_run(aircraft)
exec(open("landing_gear.py").read(), {'aircraft':aircraft})
main_control_surface(aircraft)

print("\n-----------------Summary-Stabily-&-Control----------------\n")
print(f"CG_range is from {round(aircraft.X_cg_fwd, 3)} to {round(aircraft.X_cg_aft, 3)}\n")
print(f"Sh_S is {round(aircraft.AE_Sh_S, 3)}")
print(f"Sv_V is {round(aircraft.AE_Sv_S, 3)}")
print(f"\nC_m_alpha is {round(aircraft.C_m_alpha, 3)} [1/rad]\n")


