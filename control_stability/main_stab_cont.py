import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('..')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

from horizontal_tailplane import hor_run
from vertical_tailplane import ver_run

hor_run(aircraft)
ver_run(aircraft)

print("\n-----------------Summary-------------------\n")
print(f"CG_range is from {aircraft.X_cg_fwd} to {aircraft.X_cg_aft}\n")
print(f"Sh_S is {aircraft.Sh_S}")
print(f"Sv_V is {aircraft.Sv_S}")
print(f"\nC_m_alpha is {aircraft.C_m_alpha} [1/rad]\n")

