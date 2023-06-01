import sys

sys.path.append("..")

# Start your import below this
from parameters import UAV, atmosphere
from box_drop import BOX, subdivide
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mlpatch
import matplotlib.colorbar as cbar
import time
import itertools as it

# INITIALISE

# get defaults design values
AC = UAV('aircraft')
print(AC.__dict__)
print("\n")

atm = atmosphere()
print(atm.__dict__)
print("\n")

box = BOX()
print(box.__dict__)
print("\n")

# Simulation param
dt_sim = 1 / 10E3  # [s] between frames
IT_max = 10E4
plot_traj = False
plot_bounds = True

# Maneuver
# m limits
M_V_app = [AC.V_s_min, AC.OP_V_boxlim] # [m/s]
M_n_app = [-0.5, 4.5] # [g]
M_pitch = [-45, 10] # [deg]
M_heading = [-180,180] # [deg]
M_Hmin = [AC.OP_hmin, 50] # [m]
# m combinations
N_div_M = 2
M_VAR = {
    "V_app": subdivide(M_V_app, N_div_M),
    "n_app": subdivide(M_n_app, N_div_M),
    "pitch": subdivide(M_pitch, N_div_M),
    "heading": subdivide(M_heading, N_div_M),
    "Hmin": subdivide(M_Hmin, N_div_M)
}
allNames = sorted(M_VAR)
M_COMB = list(it.product(*(M_VAR[Name] for Name in allNames)))

# Drop Cases
# limits
C_Mbox = [10, 23] # [kg]
C_Vw = [0, AC.OP_V_crosswind] # [kg]
C_w_heading = [0, +180] # [deg]
# c combinations
N_div_C = 2
C_VAR = {
    "Mbox": subdivide(D_Mbox, N_div_C),
    "Vw": subdivide(C_Vw, N_div_C),
    "w_heading": subdivide(C_w_heading, N_div_C),
}
allNames = sorted(C_VAR)
C_COMB = list(it.product(*(C_VAR[Name] for Name in allNames)))

# Uncertainties
# limits
N_div_U = 1
U_Mbox = 1 # [kg]
U_T_drop = 1 # [s]
U_T_brake = 1 # [s]
U_T_flap = 1 # [s]
U_pos = 2 # [m]
U_Vw = 5 # [m/s]
U_w_heading = 10 # [deg]
# U combinations
N_div_U = 1
U_VAR = {
    "Mbox": subdivide(U_Mbox, N_div_U),
    "T_drop": subdivide(U_T_drop, N_div_U),
    "T_brake": subdivide(U_T_brake, N_div_U),
    "T_flap": subdivide(U_T_flap, N_div_U),
    "pos": subdivide(U_pos, N_div_U),
    "Vw": subdivide(U_Vw, N_div_U),
    "w_heading": subdivide(U_w_heading, N_div_U)
}
allNames = sorted(U_VAR)
U_COMB = list(it.product(*(U_VAR[Name] for Name in allNames)))

# START TRYING COMBINATIONS
N_Case = N_Man = N_Unc = 0 # set counter
print("\n================= ================= =================",
      len(C_COMB), "Drop Cases Combinations")
for DropCase in C_COMB:
    print(N_Case, "/", len(C_COMB),"drop case", str(DropCase))
    CaseResults = {
        "best M": [],
        "avg Racc": [],
        "worst Racc": [],
        "avg Dfit": [],
        "worst Dfit": []
    }
    print("\n================= =================",
          len(M_COMB), "Maneuver Combinations")
    for ManCase in M_COMB:
        print(N_Man, "/", len(M_COMB), "maneuver", str(ManCase))
        ManResults = {
            "Racc": [],
            "Dfit": [],
            "worst V_impact": [],
            "worst a_max": []
        }
        print("\n=================",
              len(U_COMB), "Disturbance Uncertainty Combinations")
        for Drop in U_COMB:
            print(N_Unc, "/", len(U_COMB), "variation", str(Drop))
            DropResults = {
                "impact DX [m]": [],
                "impact DY [m]": [],
                "impact velocity [m/s]": [],
                "max acceleration [g]": [],
                "max velocity [m/s]": [],
                "time to land [s]": []
            }
            # define inputs
            # run simu
            N_Unc += 1
        N_Man += 1
    N_Case += 1

# Mparam = f(Case param)
# CaseResults = f(Case param)