import sys

sys.path.append("..")

# Start your import below this
from parameters import UAV, atmosphere
from Drop import BOX
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mlpatch
import matplotlib.colorbar as cbar
import time
import itertools as it
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import random
#from sympy import Point, Polygon
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, MultiPoint
import shapely

# FUNCTIONS
def SimuDrop(DropCon, DropMan, DropUnc, DropResults):
    # get states and apply uncertainties
    CR_Vw = DropCon[1] + DropUnc[1]
    CR_wheading = DropCon[2] + DropUnc[2]

    Approach = {
        "V_app": DropMan[0],
        "n_app": DropMan[1],
        "pitch": DropMan[2],
        "heading": DropMan[3],
        "Hmin": DropMan[4],
    }
    # calc AC traj s.t. h_min == 15m
    Dh = 0
    if Approach["pitch"] != 0:
        R = Approach["V_app"]  / Approach["n_app"]   # approach pull radius
        Dh = R * np.abs(np.sin(Approach["pitch"] ))
    AC_pos = [0, 0, Approach["Hmin"] + Dh]

    # disturb drop time
    time_dev_drop = DropUnc[3]
    if time_dev_drop != 0:
        if Approach["pitch"] != 0:
            Approach["pitch"] += (atm.g / Approach["n_app"]) * (Approach["n_app"] - 1) * time_dev_drop
        AC_pos = np.array([0, 0, Approach["Hmin"] + Dh])\
                 + np.array([
                    np.cos(Approach["pitch"]),
                    np.sin(Approach["pitch"]),
                    0
                ]) * Approach["V_app"] * time_dev_drop

    # apply uncertainties to box
    newbox = BOX(DropCon[-1])
    newbox.mass += DropCon[0] + DropUnc[0]
    newbox.DT_brake += DropUnc[4]
    newbox.DT_flaps += DropUnc[5]

    # get box position relative to AC
    boxDX = 0.2+0.4
    boxDY = 0.1+0.4/2
    boxDZ = 0
    box_pos = {
        0: [0 * boxDX, boxDY, boxDZ],
        1: [0 * boxDX, -boxDY, boxDZ],
        2: [5 * boxDX, boxDY, boxDZ],
        3: [5 * boxDX, -boxDY, boxDZ],
        4: [1 * boxDX, boxDY, boxDZ],
        5: [1 * boxDX, -boxDY, boxDZ],
        6: [4 * boxDX, boxDY, boxDZ],
        7: [4 * boxDX, -boxDY, boxDZ],
        8: [2 * boxDX, boxDY, boxDZ],
        9: [2 * boxDX, -boxDY, boxDZ],
        10: [3 * boxDX, boxDY, boxDZ],
        11: [3 * boxDX, -boxDY, boxDZ],
    }
    box_pos = np.array(box_pos[DropCon[3]])

    # get initial state
    rho = atm.rho0
    g0 = atm.g
    N = t = 0

    Pos = AC_pos + box_pos + np.array([DropUnc[6], DropUnc[7], 0])
    Vel = np.array([
        Approach["V_app"]*np.cos(np.pi*Approach["pitch"]/180),
        0,
        Approach["V_app"]*np.sin(np.pi*Approach["pitch"]/180)
    ])
    Acc = [0, 0, -g0]

    State = [[[Pos[0]],[Pos[1]],[Pos[2]]],
             [[Vel[0]],[Vel[1]],[Vel[2]],[magnitude(Vel)]],
             [[Acc[0]],[Acc[1]],[Acc[2]],[magnitude(Acc)]],
             [t]]

    # loop until ground is hit
    while Pos[-1] > 0 and N < IT_max:
        newbox.compute_CD(t)
        # get relative airspeed
        AirVel = np.array([
            Vel[0] - CR_Vw * np.cos(np.pi * CR_wheading / 180),
            Vel[1] - CR_Vw * np.sin(np.pi * CR_wheading / 180),
            Vel[2]
        ])
        # compute
        Acc = np.array([
            -np.sign(AirVel[0]) * (1/2) * newbox.CD * rho * newbox.S * (AirVel[0] ** 2) / newbox.mass,
            -np.sign(AirVel[1]) * (1/2) * newbox.CD * rho * newbox.S * (AirVel[1] ** 2) / newbox.mass,
            -np.sign(AirVel[2]) * (1/2) * newbox.CD * rho * newbox.S * (AirVel[2] ** 2) / newbox.mass - g0
        ])
        Vel = np.add(Vel, Acc*dt_sim)
        Pos = np.add(Pos, Vel*dt_sim)

        # append data
        State[0][0].append(Pos[0])
        State[0][1].append(Pos[1])
        State[0][2].append(Pos[2])
        State[1][0].append(Vel[0])
        State[1][1].append(Vel[1])
        State[1][2].append(Vel[2])
        State[1][3].append(magnitude(Vel))
        State[2][0].append(Acc[0])
        State[2][1].append(Acc[1])
        State[2][2].append(Acc[2])
        State[2][3].append(magnitude(Acc))
        State[3].append(t)

        t += dt_sim
        N += 1
    if N + 1 == IT_max:
        print("SIMULATION ABORTED")

    # store results
    DropResults["impact DX [m]"].append(State[0][0][-1])
    DropResults["impact DY [m]"].append(State[0][1][-1])
    DropResults["impact velocity [m/s]"].append(State[1][2][-1])
    DropResults["max acceleration [g]"].append(np.max(State[2][-1]) / g0)
    DropResults["max velocity [m/s]"].append(np.max(State[1][-1]))
    DropResults["time to land [s]"].append(State[-1][-1])

    return State

def magnitude(vector):
    mag = np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
    return mag

def subdivide(range,N):
    var = []
    if N <= 0:
        if range[0]==0: var.append(0)
        else: var.append(np.mean(range))
    elif N == 1:
        var = range
    else:
        i = 0
        DR = range[1] - range[0]
        dr = DR/N
        while i <= N:
            var.append(range[0] + i * dr)
            i += 1
    return var

def PlotTraj(states, title):

    plt.suptitle("Box drop")

    # trajectory
    plt.subplot(231)
    plt.plot(states[0][0], states[0][1])
    plt.title("trajectory XY")
    plt.subplot(232)
    plt.plot(states[0][0], states[0][2])
    plt.title("trajectory XZ")
    plt.subplot(233)
    plt.plot(states[1][1], states[1][2])
    plt.title("trajectory YZ")

    plt.subplot(234)
    plt.plot(states[3], states[0][0], label='Px')
    plt.plot(states[3], states[0][1], label='Px')
    plt.plot(states[3], states[0][2], label='Px')
    plt.title("position vs time")
    plt.legend(loc='lower center', ncol=1, bbox_to_anchor=(0.5, -0.7))

    plt.subplot(235)
    plt.plot(states[3], states[1][0], label='Vx')
    plt.plot(states[3], states[1][1], label='Vy')
    plt.plot(states[3], states[1][2], label='Vz')
    plt.title("velocity vs time")
    plt.legend(loc='lower center', ncol=1, bbox_to_anchor=(0.5, -0.7))

    plt.subplot(236)
    plt.plot(states[3], states[2][0], label='Ax')
    plt.plot(states[3], states[2][1], label='Ay')
    plt.plot(states[3], states[2][2], label='Az')
    plt.title("acceleration vs time")
    plt.legend(loc='lower center', ncol=1, bbox_to_anchor=(0.5, -0.7))

    plt.subplots_adjust(left=0.1,
                        bottom=0.2,
                        right=0.95,
                        top=0.9,
                        wspace=0.3,
                        hspace=0.4)
    plt.title("reference drop" + str(title))
    plt.show()

def PlotScatter(Drops,Stats, mancase, concase, poly1, poly2):
    points = Stats["bounds"][-1]
    cmap = plt.get_cmap('hot')

    # possible landing (combinatory uncertainties)
    plt.scatter(Drops["impact DX [m]"][1:] - Drops["impact DX [m]"][0],
                Drops["impact DY [m]"][1:] - Drops["impact DY [m]"][0],
                c=Drops["impact velocity [m/s]"][1:], ec='k', cmap=cmap,label='possible', vmin=-V_boxhit_lim, vmax=0)
    plt.colorbar(label='impact velocity [m/s]')

    # reference drop
    plt.scatter(0,0,
                marker='X', label='expected')
    # average drop
    plt.scatter(Stats["Xavg"][-1] - Drops["impact DX [m]"][0],
                Stats["Yavg"][-1] - Drops["impact DY [m]"][0],
                marker='X',label='average')

    # correct overlap
    plt.plot(*poly1.exterior.xy, color='b', label='required zone')
    plt.plot(*poly2.exterior.xy, color='r', label='landing bounds')

    plt.legend(loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.3))
    plt.subplots_adjust(left=0.1,
                        bottom=0.2,
                        right=0.95,
                        top=0.8,
                        wspace=0.3,
                        hspace=0.4)
    mc_r = []
    for var in mancase: mc_r.append(round(var))
    cc_r = []
    for var in concase: cc_r.append(round(var))
    plt.title("landing spots for:" +
            "\ncondition:"+ str(cc_r) +
            " & maneuver:"+ str(mc_r) +
            "\nSfit=" + str(round(Stats["Dfit"][-1])) + "%")
    plt.show()

# INITIALISE

# get defaults design values
AC = UAV('aircraft')
print("AC default:", AC.__dict__)

atm = atmosphere()
print("atmo default:", atm.__dict__)

box = BOX(1)
print("box default:", box.__dict__)

# Simulation param
dt_sim = 1 / 10E3  # [s] between frames
IT_max = 10E4

# Operational limits
V_boxdrop_lim = AC.OP_V_boxlim = 100/6.5 # [m/s] 100kph drop limit
V_boxhit_lim = AC.OP_Vbox_LDG = 40/6.5 # [m/s] 100kph drop limit
b_zone = 25 # [m]
h_zone = 25 # [m]
dropzone = np.array([
    [-b_zone/2,-h_zone/2],
    [-b_zone/2, h_zone/2],
    [ b_zone/2, h_zone/2],
    [ b_zone/2,-h_zone/2]
])
dropzone_poly = Polygon(dropzone)
Hmin = AC.OP_hmin # [m]

# Design limits
AC_V_s = AC.V_s_min
AC_nmax = 3

# Drop Conditions
# limits
C_Mbox = [10, 20] # [kg]
C_Vw = [0, AC.OP_V_wind] # [kg]
C_w_heading = [0, 180] # [deg]
C_Nbox = 1 # [deg]
# c combinations
N_div_C = 2
C_VAR = {
    "Mbox": subdivide(C_Mbox, N_div_C-1),
    "Vw": subdivide(C_Vw, N_div_C),
    "w_heading": subdivide(C_w_heading, N_div_C),
    "box_number": range(C_Nbox)
}
allNames = C_VAR
C_COMB = list(it.product(*(C_VAR[Name] for Name in allNames)))
print("try conditions:", C_COMB)

# Maneuver
# m limits
M_V_app = [AC_V_s, V_boxdrop_lim] # [m/s]
M_n_app = [0.5, AC_nmax] # [g]
M_pitch = [-45, 5] # [deg]
M_heading = [0,0] # [deg] not actually a free var
M_Hmin = [Hmin, Hmin+1] # [m]
# m combinations
N_div_M = 0
M_VAR = {
    "V_app": subdivide(M_V_app, N_div_M-1),
    "n_app": subdivide(M_n_app, N_div_M-1),
    "pitch": subdivide(M_pitch, N_div_M),
    "heading": subdivide(M_heading, 0),
    "Hmin": subdivide(M_Hmin, 0)
}
allNames = M_VAR
M_COMB = list(it.product(*(M_VAR[Name] for Name in allNames)))
print("try maneuvers:", M_COMB)

# Uncertainties
# limits
U_Mbox = [-0.5,0.5] # [kg]
U_T_drop = [-1,1] # [s]
U_T_brake = [-0.5,0.5] # [s]
U_T_flap = [-0.5,0.5] # [s]
U_pos = [-5,5] # [m]
U_Vw = [-5,5] # [m/s]
U_w_heading = [-25,25] # [deg]
# U combinations
N_div_U = 2
U_VAR = {
    "Mbox": subdivide(U_Mbox, 0),
    "Vw": subdivide(U_Vw, N_div_U-1),
    "w_heading": subdivide(U_w_heading, N_div_U),
    "T_drop": subdivide(U_T_drop, 0),
    "T_brake": subdivide(U_T_brake, 1),
    "T_flap": subdivide(U_T_flap, 1),
    "posX": subdivide(U_pos, 0),
    "posY": subdivide(U_pos, 0)
}
allNames = U_VAR
U_COMB = list(it.product(*(U_VAR[Name] for Name in allNames)))
print("with uncertainty:", U_COMB)

Ntot = len(M_COMB)*len(C_COMB)*len(U_COMB)
print("total simulations:", len(M_COMB), len(C_COMB), len(U_COMB), "=", Ntot, "in about", (187/384)*Ntot/60, "[min]")
time_start = time.time()

# START TRYING COMBINATIONS
print("\n================= ================= =================",
      len(C_COMB), "Drop Conditions Combinations")
ConditionResults = {
    "best M": [],
    "avg Vimpact": [],
    "avg Sfit": [],
    "worst Sfit": []
}
N_Con = 0
for DropCon in C_COMB:
    print("\n================= =================")
    print(N_Con, "/", len(C_COMB), "drop condition", str(DropCon))
    ManResults = {
        "Xavg": [],
        "Yavg": [],
        "DXdev": [],
        "DYdev": [],
        "Dfit": [],
        "Sfit": [],
        "worst V_impact": [],
        "PVi95": [], # counts boxes within and out of landing speed limit
        "worst a_max": [],
        "Racc": [],
        "bounds": []
    }
    N_Man = 0
    for ManCase in M_COMB:
        print("\n=================")
        print(N_Man, "/", len(M_COMB), "maneuver", str(ManCase), "\n")
        DropResults = {
            "impact DX [m]": [],
            "impact DY [m]": [],
            "impact velocity [m/s]": [],
            "max acceleration [g]": [],
            "max velocity [m/s]": [],
            "time to land [s]": []
        }
        N_Unc = 0

        # simulate reference drop
        DropUncRef = np.zeros(len(U_COMB[0]))
        state = SimuDrop(DropCon, ManCase, DropUncRef, DropResults)
        #PlotTraj(state, ManCase)

        # simulate uncertainty drops
        Racc = []
        Pvi = 0
        for DropUnc in U_COMB:
            print(N_Unc, "/", len(U_COMB), "variation", str(DropUnc))
            SimuDrop(DropCon, ManCase, DropUnc, DropResults)
            Racc.append(np.sqrt(DropResults["impact DX [m]"][-1] ** 2 + DropResults["impact DY [m]"][-1] ** 2))
            if DropResults["impact velocity [m/s]"][-1] < V_boxhit_lim:
                Pvi += 1
            N_Unc += 1

        # get maneuver evaluation
        points = np.array(list(zip( DropResults["impact DX [m]"] - DropResults["impact DX [m]"][0],
                                    DropResults["impact DY [m]"] - DropResults["impact DY [m]"][0])))
        mpt = MultiPoint(points)
        hull_poly = shapely.wkt.loads(mpt.convex_hull.wkt)
        if hull_poly.intersection(dropzone_poly):
            Sfit = (hull_poly.intersection(dropzone_poly).area/hull_poly.area)*100
            Dfit = (dropzone_poly.intersection(hull_poly).area / dropzone_poly.area) * 100
        else: Sfit = Dfit =0
        ManResults["Dfit"].append(Dfit)
        ManResults["Sfit"].append(Sfit)

        ManResults["PVi95"].append(100*Pvi/N_Unc)
        ManResults["Racc"].append(max(Racc))
        ManResults["bounds"].append(points)

        ManResults["Xavg"].append(
            np.mean(DropResults["impact DX [m]"]))
        ManResults["Yavg"].append(
            np.mean(DropResults["impact DY [m]"]))
        ManResults["DXdev"].append(
            max(DropResults["impact DX [m]"])-DropResults["impact DX [m]"][0])
        ManResults["DYdev"].append(
            max(DropResults["impact DY [m]"])-DropResults["impact DY [m]"][0])
        ManResults["worst a_max"].append(
            max(DropResults["max acceleration [g]"]))
        ManResults["worst V_impact"].append(
            max(DropResults["impact velocity [m/s]"]))

        PlotScatter(DropResults, ManResults, ManCase, DropCon, dropzone_poly, hull_poly)
        N_Man += 1

    # get case evaluation, find best maneuver
    print(ManResults)
    N_Con += 1

# Mparam = f(Case param)
# CaseResults = f(Case param)
time_end = time.time()
print(Ntot, "simulations took:",
      (time_end-time_start), "=[s] or",
      (time_end-time_start)/3600,  "[h]")
