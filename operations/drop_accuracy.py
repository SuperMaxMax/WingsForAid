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
import winsound

# FUNCTIONS
def SimuDrop(DropCon, DropMan, DropUnc, DropResults):
    # TODO adjust box param s.t. it lands as expected with default maneuver( current design)
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
                    0,
                    np.sin(Approach["pitch"])
                ]) * Approach["V_app"] * time_dev_drop

    # apply uncertainties to box
    newbox = BOX(DropCon[-1])
    newbox.mass += DropCon[0] + DropUnc[0]
    newbox.DT_brake += DropUnc[4]
    newbox.DT_flaps += DropUnc[5]

    # get box position relative to AC
    boxDX = box_layout[0]
    boxDY = box_layout[1]
    boxDZ = box_layout[2]
    box_pos = {
        0: [0, 0, 0], # reference
        1: [0 * boxDX, boxDY, boxDZ], #diagonal front left
        2: [5 * boxDX, -boxDY, boxDZ], #diagonal back right
        3: [0 * boxDX, -boxDY, boxDZ],
        4: [5 * boxDX, boxDY, boxDZ],
        5: [1 * boxDX, boxDY, boxDZ], #diagonal front left
        6: [4 * boxDX, -boxDY, boxDZ], #diagonal back right
        7: [1 * boxDX, -boxDY, boxDZ],
        8: [4 * boxDX, boxDY, boxDZ],
        9: [2 * boxDX, boxDY, boxDZ], #diagonal front left
        10: [3 * boxDX, -boxDY, boxDZ], #diagonal back right
        11: [2 * boxDX, -boxDY, boxDZ],
        12: [3 * boxDX, boxDY, boxDZ],
    }
    box_pos = np.array(box_pos[DropCon[-1]])

    # get initial state
    rho = atm.rho0
    g0 = atm.g
    N = t = 0

    Pos = AC_pos + box_pos
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

    #plt.suptitle("Box drop")
    plt.suptitle("reference drop" + str(title))

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
    plt.show()

def PlotScatter(Drops,Stats, mancase, concase, poly1, poly2, score=0, plot=True):
    cmap = plt.get_cmap('hot')
    plt.grid(visible=None, which='major', axis='both', color='black', linestyle='-', linewidth=0.1)

    # possible landing (combinatory uncertainties)
    Vmax = 30 # AC.OP_Vbox_LDG
    plt.scatter(Drops["impact DX [m]"][1:] - Drops["impact DX [m]"][0],
                Drops["impact DY [m]"][1:] - Drops["impact DY [m]"][0],
                c=Drops["impact velocity [m/s]"][1:], ec='k', cmap=cmap,label='possible', vmin=-Vmax, vmax=0)
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
                        bottom=0.1,
                        right=0.95,
                        top=0.85,
                        wspace=0.3,
                        hspace=0.4)
    mc_r = []
    for var in mancase: mc_r.append(np.round(var, 2))
    cc_r = []
    for var in concase: cc_r.append(np.round(var, 2))
    if plot:
        plt.title("landing spots for:" +
                "\ncondition:"+ str(cc_r) +
                " & maneuver:"+ str(mc_r) +
                "\nSfit=" + str(round(Stats["Sfit"][-1])) + "%" +
                  " and Dfit=" + str(round(Stats["Dfit"][-1]))
                  + "% | Score:" + str(np.round(score,1))
        )
    else:
        plt.title("landing spots for all conditions, single maneuver, all uncertainties")

    plt.xlabel('length [m]')
    plt.ylabel('width [m]')
    plt.show()

def PlotBounds(ConDrops, droppoly):
    cmap = plt.get_cmap('hot')
    plt.grid(visible=None, which='major', axis='both', color='black', linestyle='-', linewidth=0.1)

    for cat, var in ConDrops.items():
        NDrops = len(var)
    NDrops -= 1
    Vmax = 30
    while NDrops >= 0:
        # average drop
        plt.scatter(ConDrops["Xavg"][NDrops] - ConDrops["Xref"][NDrops],
                    ConDrops["Yavg"][NDrops] - ConDrops["Yref"][NDrops],
                    marker='x',
                    c=ConDrops["Vimax"][NDrops], ec='k', cmap=cmap, vmin=-Vmax, vmax=0
                    )
        # bounds
        points = ConDrops["bounds"][NDrops]
        mpt = MultiPoint(points)
        hull_poly = shapely.wkt.loads(mpt.convex_hull.wkt)
        plt.plot(*hull_poly.exterior.xy, linewidth=1,
                 #c=ConDrops["Vimax"][NDrops], ec='k', cmap=cmap, vmin=-Vmax, vmax=0
                 )
        NDrops -= 1

    # dropzone
    plt.plot(*droppoly.exterior.xy, color='b', label='required zone')

    # plot surronding thingies
    plt.colorbar(label='worst impact velocity [m/s]')
    plt.legend(loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.3))
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.95,
                        top=0.85,
                        wspace=0.3,
                        hspace=0.4)
    plt.title("landing spots for all conditions, all uncertainties")
    plt.xlabel('length [m]')
    plt.ylabel('width [m]')
    plt.show()

# INITIALISE

def PlotBar(ManTPM, TPMw, Labels, Title, scores):

    labels = []
    for lab in Labels:
        lab_name = str()
        for val in lab:
            lab_name += str(round(val)) + " "
        labels.append(lab_name)

    fig, ax = plt.subplots()
    N = len(labels)
    bottom = np.zeros(N)
    width = 0.5

    w_counter = 0
    for boolean, tpm in ManTPM.items():
        tpm = np.array(tpm)
        w = TPMw[w_counter]
        bar = tpm * w
        p = ax.bar(labels, bar, width, bottom=bottom,
                   label=str(boolean+"|"+str(np.round(w, 2))))
        bottom += bar
        w_counter += 1

    dy = 0.25
    ax.set_title("TPM per case for condition " + str(Title))
    plt.legend(loc='lower center', ncol=3, bbox_to_anchor=(0.5, 1+dy))
    plt.subplots_adjust(left=0.1,
                        bottom=dy,
                        right=0.95,
                        top=0.7,
                        wspace=0.3,
                        hspace=0.4)
    ax.tick_params(axis='x', labelrotation=-75)
    plt.show()

    return scores

def unique(list1): # get unique elements of list
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

# get defaults design values
AC = UAV('aircraft')
print("AC default:", AC.__dict__)

atm = atmosphere()
print("atmo default:", atm.__dict__)

box = BOX(1)
print("box default:", box.__dict__)

# Simulation param
dt_sim = 5 / 10E3  # [s] between simulation frames
IT_max = 10E4
SingleTry = True # 1 maneuver, full uncertainties
SingleTryGNC = False # 1 maneuver, all conditions with some uncertainty
DoVerif = False

# Operational limits & requirements
V_boxdrop_lim = 100/3.6 # [m/s] 100kph drop limit
V_boxhit_lim = 40/3.6 # [m/s] 40kph drop limit
box_layout = [AC.boxDX, AC.boxDY, AC.boxDZ]
Hmin = AC.OP_hmin # [m]

# Design limits
AC_V_s = AC.V_s_min
AC_nmax = AC.OP_n_app_max

# Selection criteria
cat_w = [1, 0.5, 0.1, 0.4]
TPM_weights = np.array([
    cat_w[0]*1, # REQ
    cat_w[1]*0.5, cat_w[1]*0.5,  # REQ: Sfit, Pvi
    cat_w[2]*0.2, cat_w[2]*0.8, # OPS: Racc, amax
    cat_w[3]*0.7, cat_w[3]*0.15, cat_w[3]*0.15 # AC: n,pitch,V
])

# Drop Conditions
# limits
C_Mbox = [10, 20] # [kg]
C_Vw = [0, AC.OP_V_wind] # [kg]
C_w_heading = [0, 180] # [deg]
C_Nbox = 1 # [#]
# c combinations
N_div_C = 2
if SingleTry:
    C_Nbox = 12
    C_VAR = {
        "Mbox": subdivide(C_Mbox, 0),
        "Vw": subdivide([0,0], 0),
        "w_heading": subdivide([0,0], 0),
        "box_number": [C_Nbox]
    }
elif SingleTryGNC:
    C_Nbox = 1
    C_VAR = {
        "Mbox": subdivide(C_Mbox, 1),
        "Vw": subdivide(C_Vw, 2),
        "w_heading": subdivide(C_w_heading, 4),
        "box_number": [C_Nbox]
    }
else:
    C_VAR = {
        "Mbox": subdivide(C_Mbox, min(N_div_C,1)),
        "Vw": subdivide(C_Vw, min(N_div_C,1)),
        "w_heading": subdivide(C_w_heading, N_div_C),
        "box_number": [C_Nbox]
    }
    if N_div_C == 0: # case with no distrubances, just box variations
        C_Nbox = 1
        C_VAR = {
            "Mbox": subdivide([10,10], 0), # use min box == max var
            "Vw": subdivide([0,0], 0),
            "w_heading": subdivide([0,0], 0),
            "box_number": [6]
        }
allNames = C_VAR
C_COMB = list(it.product(*(C_VAR[Name] for Name in allNames)))
C_COMB_red = []
for c_comb in C_COMB:
    c_comb = list(c_comb)
    if c_comb[1] == 0:c_comb[2] = 0 # no wing heading difference if no wind
    C_COMB_red.append(c_comb)
C_COMB = unique(C_COMB_red)
print("try conditions:", C_COMB)

# Maneuver
# m limits
M_V_app = [AC_V_s, V_boxdrop_lim] # [m/s]
M_n_app = [0, AC_nmax] # [g]
M_pitch = [-45, 0] # [deg]
M_heading = [0, 0] # [deg] not actually a free var
M_Hmin = [Hmin, 2*Hmin] # [m]
# m combinations
N_div_M = 2
M_VAR = {
    "V_app": subdivide(M_V_app, min(N_div_M, 1)),
    "n_app": subdivide(M_n_app, N_div_M),
    "pitch": subdivide(M_pitch, N_div_M),
    "heading": subdivide(M_heading, 0),
    "Hmin": subdivide(M_Hmin, N_div_M)
}
if N_div_M == 0: # use standard W4A maneuver
    M_VAR = {
        "V_app": [AC_V_s],
        "n_app": [0],
        "pitch": [0],
        "heading": [0],
        "Hmin": [Hmin]
    }
allNames = M_VAR
M_COMB = list(it.product(*(M_VAR[Name] for Name in allNames)))
M_COMB_red = []
for m_comb in M_COMB:
    m_comb = list(m_comb)
    if m_comb[1]==0 or m_comb[2]==0:
        m_comb[1] = 0 # no load factor if angle is 0
        m_comb[2] = 0
    M_COMB_red.append(m_comb)
M_COMB = unique(M_COMB_red)
if SingleTry or SingleTryGNC:
    M_COMB = [[AC_V_s, 0, 0, 0, Hmin]]
print("try maneuvers:", M_COMB)

# Uncertainties
# limits
U_dt = 0.1 # [s] default time deviation
U_Mbox = [-0.5,0.5] # [kg]
U_T_drop = [U_dt, U_dt] # [s]
U_T_brake = [-U_dt, U_dt] # [s]
U_T_flap = [-U_dt, U_dt] # [s]
U_pos = [0.25, 0.25] # [m]
U_Vw = [-5, 5] # [m/s]
U_w_heading = [-10, 10] # [deg]
# U combinations
N_div_U = 1
if SingleTry:
    U_Mbox = [-5, 5]
    U_Vw = [-AC.OP_V_wind, AC.OP_V_wind]
    U_w_heading = [10, 170]
    U_VAR = {
        "Mbox": subdivide(U_Mbox, 1),
        "Vw": subdivide(U_Vw, 2),
        "w_heading": subdivide(U_w_heading, 8),
        "T_drop": subdivide(U_T_drop, 0),
        "T_brake": subdivide(U_T_brake, 1),
        "T_flap": subdivide(U_T_flap, 1),
    }
else:
    U_VAR = {
        "Mbox": subdivide(U_Mbox, 1),
        "Vw": subdivide(U_Vw, 1),
        "w_heading": subdivide(U_w_heading, N_div_U+1),
        "T_drop": subdivide(U_T_drop, 0),
        "T_brake": subdivide(U_T_brake, 0),
        "T_flap": subdivide(U_T_flap, 0),
    }
allNames = U_VAR
U_COMB = list(it.product(*(U_VAR[Name] for Name in allNames)))
U_COMB = unique(U_COMB)
print("with uncertainty:", U_COMB)

# re-assess dropzone to include position inaccuracy
b_zone = AC.OP_b_dropzone - U_pos[0] # [m]
h_zone = AC.OP_l_dropzone - U_pos[1] # [m]
dropzone = np.array([
    [-b_zone/2,-h_zone/2],
    [-b_zone/2, h_zone/2],
    [ b_zone/2, h_zone/2],
    [ b_zone/2,-h_zone/2]
])
dropzone_poly = Polygon(dropzone)

# Trial summary
Ntot = len(M_COMB)*len(C_COMB)*len(U_COMB)*C_Nbox
print("total simulations:", len(C_COMB), len(M_COMB), len(U_COMB), C_Nbox,"=", Ntot, "in about", ((136.966/60)/1440)*Ntot, "[min]")
time_start = time.time()

# START TRYING COMBINATIONS
print("\n================= ================= =================")
ConResults = {
    "Condition-Maneuver": [], # parameters of chosen maneuver
    "Score": [], # weighted TPM sum
    "Xavg": [],
    "Yavg": [],
    "Xref": [],
    "Yref": [],
    "bounds": [], # better bounds
    "Vimax": []
}
ConditionTPM = {
    "REQ passed": [],
    "REQ_S": [],
    "REQ_Vi": [],
    "OPS_R": [],
    "OPS_am": [],
    "AC_na": [],
    "AC_pa": [],
    "AC_Va": []
}

N_Con = 1
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
    ManeuverTPM = {
        "REQ passed": [],
        "REQ_S": [],
        "REQ_Vi": [],
        "OPS_R": [],
        "OPS_am": [],
        "AC_na": [],
        "AC_pa": [],
        "AC_Va": []
    }
    Scores = []
    TPMs = []

    N_Man = 1
    for ManCase in M_COMB:
        #print("\n=================")
        print(N_Man, "/", len(M_COMB), "maneuver", str(ManCase), "\n")
        DropResults = {
            "impact DX [m]": [],
            "impact DY [m]": [],
            "impact velocity [m/s]": [],
            "max acceleration [g]": [],
            "max velocity [m/s]": [],
            "time to land [s]": []
        }
        N_Unc = 1

        # simulate reference drop
        DropUncRef = np.zeros(len(U_COMB[0]))
        DropCon[-1] = 0
        state = SimuDrop(DropCon, ManCase, DropUncRef, DropResults)
        if SingleTry or True: PlotTraj(state, ManCase)

        # simulate drop uncertainties
        Racc = []
        Pvi = 0
        for DropUnc in U_COMB:
            DropUnc = list(DropUnc)
            if SingleTry: print(N_Unc, "/", len(U_COMB), "variation", str(DropUnc))
            DropCon[-1] = C_Nbox
            DropUnc[3] = 0
            while DropCon[-1] > 0:
                SimuDrop(DropCon, ManCase, DropUnc, DropResults)
                Racc.append(np.sqrt(DropResults["impact DX [m]"][-1] ** 2 + DropResults["impact DY [m]"][-1] ** 2))
                if DropResults["impact velocity [m/s]"][-1] < V_boxhit_lim:
                    Pvi += 1
                # try for N boxes
                DropCon[-1] -= 1
                DropUnc[3] += 0.5 # time delay between each box
            N_Unc += 1

        # get maneuver evaluation
        points = np.array(list(zip( DropResults["impact DX [m]"] - DropResults["impact DX [m]"][0],
                                    DropResults["impact DY [m]"] - DropResults["impact DY [m]"][0])))
        mpt = MultiPoint(points)
        hull_poly = shapely.wkt.loads(mpt.convex_hull.wkt)
        if hull_poly.intersection(dropzone_poly):
            Sfit = (hull_poly.intersection(dropzone_poly).area/hull_poly.area)
            Dfit = (dropzone_poly.intersection(hull_poly).area / dropzone_poly.area)
        else:
            Sfit = 0
            Dfit = 0
        ManResults["Dfit"].append(Dfit)
        ManResults["Sfit"].append(Sfit)

        Pvi *= 1/N_Unc
        ManResults["PVi95"].append(Pvi)
        ManResults["Racc"].append(max(Racc))
        ManResults["bounds"].append(points)

        ManResults["Xavg"].append(
            np.mean(DropResults["impact DX [m]"]))
        ManResults["Yavg"].append(
            np.mean(DropResults["impact DY [m]"]))
        ManResults["DXdev"].append(
            DropResults["impact DX [m]"][0])
        ManResults["DYdev"].append(
            DropResults["impact DY [m]"][0])
        ManResults["worst a_max"].append(
            max(DropResults["max acceleration [g]"]))
        ManResults["worst V_impact"].append(
            max(DropResults["impact velocity [m/s]"]))

        # normalize & prep TPM
        if Sfit >= 100:
            ManeuverTPM["REQ passed"].append(1)
        if Pvi >= 95:
            ManeuverTPM["REQ passed"].append(1)
        else:
            ManeuverTPM["REQ passed"].append(0)
        ManeuverTPM["REQ_S"].append(round(Sfit))
        ManeuverTPM["REQ_Vi"].append(round(Pvi))
        ManeuverTPM["OPS_R"].append(round(1-Dfit))
        ManeuverTPM["OPS_am"].append(1-(ManResults["worst a_max"][-1]/AC_nmax))
        ManeuverTPM["AC_na"].append(1 - (ManCase[1] / AC_nmax))
        ManeuverTPM["AC_pa"].append(1 - (abs(ManCase[2]) / max(abs(M_pitch[0]), abs(M_pitch[1]))))
        ManeuverTPM["AC_Va"].append((ManCase[0]-AC_V_s)/AC_V_s)

        TPM = []
        for boolean, tpm in ManeuverTPM.items():
            TPM.append(100*max(tpm[-1], 0))
        TPMs.append(TPM)

        Score = np.sum(np.vdot(
            np.array(TPM),
            np.array(TPM_weights)))
        Scores.append(Score)
        if SingleTry or SingleTryGNC:
            PlotScatter(DropResults, ManResults, ManCase, DropCon, dropzone_poly, hull_poly, Score)

        N_Man += 1

    index = np.argsort(Scores)

    ConResults["Condition-Maneuver"].append([DropCon, M_COMB[index[-1]]])
    ConResults["Score"].append(Scores[index[-1]])
    ConResults["Vimax"].append(ManResults["worst V_impact"][index[-1]])
    ConResults["Xavg"].append(ManResults["Xavg"][index[-1]])
    ConResults["Yavg"].append(ManResults["Yavg"][index[-1]])
    ConResults["Xref"].append(ManResults["DXdev"][index[-1]])
    ConResults["Yref"].append(ManResults["DYdev"][index[-1]])
    ConResults["bounds"].append(ManResults["bounds"][index[-1]])
    i = 0
    for boolean, tpm in ConditionTPM.items():
        ConditionTPM[boolean].append(TPMs[index[-1]][i])
        i += 1

    # bar plot
    if not (SingleTry or SingleTryGNC) and N_div_M > 0:
        #print("Results from ", DropCon, "are:\n", ConResults)
        PlotBar(ManeuverTPM, TPM_weights, M_COMB, DropCon, Scores)
    N_Con += 1

# plot ref drop overlay
# plot scatter overlay
PlotBounds(ConResults, dropzone_poly)
# plot bar (cond)
if not (SingleTry or SingleTryGNC) and N_div_C > 0 and False:
    PlotBar(ConditionTPM, TPM_weights, ConResults["Condition-Maneuver"][0], "all conditions", ConResults["Score"])

# Sensitivity & Causality
# Mparam = f(Case param)
# CaseResults = f(Case param)

time_mid = time.time()
print(Ntot, "simulations took:",
      (time_mid - time_start), "[s] or",
      (time_mid - time_start) / 60, "[min]")

if DoVerif: # optional Verification
    print("\n================= ================= =================")
    # rerun all winning maneuvers through all conditions on better resolution
    U_VAR = {
        "Mbox": subdivide(U_Mbox, 1),
        "Vw": subdivide(U_Vw, 2),
        "w_heading": subdivide(U_w_heading, 4),
        "T_drop": subdivide(U_T_drop, 1),
        "T_brake": subdivide(U_T_brake, 1),
        "T_flap": subdivide(U_T_flap, 1),
    }
    allNames = U_VAR
    U_COMB = list(it.product(*(U_VAR[Name] for Name in allNames)))
    U_COMB = unique(U_COMB)
    print(len(U_COMB), "new uncertainties")

    dt_sim = 1 / 10E3  # [s] between simulation frames
    IT_max = 10E5

    DropResults = {
        "impact DX [m]": [],
        "impact DY [m]": [],
        "impact velocity [m/s]": [],
        "max acceleration [g]": [],
        "max velocity [m/s]": [],
        "time to land [s]": []
    }
    N_Man = N_Con = 1
    for ManConCase in ConResults["Condition-Maneuver"]:
        # Rerun set maneuver on all conditions
        ManCase = ManConCase[1]
        DropCon = ManConCase[0]
        print(N_Con, "/", len(ConResults["Condition-Maneuver"]), "drop condition", str(DropCon), " & ", str(ManCase), "maneuver")

        DropResults = {
            "impact DX [m]": [],
            "impact DY [m]": [],
            "impact velocity [m/s]": [],
            "max acceleration [g]": [],
            "max velocity [m/s]": [],
            "time to land [s]": []
        }

        # simulate reference drop
        DropUncRef = np.zeros(len(U_COMB[0]))
        ConResults["states"].append(SimuDrop(DropCon, ManCase, DropUncRef, DropResults))

        # simulate drop uncertainties
        for DropUnc in U_COMB:
            print(N_Unc, "/", len(U_COMB), "variation", str(DropUnc))
            SimuDrop(DropCon, ManCase, DropUnc, DropResults)
            N_Unc += 1

        # get maneuver evaluation
        points = np.array(list(zip( DropResults["impact DX [m]"] - DropResults["impact DX [m]"][0],
                                    DropResults["impact DY [m]"] - DropResults["impact DY [m]"][0])))
        ConResults["bounds"].append(points)
        mpt = MultiPoint(points)
        hull_poly = shapely.wkt.loads(mpt.convex_hull.wkt)
        if hull_poly.intersection(dropzone_poly):
            if hull_poly.area == 0: Sfit = 100
            else:Sfit = (hull_poly.intersection(dropzone_poly).area / hull_poly.area)*100
            Dfit = (dropzone_poly.intersection(hull_poly).area / dropzone_poly.area) * 100
        else: Sfit = Dfit =0
        ConResults["Dfit"].append(Dfit)
        ConResults["Sfit"].append(Sfit)
        ConResults["Xavg"].append(
            np.mean(DropResults["impact DX [m]"]))
        ConResults["Yavg"].append(
            np.mean(DropResults["impact DY [m]"]))

        # plots
        PlotScatter(DropResults, ConResults, ManCase, DropCon, dropzone_poly, hull_poly, ConResults["Score"][N_Man]) # TODO don't draw poitns
        PlotTraj(ConResults["states"][N_Man], ManCase) # TODO overlay all ref

        N_Man += 1
        N_Con += 1

# PROGRAM END
print("\n================= ================= =================")
time_end = time.time()
frequency = 5000  # Set Frequency To 2500 Hertz
duration = 1*1000  # Set Duration To 1000 ms == 1 second
winsound.Beep(frequency, duration)
print(Ntot, "simulations took:",
      (time_end-time_start), "[s] or",
      (time_end-time_start)/60,  "[min]")
