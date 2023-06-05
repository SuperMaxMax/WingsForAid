import sys

import matplotlib.collections

sys.path.append("..")

# Start your import below this
from parameters import UAV,atmosphere
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mlpatch
import matplotlib.colorbar as cbar
import time
import itertools as it

class BOX:
    def __init__(self, boxN, AC):
        # parameters
        self.CD = 0
        self.CD0 = 	2.1 #[-]
        self.DCD_brake = 0.5 #[-]
        self.DCD_flaps = 1.28 #[-]
        self.DT_brake = 0.5 # [s]
        self.DT_flaps = 2  # [s]
        self.T_deployment = 0.5 # [s]

        self.S = 0.4*0.4 # [m2]

        self.mass = 3 # [kg]
        self.compute_PosAC(boxN)

        self.t = 0 # [s] time after hatch opening

        # initial state
        self.Pos = np.array(AC.pos) + self.Pos_AC
        self.Vel = np.array([
            AC.OP_drop_V * np.cos(AC.OP_drop_angle),
            0,
            AC.OP_drop_V * np.sin(AC.OP_drop_angle),
        ])
        self.AirVel = self.Vel
        self.Acc = np.array([0, 0, (1 + AC.OP_drop_Dn) * atm.g])

        self.LogState = [[[self.Pos[0]], [self.Pos[1]], [self.Pos[2]]],
                         [[self.Vel[0]], [self.Vel[1]], [self.Vel[2]], [self.magnitude(self.Vel)]],
                         [[self.Acc[0]], [self.Acc[1]], [self.Acc[2]], [self.magnitude(self.Acc)]],
                         [self.t]]

    def magnitude(self, vector):
        mag = np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
        return mag
    def update(self, dt, atm):
        # get parameters
        rho = atm.rho0
        g0 = atm.g
        self.t += dt
        self.compute_CD()

        # compute
        self.AirVel = np.array([
            self.Vel[0] - atm.V_tailwind,
            self.Vel[1] - atm.V_crosswind,
            self.Vel[2]
        ])

        self.Acc = np.array([
            -np.sign(self.AirVel[0]) * (1 / 2) * self.CD * rho * self.S * (self.AirVel[0] ** 2) / self.mass,
            -np.sign(self.AirVel[1]) * (1 / 2) * self.CD * rho * self.S * (self.AirVel[1] ** 2) / self.mass,
            -np.sign(self.AirVel[2]) * (1 / 2) * self.CD * rho * self.S * (self.AirVel[2] ** 2) / self.mass - g0
        ])
        self.Vel = np.add(self.Vel, self.Acc * dt)
        self.Pos = np.add(self.Pos, self.Vel * dt)

        # append
        self.LogState[0][0].append(self.Pos[0])
        self.LogState[0][1].append(self.Pos[1])
        self.LogState[0][2].append(self.Pos[2])
        self.LogState[1][0].append(self.Vel[0])
        self.LogState[1][1].append(self.Vel[1])
        self.LogState[1][2].append(self.Vel[2])
        self.LogState[1][-1].append(self.magnitude(self.Vel))
        self.LogState[2][0].append(self.Acc[0])
        self.LogState[2][1].append(self.Acc[1])
        self.LogState[2][2].append(self.Acc[2])
        self.LogState[2][-1].append(self.magnitude(self.Acc))
        self.LogState[3].append(self.t)

    def compute_CD(self):
        self.CD = self.CD0
        self.S = 0.4**2
        if self.t > self.DT_brake:
            if self.t < self.DT_brake + self.T_deployment:
                ease = (self.t-self.DT_brake)/self.T_deployment
            else: ease = 1
            self.CD += self.DCD_brake * ease

        if self.t > self.DT_flaps:
            if self.t < self.DT_flaps + self.T_deployment:
                ease = (self.t-self.DT_flaps)/self.T_deployment
            else: ease = 1
            self.CD += self.DCD_flaps * ease
            self.S += 4*(0.6*0.4) * ease

    def compute_PosAC(self,boxN):
        dx = 0.2+0.4/2
        dy = 0.1+0.4/2
        dz = 0
        box_pos = {
            1: [0 * dx, dy, dz],
            2: [0 * dx, -dy, dz],
            3: [5 * dx, dy, dz],
            4: [5 * dx, -dy, dz],
            5: [1 * dx, dy, dz],
            6: [1 * dx, -dy, dz],
            7: [4 * dx, dy, dz],
            8: [4 * dx, -dy, dz],
            9: [2 * dx, dy, dz],
            10: [2 * dx, -dy, dz],
            11: [3 * dx, dy, dz],
            12: [3 * dx, -dy, dz],
        }
        self.Pos_AC = np.array(box_pos[boxN])  # [m] relative to AC c.g.

# def print_name(obj):
#     # print aircraft name
#     print(obj.name)
#
# def overwrite_value(obj):
#     # overwrite aircraft name
#     obj.name = 'new name'
#
# print_name(aircraft)
# overwrite_value(aircraft)
# print_name(aircraft)

def subdivide(range,N):
    var = []
    if N <= 0:
        var.append(np.mean(range))
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

if True:
    AC = UAV('aircraft')
    print(AC.__dict__)
    print("\n")

    atm = atmosphere()

    # simulation param
    dt_sim = 1/10E3 # [s] between frames
    IT_max = 10E4
    plot_traj = False
    plot_scatter = True
    plot_bounds = True
    N_divisions_maneuver = 1
    N_divisions_disturbance = 1
    N_divisions_timing = 0

    # global parameters
    AC.n_drops = 1  # [-]
    AC.n_boxes = 1  # [-]

    # drop parameters limits
    AC.OP_app_V = 25  #[m/s]            # approach
    AC.OP_app_Dn = 1  # [g0]
    AC.OP_drop_V = 25  # [m/s]          # drop
    AC.OP_drop_angle = 45*(np.pi/180)  # [rad]
    AC.OP_drop_Dn = 0  # [g0]
    AC.OP_exit_Dn = 1  # [g0]           # exit
    H_min = 15  # [m]

    # disturbances limits
    atm.V_crosswind = 10  # [m/s]
    atm.V_headwind = 10  # [m/s]
    atm.V_tailwind = -atm.V_headwind  # m/s]
    AC.PL_per_box = 20 # [kg]
    AC.minW_per_box = 10 # [kg]
    time_dev_drop = 0.5 # [s]
    time_dev_brake = 0.5  # [s]
    time_dev_flap = 1  # [s]

    # START CASE EVALUATION

    # Generate drop maneuver variations
    Maneuver = {
        "approach speed": subdivide([20, AC.V_s_min], N_divisions_maneuver), # add actual minimum stall speed
        "approach load factor": subdivide([0, 3], N_divisions_maneuver),
        "drop angle": subdivide([-45, 15], N_divisions_maneuver)
    }

    VarResults = {
        "mean DX [m]": [],
        "mean DY [m]": [],
        "deviation DX [m]": [],
        "deviation DY [m]": [],
        "worst-case impact velocity [m/s]": [],
        "worst-case acceleration [g]": []
    }

    allNames = sorted(Maneuver)
    Mancombinations = list(it.product(*(Maneuver[Name] for Name in allNames)))
    print("================= =================\n", len(Mancombinations), "maneuvers combinations\n")
    #print(Mancombinations)

    Maneuver_ID = 1
    for maneuver_case in Mancombinations:
        print("\n", Maneuver_ID, "/", len(Mancombinations),
              "maneuvers case ", str(maneuver_case))
        AC.OP_app_V = maneuver_case[0]
        AC.OP_drop_Dn = maneuver_case[1]
        AC.OP_drop_angle = maneuver_case[2]

        # calc AC traj s.t. h_min == 15m
        Dh = 0
        if AC.OP_drop_angle != 0:
            R = AC.OP_app_V/AC.OP_app_Dn # approach pull radius
            Dh = R * np.abs(np.sin(AC.OP_drop_angle))
        AC.pos = [0, 0, H_min+Dh]

        # Generate disturbance variations

        DropResults = {
            "impact DX [m]": [],
            "impact DY [m]": [],
            "impact velocity [m/s]": [],
            "impact deceleration [g]": [],
            "max acceleration [g]": [],
            "max velocity [m/s]": [],
            "time to land [s]": []
        }

        Disturbances = {
            "payload mass [kg]": subdivide([AC.minW_per_box, AC.PL_per_box], N_divisions_disturbance),
            "crosswind [m/s]": subdivide([-atm.V_crosswind, atm.V_crosswind], N_divisions_disturbance),
            "headwind [m/s]": subdivide([-atm.V_headwind, atm.V_headwind], N_divisions_disturbance),
            "timing deviation drop [s]": subdivide([-time_dev_drop ,time_dev_drop],N_divisions_timing),
            "timing deviation brake [s]": subdivide([-time_dev_brake, time_dev_brake], N_divisions_timing),
            "timing deviation flaps [s]": subdivide([-time_dev_flap, time_dev_flap], N_divisions_timing)
        }

        allNames = sorted(Disturbances)
        Varcombinations = list(it.product(*(Disturbances[Name] for Name in allNames)))
        print("\n=================\n", len(Varcombinations), "disturbance combinations")
        #print(Varcombinations)

        Disturbance_ID = 1
        for disturbance_case in Varcombinations:
            print(Disturbance_ID, "/", len(Varcombinations),
                  "disturbance case ", str(disturbance_case))
            atm.V_crosswind = disturbance_case[0]  # [m/s]
            atm.V_headwind = disturbance_case[1]  # [m/s]
            atm.V_tailwind = -atm.V_headwind  # m/s]
            AC.PL_mass = disturbance_case[2] # [kg]
            time_dev_brake = disturbance_case[4]  # [s]
            time_dev_drop = disturbance_case[3] # [s]
            time_dev_flap = disturbance_case[5]  # [s]

            # disturb drop time
            if time_dev_drop != 0:
                AC.OP_drop_angle += (atm.g / AC.OP_app_Dn) * (AC.OP_app_Dn - 1) * time_dev_drop
                AC.pos = np.array([0, 0, H_min + Dh])
                + np.array(
                    [np.cos(AC.OP_drop_angle), np.sin(AC.OP_drop_angle), 0]
                ) * AC.OP_app_V * time_dev_drop

            # get drop state for (each) box with drop order and box layout
            box_number = 1
            while box_number <= AC.n_boxes:

                # initialise box pos & vel
                box = BOX(box_number, AC)
                box.mass += AC.PL_mass # payload

                # disturb box timings
                box.DT_brake += time_dev_brake
                box.DT_flaps += time_dev_flap

                # get trajectory of (each) box from initial state & CD(t)
                N = 0
                while box.Pos[-1] > 0 and N < IT_max: #until touchdown
                    box.update(dt_sim, atm)
                    N += 1

                # get ground impact deceleration
                box.crumple_size = 0.5 # [m]
                V_impact = box.LogState[1][-1][-1]
                a = V_impact**2/box.crumple_size

                # store results
                if N + 1 == IT_max: print("SIMULATION ABORTED")
                else:
                    DropResults["impact DX [m]"].append(box.LogState[0][0][-1])
                    DropResults["impact DY [m]"].append(box.LogState[0][1][-1])
                    DropResults["impact velocity [m/s]"].append(box.LogState[1][2][-1])
                    DropResults["impact deceleration [g]"].append(a / atm.g)
                    DropResults["max acceleration [g]"].append(np.max(box.LogState[2][-1])/atm.g)
                    DropResults["max velocity [m/s]"].append(np.max(box.LogState[1][-1]))
                    DropResults["time to land [s]"].append(box.LogState[-1][-1])

                # show DropResults
                if plot_traj == box_number:
                    # print results
                    print("\ndisturbance case ", str(disturbance_case), "box ", str(box_number), ":")
                    for var in DropResults:
                        print(var, DropResults[var][-1])

                    # plot
                    plt.suptitle("Box drop")

                    plt.subplot(231)
                    plt.plot(box.LogState[0][0], box.LogState[0][1])
                    plt.title("trajectory XY")

                    plt.subplot(232)
                    plt.plot(box.LogState[0][0], box.LogState[0][2])
                    plt.title("trajectory XZ")

                    plt.subplot(233)
                    plt.plot(box.LogState[0][1], box.LogState[0][2])
                    plt.title("trajectory YZ")

                    plt.subplot(234)
                    plt.plot(box.LogState[-1], box.LogState[0][-1])
                    plt.title("height over time")

                    plt.subplot(235)
                    plt.plot(box.LogState[-1], box.LogState[1][-1])
                    plt.title("speed vs time")

                    plt.subplot(236)
                    plt.plot(box.LogState[-1], box.LogState[2][-1])
                    plt.title("acceleration vs time")

                    plt.show()

                # iterate
                box_number += 1

            Disturbance_ID += 1

        # process results
        VarResults["mean DX [m]"].append(
            np.mean(DropResults["impact DX [m]"]))
        VarResults["mean DY [m]"].append(
            np.mean(DropResults["impact DY [m]"]))
        VarResults["deviation DX [m]"].append(
            abs(max(DropResults["impact DX [m]"])-min(DropResults["impact DX [m]"])))
        VarResults["deviation DY [m]"].append(
            abs(max(DropResults["impact DY [m]"]) - min(DropResults["impact DY [m]"])))
        VarResults["worst-case impact velocity [m/s]"].append(
            max(DropResults["impact velocity [m/s]"]))
        VarResults["worst-case acceleration [g]"].append(
            max(DropResults["impact deceleration [g]"])) # ignore ground impact accel as indicated by impact speed already
        # calculate expected spread TODO

        # show VarResult
        if plot_scatter:
            # print results
            print("\ndrop case with", len(Varcombinations), "combinations and", AC.n_boxes, "boxes :")
            for var in VarResults:
                print(var, VarResults[var][-1])

            # landing scatter
            cmap = plt.get_cmap('hot')
            plt.scatter(DropResults["impact DX [m]"]-VarResults["mean DX [m]"][-1],
                        DropResults["impact DY [m]"]-VarResults["mean DY [m]"][-1],
                        c=DropResults["impact velocity [m/s]"], ec='k', cmap=cmap)
            plt.title("landing spots for maneuver "+str(maneuver_case))
            plt.colorbar()
            plt.show()

        Maneuver_ID += 1

    if plot_bounds:

        # variables to plot
        vs   = VarResults["worst-case impact velocity [m/s]"]
        xavg = VarResults["mean DX [m]"]
        yavg = VarResults["mean DY [m]"]
        xdev = VarResults["deviation DX [m]"]
        ydev = VarResults["deviation DY [m]"]
        REQ_LDG = 25  # [m] TODO link  to val

        # initialize plot area
        fig, ax = plt.subplots(1)
        Rlim = 0.6*max(max(xdev), max(ydev), REQ_LDG)
        plt.ylim(-Rlim, Rlim)
        plt.xlim(-Rlim, Rlim)
        handles, labels = ax.get_legend_handles_labels()

        # add REQ limit
        rect = mlpatch.Rectangle((-REQ_LDG / 2, -REQ_LDG / 2), REQ_LDG, REQ_LDG,
                                 edgecolor='black',
                                 linewidth=1,
                                 hatch='//', fill=False,
                                 label="REQ limit"
                                 )
        ax.add_patch(rect)
        handles.append(rect)

        # initialize colormap
        N = len(vs)
        colors = (vs+max(vs))/min(vs) # negative normal to 0-1
        normal = plt.Normalize(min(vs), max(vs))
        cmap = plt.get_cmap('hot')
        c = cmap(colors)
        print(cmap)
        print(colors)
        print(c)

        # draw landing zones
        for i in range(N):
            rect = mlpatch.Rectangle((-xdev[i]/2,-ydev[i]/2),xdev[i],ydev[i],
                                    edgecolor=c[i],
                                    linewidth = 2,
                                    fill=False,
                                    label = str(Mancombinations[i])
                                    )
            ax.add_patch(rect)
            handles.append(rect)

        cax, _ = cbar.make_axes(ax)
        cb2 = cbar.ColorbarBase(cax, cmap=cmap, norm=normal)

        # add legend
        plt.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0, -0.15),
                   ncol=3)

        plt.show()