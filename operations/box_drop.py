import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV,atmosphere
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import time

class box:
    def __init__(self, rel_pos, AC):
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
        self.Pos_AC = np.array([0,0,0]) #[m] relative to AC c.g.

        self.t = 0 # [s] time after hatch opening

        # initial state
        self.Pos = np.array(AC.pos) + self.Pos_AC
        self.Vel = np.array([
            AC.OP_drop_V * np.cos(AC.OP_drop_angle),
            0,
            AC.OP_drop_V * np.sin(AC.OP_drop_angle),
        ])
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
        self.Acc = np.array([
            -np.sign(self.Vel[0]) * (1 / 2) * self.CD * rho * self.S * (self.Vel[0] ** 2) / self.mass,
            -np.sign(self.Vel[1]) * (1 / 2) * self.CD * rho * self.S * (self.Vel[1] ** 2) / self.mass,
            -np.sign(self.Vel[2]) * (1 / 2) * self.CD * rho * self.S * (self.Vel[2] ** 2) / self.mass - g0
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

def AC_trajectory(aircraft):
    pass

def drop_condition(aircraft):
    pass

def drop_state(aircraft):
    pass


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

if True:
    AC = UAV('aircraft')
    print(AC.__dict__)

    atm = atmosphere()

    # simulation param
    dt_sim = 1/10E2 # [s] between frames
    IT_max = 10E4

    # init parameters
    AC.n_drops = 1  # [-]
    AC.n_boxes = 1  # [-]
    AC.PL_per_box = 10 # [kg]

    AC.OP_app_angle = 0  # [rad]       # approach
    AC.OP_app_V = 25  #[m/s]
    AC.OP_exit_angle = 0  # [rad]     # exit
    AC.OP_exit_Dn = 1  # [g0]
    AC.OP_drop_V = 25  # [m/s]       # drop
    AC.OP_drop_angle = 45*(np.pi/180)  # [rad]
    AC.OP_drop_Dn = 0  # [g0]

    #init disturbances
    H_min = 15  # [m]
    V_crosswind = 0  # [m/s]
    V_tailwind = 0  # m/s]
    V_headwind = 0  # [m/s]

    # calc AC traj s.t. h_min == 15m given approach param in AC_pos(t), AC_vel(t)
    AC.pos = [0, 0, H_min]

    # find condition point : drop order

    # get drop state for (each) box with drop order and box layout
    box = box([1, 1], AC)
    box.mass += 20 # payload
    print(box.__dict__)
    N = 0

    # get trajectory of (each) box from initial state & CD(t)
    while box.Pos[-1] > 0 and N < IT_max: #until touchdown
        box.update(dt_sim, atm)
        N += 1

    # store results
    if N + 1 == IT_max: print("SIMULATION ABORTED")
    else:
        print("timesteps:", N)
        print("impact DX:", box.LogState[0][0][-1], "[m]")
        print("impact DY:", box.LogState[0][1][-1], "[m]")
        print("impact velocity:", box.LogState[1][-1][-1], "[m/s]")
        print("max acceleration:", np.max(box.LogState[2][-1])/atm.g, "[g0]")
        print("time to end:", box.LogState[-1][-1], "[s]")

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