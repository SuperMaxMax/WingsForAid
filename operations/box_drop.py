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
    def __init__(self, slot):
        self.CD = 0
        self.CD0 = 	2.1 #[-]
        self.DCD_brake = 0.5 #[-]
        self.CDC_flaps = 2 #[-]
        self.DT_brake = 1 # [s]
        self.DT_flaps = 2  # [s]

        self.S = 0.4*0.4 # [m2]

        self.mass = 3 # [kg]
        self.Pos_AC = [] #[m] relative to AC c.g.
        self.Poslog = []
        self.Pos = []
        self.Vel = []
        self.Acc = []

        self.t = 0 # [s] time after hatch opening

    def update(self, dt, atm):
        rho = atm.rho0
        g0 = atm.g

        self.t += dt
        self.CD()

        self.Acc = [
            -(1 / 2)*self.CD*rho*self.Vel[0]**2*self.S,
            -(1 / 2) * self.CD * rho * self.Vel[1] ** 2 * self.S,
            g0 - (1 / 2) * self.CD * rho * self.Vel[2] ** 2 * self.S
        ]
        self.Vel += self.Acc * dt
        self.Pos += self.Vel * dt

    def CD(self):
        self.CD = self.CD0
        if self.t > self.DT_brake: self.CD += self.DCD_brake
        if self.t > self.DT_flaps: self.CD += self.CDC_flaps

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

    box = box([1,1])
    print(box.__dict__)

    atm = atmosphere()

    # simulation param
    dt_sim = 0.1 # [s] between frames

    # init parameters
    AC.n_drops = 1  # [-]
    AC.n_boxes = 1  # [-]

    AC.OP_app_angle = 0  # [deg]       # approach
    AC.OP_app_V = 100  #[m/s]
    AC.OP_exit_angle = 0  # [deg]     # exit
    AC.OP_exit_Dn = 1  # [g0]
    AC.OP_drop_V = 100  # [m/s]       # drop
    AC.OP_drop_angle = 0  # [deg]
    AC.OP_drop_Dn = 1  # [g0]

    #init disturbances
    H_min = 15  # [m]
    V_crosswind = 0  # [m/s]
    V_tailwind = 0  # m/s]
    V_headwind = 0  # [m/s]

    # calc AC traj s.t. h_min == 15m given approach param in AC_pos(t), AC_vel(t)
    AC.pos = [0,0,H_min]

    # find condition point : drop order

    # get drop state for (each) box with drop order and box layout
    box.Pos = AC.pos
    box.Vel = [
        AC.OP_drop_V * np.cos(AC.OP_drop_angle),
        0,
        AC.OP_drop_V * np.sin(AC.OP_drop_angle),
    ]
    box.Acc = [0, 0, (1+AC.OP_drop_Dn)*atm.g]

    # get trajectory of (each) box from inital state & CD(t)
    while box.Pos[-1] > 0: #until touchdown
        box.update(dt_sim, atm)
        box.Poslog.append(box.Pos)

    # store results

    # plot