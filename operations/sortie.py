import sys

import matplotlib.collections

sys.path.append("..")

# Start your import below this

import numpy as np
from parameters import UAV
import pandas as pd
import matplotlib.pyplot as plt
#import flight_performance.simulation.py as fpsim

def flight_profile(ac): # create flight profile based on REQ, design and CHOICE (opti)
    hp_to_watt = 745.699872
    FP = {
        "time [s]": [],
        "power [W]": [],
        "altitude [m]": [],
        "distance [m]": [],
    }

    # 0:taxi
    FP["time [s]"].append(10 * 60)
    FP["power [W]"].append(0.2 * ac.power * hp_to_watt)
    FP["altitude [m]"].append(ac.h_TO)
    FP["distance [m]"].append(0)

    # 1:TO
    FP["time [s]"].append(10 * 60)
    FP["power [W]"].append(1 * ac.power * hp_to_watt)
    FP["altitude [m]"].append(ac.h_TO)
    FP["distance [m]"].append(ac.RFL)

    # per drop
    N_drops = ac.OP_n_drops
    R_intercruise = 1000*ac.OP_Range/(N_drops-1)
    while N_drops > 0:
        # 2:cruise
        FP["time [s]"].append(R_intercruise / ac.V_cruise)
        FP["power [W]"].append(ac.BHP_cruise * hp_to_watt)
        FP["altitude [m]"].append(ac.h_cruise)
        FP["distance [m]"].append(R_intercruise)

        # 3:descent
        FP["time [s]"].append(45 * 60)
        FP["power [W]"].append(0.1 * ac.power * hp_to_watt)
        FP["altitude [m]"].append(ac.OP_h_loiter)
        FP["distance [m]"].append(1000)

        # 4:loiter
        FP["time [s]"].append(20 * 60)
        FP["power [W]"].append(0.3 * ac.power * hp_to_watt)
        FP["altitude [m]"].append(ac.OP_h_loiter)
        FP["distance [m]"].append(0)

        # 5:drop
        FP["time [s]"].append(10 * 60)
        FP["power [W]"].append(0.8 * ac.power * hp_to_watt)
        FP["altitude [m]"].append(ac.h_TO)
        FP["distance [m]"].append(500)
        print(len(FP["time [s]"]))

        # 6:climb
        FP["time [s]"].append(45 * 60)
        FP["power [W]"].append(0.9 * ac.power * hp_to_watt)
        FP["altitude [m]"].append(1000*ac.h_cruise)
        FP["distance [m]"].append(3000)

        N_drops -= 1

    # cruise (back)
    FP["time [s]"].append(ac.OP_Range / ac.V_cruise)
    FP["power [W]"].append(ac.BHP_cruise * hp_to_watt)
    FP["altitude [m]"].append(ac.h_cruise)
    FP["distance [m]"].append(ac.OP_Range)

    # descent
    FP["time [s]"].append(30 * 60)
    FP["power [W]"].append(0.1 * ac.power * hp_to_watt)
    FP["altitude [m]"].append(ac.OP_h_loiter)
    FP["distance [m]"].append(1000)

    # LDG
    FP["time [s]"].append(15 * 60)
    FP["power [W]"].append(0.2 * ac.power * hp_to_watt)
    FP["altitude [m]"].append(ac.h_TO)
    FP["distance [m]"].append(ac.LDG_dist)

    return FP

def operations_eval(ac):
    Flight = flight_profile(ac)
    for boolean, tpm in Flight.items():
        print(str(boolean+":"+str(tpm)))
        tpm = np.array(tpm)

    # REQ-check: time to first delivery
    TTFD = ac.OP_TTFS + ac.OP_T_sortie_gnd + np.sum(Flight["time [s]"][0:5])/3600 # sum time until first drop (5th time entry)
    print(str("TTFD: "+str(np.round(TTFD, 1))+"[hr]"))

    # REQ-check: sortie per day
    T_sortie = ac.OP_T_sortie_gnd + np.sum(Flight["time [s]"]) # total per AC
    N_sortie_per_day = 24/(T_sortie/3600)
    print(str("Ns/day: "+str(np.round(N_sortie_per_day,2))+"[#]"))

    # REQ-check: payload rate
    # enforce PL-rate:
    PL_per_AC = ac.OP_N_boxes_per_sortie * ac.OP_PL_per_box
    PL_per_ACday = PL_per_AC * N_sortie_per_day
    ac.OP_AC_per_op = np.round(ac.OP_MR_PL / PL_per_ACday, 1)
    print(str("AC required: " + str(np.round(ac.OP_AC_per_op, 1)) + "[#AC]"))
    # eval PL-rate
    PL_per_day = PL_per_AC * (ac.OP_AC_per_op * N_sortie_per_day)
    print(str("PL per day: "+str(np.round(PL_per_day,1))+"[kg/24h]"))

    # REQ-chack: cost per kg
    N_sortie_tot = (N_sortie_per_day * AC.OP_T_ops) * ac.OP_AC_per_op * ac.OP_N_ops
    # payload total
    PL_tot =  N_sortie_tot * PL_per_AC
    # fuel volume per sortie
    M_Fs = np.vdot(Flight["time [s]"], Flight["power [W]"]) * ac.SFC # (t*P)*SFC [Ws*kg/J] = kg fuel required
    Vol_Fs = M_Fs / (AC.fueldensity)# kg/(kg/L) = L
    print("fuel per sortie:"+str(np.round(Vol_Fs,1))+"[L]")
    # fuel cost
    CST_Fs = Vol_Fs * AC.OP_fuelprice
    CST_Ftot = CST_Fs * N_sortie_tot
    print("total fuel cost:" + str(np.round(CST_Ftot, 1)) + "[euro]")
    #euro/kg (over tot operations, incl. fuel)
    CST_TOT = AC.OP_CST_nofuel + CST_Ftot
    CST_per_kg = CST_TOT / PL_tot
    print("cost per kg:"+str(np.round(CST_per_kg,3))+"[euro/kg]")

if __name__ == '__main__':
    AC = UAV('aircraft')
    print("AC default:", AC.__dict__)
    operations_eval(AC)