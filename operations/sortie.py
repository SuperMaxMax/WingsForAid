import sys

import matplotlib.collections

sys.path.append("..")

# Start your import below this

import numpy as np
from parameters import UAV, atmosphere
import pandas as pd
import matplotlib.pyplot as plt
import flight_performance.simulation as fpsim

def operations_eval(ac):
    print("\noperations evaluation")
    atm = atmosphere()
    TTFD_s, T_s, Wf = fpsim.fuelusesortie(ac, atm, ac.n_boxes, ac.OP_n_drops, ac.h_cruise, ac.W_F, ac.V_cruise, ac.OP_Range)

    # print("TTFD(sortie) [fh]", np.round(TTFD_s/3600,3))
    # print("T(sortie) [fh]", np.round(T_s/3600,3))
    # print("Wf [kg]", np.round(Wf, 2))
    # print("\n")

    # REQ-check: time to first delivery
    TTFD = ac.OP_TTFS + ac.OP_T_ground + TTFD_s/3600 # sum time until first drop (5th time entry)
    print(str("TTFD: "+str(np.round(TTFD, 1))+"[hr]"))

    # REQ-check: sortie per day
    T_sortie = ac.OP_T_ground + ac.OP_T_pilot + T_s/3600 # total per AC
    N_sortie_per_day = 24/(T_sortie)
    print(str("Ns/day: "+str(np.round(N_sortie_per_day,2))+"[#]"))

    # REQ-check: payload rate
    # enforce PL-rate:
    PL_per_AC = ac.OP_N_boxes_per_sortie * ac.OP_PL_per_box
    PL_per_ACday = PL_per_AC * N_sortie_per_day
    ac.OP_AC_per_op = np.round(ac.OP_MR_PL / PL_per_ACday, 0)
    print(str("AC required: " + str(ac.OP_AC_per_op) + "[#AC]"))
    # eval PL-rate
    PL_per_day = PL_per_AC * (ac.OP_AC_per_op * N_sortie_per_day)
    print(str("PL per day: "+str(np.round(PL_per_day,1))+"[kg/24h]"))

    # REQ-chack: cost per kg
    N_sortie_tot = (N_sortie_per_day * AC.OP_T_ops) * ac.OP_AC_per_op * ac.OP_N_ops
    # payload total
    PL_tot =  N_sortie_tot * PL_per_AC
    # fuel volume per sortie
    M_Fs = Wf # (t*P)*SFC [Ws*kg/J] = kg fuel required
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

    ac.CST_PL_per_kg = CST_per_kg
    ac.N_sortie_per_day = N_sortie_per_day
    ac.OP_TTFD = TTFD
    ac.T_sortie = T_sortie
    return CST_per_kg, N_sortie_per_day, TTFD

if __name__ == '__main__':
    AC = UAV('aircraft')
    print("AC default:", AC.__dict__)
    operations_eval(AC)