import numpy as np
from math import *
from parameters import *
import matplotlib.pyplot as plt

#########################################################################
"CLASS I WEIGHT ESTIMATION"
#########################################################################


def weight_take_off(W_OE, W_F, W_PL):
    # take_off weight is a function of operational empty weight, fuel weight and payload weight
    W_TO = W_OE + W_F + W_PL
    return W_TO


def weight_empty_operational(W_TO):
    # Operational empty weigth is a function of the empty weight and the trapped fuel weight, which are both a function of take_off weight
    W_E = lin_par1 * W_TO + lin_par2  # relation from Roskam
    W_tfo = 0.002 * W_TO
    W_OE = W_E + W_tfo
    return W_OE


def L_D_cruise():
    # Lift over drag calculation
    L_D = np.sqrt(np.pi * A * e / (4 * CD0))
    return L_D


def W5W4_calculation():
    # Cruise fuel fraction calculation
    L_D = L_D_cruise() # Lift over drag ratio from method above
    W5W4 = 1/exp(R * g * c_p / (prop_eff * L_D)) #Breguet's range equation
    return W5W4


def W7W5_calculation():
    # Dropping manoeuvre
    W7W5 = W4W3 * W10W9
    return W7W5


def Mff_calculation():
    W7W5 = W7W5_calculation()
    W5W4 = W5W4_calculation()

    Mff = W1W_TO * W2W1 * W3W2 * W4W3 * W5W4 * W7W5**n_drops * W10W9 * WfinalW10
    return Mff


def weight_fuel_used(W_TO):
    Mff = Mff_calculation()

    M_f_used = (1 - Mff) * W_TO
    return M_f_used


def weight_fuel(W_TO):
    M_f_used = weight_fuel_used(W_TO)

    W_F = (1 + M_res) * M_f_used
    return W_F


def iteration(W_TO):
    it = True
    while it:
        W_F = weight_fuel(W_TO)                        # perform calculation of fuel weight
        W_OE = weight_empty_operational(W_TO)           # perform calculation of operational empty weight
        W_TO_new = weight_take_off(W_OE, W_F, W_PL)         # combines weights to find total weight
        change = (W_TO_new - W_TO)/W_TO
        if abs(change) < 0.001:
            W_TO = W_TO_new    # change between iteration is smaller than 0.1 percent
            it = False
        else:
            W_TO = W_TO_new

    print(f"W_OE:{W_OE}")
    print(f"W_F:{W_F}")
    print(f"W_PL:{W_PL}")
    print(f"W_TO:{W_TO}")

    return W_OE, W_F, W_PL, W_TO


if __name__ == "__main__":
    W_TO = 600
    iteration(W_TO)
