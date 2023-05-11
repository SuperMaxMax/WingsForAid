import numpy as np
from math import *
import parameters as para


class Weight:
    def __init__(self, para):
        # parameters.py is a file containing all constant parameters of the aircraft
        self.para = para
        # take_off weight that is assumed initially
        self.W_TO = 750

    #########################################################################
    "CLASS II WEIGHT ESTIMATION"
    #########################################################################

    def weight_take_off(self):
        # take_off weight is a function of operational empty weight, fuel weight and payload weight
        W_TO = self.W_OE + self.W_F + self.W_PL
        return W_TO

    def weight_empty_operational(self):
        # Operational empty weigth is a function of the empty weight and the trapped fuel weight, which are both a function of take_off weight
        W_E = para.lin_par1 * self.W_TO + para.lin_par2  # relation from Roskam
        W_tfo = 0.002 * self.W_TO
        self.W_OE = W_E + W_tfo

    def L_D_cruise(self):
        # Lift over drag calculation
        self.L_D = np.sqrt(np.pi * self.para.A * self.para.e / (4 * self.para.CD0))

    def W5W4(self):
        # Cruise fuel fraction calculation
        weight.L_D_cruise() # Lift over drag ratio from method above
        W5W4 = 1/exp(self.R * self.para.g * self.para.c_p / (self.para.prop_eff * self.L_D))
        return(W5W4)

    def W7W5(self):
        # Dropping manoeuvre
        W7W5 = self.para.W4W3 * self.para.W10W9
        return(W7W5)

    def Mff_calculation(self):
        self.R = self.para.R/self.para.n_drops  # flight range between drops
        W7W5 = weight.W7W5()
        W5W4 = weight.W5W4()

        self.Mff = self.para.W1W_TO * self.para.W2W1 * self.para.W3W2 * self.para.W4W3 * W5W4 * (W7W5 * W5W4)**self.para.n_drops * self.para.W10W9 * self.para.WfinalW10

    def weight_fuel_used(self):
        weight.Mff_calculation()

        self.M_f_used = (1 - self.Mff) * self.W_TO

    def weight_fuel(self):
        weight.weight_fuel_used()

        self.W_F = (1 + self.para.M_res) * self.M_f_used

    def iteration(self):
        it = True

        while it:
            weight.weight_fuel()
            weight.weight_empty_operational()
            self.W_PL = self.para.W_PL
            W_TO_new = weight.weight_take_off()
            change = (W_TO_new - self.W_TO)/self.W_TO
            if change < 0.01:
                it = False
            else:
                self.W_TO = W_TO_new
                print(self.W_TO)

        print(f"W_OE:{self.W_OE}")
        print(f"W_F:{self.W_F}")
        print(f"W_PL:{self.W_PL}")
        print(f"W_TO:{self.W_TO}")

        return self.W_OE, self.W_F, self.W_PL, self.W_TO


if __name__ == "__main__":
    weight = Weight(para)
    weight.iteration()