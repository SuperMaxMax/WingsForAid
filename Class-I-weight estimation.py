import numpy as np
from math import *
import parameters as para


class Weight:
    def __init__(self, para):
        self.para = para
        self.W_TO = 750
    #########################################################################
    "CLASS I WEIGHT ESTIMATION"
    #########################################################################

    def weight_take_off(self):

        W_TO = self.W_OE + self.W_F + self.W_PL
        return(W_TO)


    def weight_empty_operational(self):

        self.W_E = 0.5482 * self.W_TO + 486.68
        self.W_tfo = 0.002 * self.W_TO
        self.W_OE = self.W_E + self.W_tfo

    def L_D_cruise(self):
        self.L_D = np.sqrt(np.pi * self.para.A * self.para.e / self.para.CD0)
        return(self.L_D)

    def W5W4(self):
        # Cruise fuel fraction
        weight.L_D_cruise()
        W5W4 = 1/exp(self.R * self.para.g * self.para.c_p / self.L_D)
        return(W5W4)

    def W7W5(self):
        # Dropping manoeuvre
        W7W5 = self.para.W4W3 * self.para.W10W9
        return(W7W5)

    def Mff(self):
        self.R = self.para.R/self.para.n_drops
        W7W5 = weight.W7W5()
        print(W7W5)
        W5W4 = weight.W5W4()
        print(W5W4)
        self.Mff = self.para.W1W_TO * self.para.W2W1 * self.para.W3W2 * self.para.W4W3 * W5W4 * W7W5 * W5W4 * self.para.W10W9 * self.para.WfinalW10

    def weight_fuel_used(self):
        weight.Mff()

        self.M_f_used = (1 - self.Mff) * self.W_TO

    def weight_fuel(self):
        weight.weight_fuel_used()
        self.W_F = (1 + self.para.M_res) * self.M_f_used


    def iteration(self):
        it = True
        while it == True:
            weight.weight_fuel()
            weight.weight_empty_operational()
            self.W_PL = self.para.W_PL
            W_TO_new = weight.weight_take_off()
            print(W_TO_new)
            change = (W_TO_new - self.W_TO)/self.W_TO
            if change < 0.01:
                it = False
            else:
                self.W_TO = W_TO_new
                print(self.W_TO)
        print("WTO",self.W_TO)

if __name__ == "__main__":
    weight = Weight(para)

    #weight.weight_wing()
    #weight.weight_empennage()
    #weight.Mff()
    weight.iteration()