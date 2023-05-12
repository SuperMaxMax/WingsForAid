import subprocess
import os
import Class_I_weight_estimation as c1
import geometry_determination as geo
import parameters as para
import numpy as np

W_TO = 750

class initial_sizing:
    def __init__(self):
        pass

    def iteration(self):
        #self.W_TO = 750 # kg

        #weight = c1.Weight(para)
        #weight.iteration()
            # Loop over class 1 and geometry determination

            # Class 1 iteration

        os.system("python Class_I_weight_estimation.py")
            # Geometry determination

            # Class II estimation



if __name__ == "__main__":
    initsiz = initial_sizing()
    initsiz.iteration()

