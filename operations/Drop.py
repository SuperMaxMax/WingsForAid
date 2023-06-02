import numpy as np

class BOX:
    def __init__(self, box_number):
        # parameters
        self.CD = 0
        self.CD0 = 	2.1 #[-]
        self.DCD_brake = 0.5 #[-]
        self.DCD_flaps = 1.28 #[-]
        self.DT_brake = 0.5 # [s]
        self.DT_flaps = 1  # [s]
        self.T_deployment = 0.5 # [s] ease drag transition

        self.S = 0.4*0.4 # [m2]

        self.mass = 3 # [kg]

    def compute_CD(self, t):
        self.CD = self.CD0
        self.S = 0.4 ** 2
        if t > self.DT_brake:
            if t < self.DT_brake + self.T_deployment:
                ease = (t - self.DT_brake) / self.T_deployment
            else:
                ease = 1
            self.CD += self.DCD_brake * ease

        if t > self.DT_flaps:
            if t < self.DT_flaps + self.T_deployment:
                ease = (t - self.DT_flaps) / self.T_deployment
            else:
                ease = 1
            self.CD += self.DCD_flaps * ease
            self.S += 4 * (0.6 * 0.4) * ease