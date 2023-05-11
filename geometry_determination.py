import  numpy as np
from    parameters import *

# W/S and W/P diagrams
# drag polar
class Geometry:
    def __init__(self):
        pass

    def dragpolar(self, CL, CD0, e, A):
        CD = CD0 + (CL**2/(np.pi*A*e))
        return CD

    def stallWS(self, V, rho, CL_max):
        WoS = 1/2 * rho * V**2 * CL_max
        return WoS

    def TOP(self, WS, sigma, CL_TO, BHP, W_TO):
        TOP = WS/(sigma*CL_TO*(BHP/W_TO))
        return TOP

    def altitude_effects(self, h, Lambda, R, g, T0, rho0, BHP0):
        rho     = rho0*(1+ ((Lambda*h)/T0))**(-((g/(R*Lambda))+1))
        sigma   = rho/rho0
        BHP     = BHP0*(sigma)**(3/4)
        return rho, sigma, BHP