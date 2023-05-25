import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

# Start your import below this
# import packages
import numpy as np
import scipy as sp
import pandas as pd
# Import other files
from parameters import UAV, airport, atmosphere

aircraft = UAV("aircraft")
airfield = airport("Sudan")
atm      = atmosphere()

def takeoffweight(obj, W_F):
    ZFW = obj.W_OE + obj.n_boxes*obj.boxweight
    TOW = ZFW + W_F
    return TOW

def atm_parameters(obj, h):
    T    = obj.T0 + obj.lambd * h
    rho  = obj.rho0*np.power((1+(T/obj.T0)), (-((obj.g / (obj.lambd * obj.R))+1)))
    p    = obj.p0*np.power((1+(T/obj.T0)), (-(obj.g / (obj.lambd * obj.R))))
    a    = np.sqrt(obj.gamma*obj.R*T)
    return p, T, rho, a

def propthrust(obj, h, ):
    p, T, rho, a = atm_parameters(atm, h)

