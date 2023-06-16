# ------- Import statements -------
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

# -------- Import packages --------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from scipy.integrate import quad
# ------- Import other files -------
from parameters import UAV, airport, atmosphere
from flight_performance.simulation import atm_parameters

# -------- Create objects --------
aircraft = UAV("Pel")
airfield = airport("Woensel")
atm      = atmosphere()

def Chordlength(y):
    return aircraft.rootchord * (1 - ((1-aircraft.taper) / (aircraft.b / 2)) * y)

def HLDforces(ac_obj, atm_obj, V_extend, h_extend, n_boxes, W_F):
    # ------ Control surface areas ------
    FlapArea    = quad(Chordlength*(1 - ac_obj.xc_aft_spar), 0, ac_obj.yend_flap)[0]
    # ---------- Deflections ------------
    MaxFlapDeflection = 60 * np.pi/180                          # [rad]
    # ------ Find aircraft weight -------
    W = ac_obj.W_OE + ac_obj.boxweight * n_boxes + W_F
    # --- Find atmospheric properties ---
    p, rho = atm_parameters(atm_obj, h_extend)[0], atm_parameters(atm_obj, h_extend)[2]
    # ------------ Find AOA -------------
    CL = 2*W*atm_obj.g / (rho * V_extend**2 * ac_obj.Sw)
    alpha_w = CL/ac_obj.CL_a_w + ac_obj.af_alpha0
    # -- Find lift with deflected flap --
    L_ext   = 1/2 * rho * V_extend**2 * ac_obj.Sw * ac_obj.FP_CL_max_land
    F_flaps = L_ext - W
    c_p_flap = 0.5
    c_flap_hinge = 0.1
    

def Aileronforces(ac_obj, atm_obj):
    c_ail_hinge = 0.1                                                           # hinge location as fraction of chord
    cp_ail      = 0.5                                                           # cp location as fraction of chord
    F_ail = 1/2 * atm.rho0 * (ac_obj.V_A)**2 * aircraft.Sw * ac_obj.Clmax_ail   # Force experienced by a single aileron
    M_h_ail = (cp_ail - c_ail_hinge) * Chordlength((ac_obj.yend_ail + ac_obj.ystart_ail)/2) * F_ail
    return M_h_ail

print(Aileronforces(aircraft, atm))