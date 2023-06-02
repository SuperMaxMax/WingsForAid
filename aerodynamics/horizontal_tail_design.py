import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate

from parameters import UAV
aircraft = UAV('aircraft')

def required_lift():
    #load in parameters from current design
    AR = aircraft.A
    sweep_c4 = aircraft.lambda_co4
    taper = aircraft.taper
    C_m_af = aircraft.AE_cm0
    alpha_t = aircraft.wing_twist
    V_c = aircraft.V_cruise
    rho_c = aircraft.rho_cruise
    W_TO = aircraft.W_TO * 9.80665
    S = aircraft.Sw
    h0 = aircraft.MAC_ac
    MAC = aircraft.MAC_length
    x_LEMAC = aircraft.X_LEMAC
    xcg_fwrd = aircraft.X_cg_fwd
    xcg_aft = aircraft.X_cg_aft

    #test to verify code
    AR = 28
    sweep_c4 = 0
    taper = 0.8
    C_m_af = -0.013
    alpha_t = -1.1 * np.pi/180
    V_c = 48.872222
    rho_c = 0.905
    W_TO = 850 * 9.81
    S = 18
    h0 = 0.23
    MAC = 0.8
    x_LEMAC = 1.61416
    xcg_fwrd = 1.86816
    xcg_aft = 1.86816

    #temp values from book
    V_H = 0.6 #table 6.4 "aircraft design synthesis a systems engineering approach"
    #V_H = Sh*lh / Sw * MAC

    #inbetween calculutions 
    C_L = 2*W_TO/(rho_c * (V_c**2) * S) #lift in cruise
    sweep_LE = np.arctan(np.tan(sweep_c4) - (4/AR) * (25/100 * (1-taper)/(1+taper))) #leading edge sweep
    sweep_LE = 8 * np.pi/180 #test to verify code
    C_m_0_wf = C_m_af * (AR * (np.cos(sweep_LE) ** 2)) / (AR + 2 * np.cos(sweep_LE)) + 0.01 * alpha_t
    x_ac_nose = x_LEMAC + h0*MAC
    most_extreme_cg = [xcg_fwrd, xcg_aft]


    C_L_h = 0 
    for xcg_wf in most_extreme_cg:
        x_cg_wf_nose = xcg_wf*MAC + x_LEMAC
        x_cg_wf_mac_difference = x_cg_wf_nose - x_ac_nose
        x_cg = (h0 * MAC + x_cg_wf_mac_difference)
        h = x_cg/MAC
   
        C_L_h_new = (C_m_0_wf + C_L * (h - h0)) / (V_H)
        if abs(C_L_h_new) > abs(C_L_h):
            C_L_h = C_L_h_new
    
    print(C_L_h)
    


required_lift()