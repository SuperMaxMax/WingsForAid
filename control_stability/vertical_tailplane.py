
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('..')

from parameters import atmosphere

atm      = atmosphere()

def lateral_coefficients(aircraft, lcg):

    """Torenbeek Chapter 9, page 337 and 338"""
    K_beta = 0.3 * lcg/aircraft.l_f + 0.75 * aircraft.h_out/aircraft.l_f - 0.105     # (9-65)

    C_n_beta_fus = -K_beta * (aircraft.l_f*aircraft.h_out) * aircraft.l_f / (aircraft.Sw * aircraft.b) * (aircraft.h_out/aircraft.h_out)**(0.5) * (aircraft.w_out/aircraft.w_out)**(1/3)   #(9-64)

    C_n_beta_prop = -0.053 * aircraft.CS_n_blades * (aircraft.X_LEMAC) * (2 * aircraft.prop_radius)**2 / (aircraft.Sw * aircraft.b)

    C_n_beta_i = -0.017 #(9.67, high wing)

    C_n_beta_spec = 0.065  # This is the C_n_beta from the Cessna-152 and as we have to be at least as stable as them this is our requirements
    C_Y_v_alpha = 3 # This is an assumed value in [1/rad]. This depends on the chosen airfoil for the vertical tailplane. 
    # l_v = aircraft.l_h # It was assumed the horizontal tailplane was positioned at the same position as the vertical tailplane.
    aircraft.y_mac_v = 1/3*(aircraft.AE_b_v/2)*(1+2*0.7)/(1+0.7)
    l_v = (aircraft.l_f - aircraft.l_fus_tail_cone + aircraft.l_f_boom - 1/2 * aircraft.AE_rootchord_h) + (1/4 * aircraft.AE_rootchord_v + aircraft.y_mac_v * np.cos(aircraft.sweep_co4_v)) - (aircraft.X_LEMAC+ aircraft.X_cg_aft*aircraft.MAC_length)
    aircraft.AE_l_v = l_v

    Sv_S = aircraft.b/l_v * (C_n_beta_spec - (C_n_beta_fus + C_n_beta_prop + C_n_beta_i)) / (C_Y_v_alpha * (aircraft.AE_Vv_V)**2)
    aircraft.Sv_S = Sv_S
    aircraft.AE_Sv = Sv_S * aircraft.Sw



def ver_run(aircraft):
    lcg = aircraft.X_LEMAC + aircraft.X_cg_aft * aircraft.MAC_length
    lateral_coefficients(aircraft, lcg)



