
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('..')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

def lateral_coefficients(aircraft, lcg):

    """Torenbeek Chapter 9, page 337 and 338"""
    K_beta = 0.3 * lcg/aircraft.l_f + 0.75 * aircraft.h_out/aircraft.l_f - 0.105     # (9-65)

    C_n_beta_fus = -K_beta * (aircraft.l_f*aircraft.h_out) * aircraft.l_f / (aircraft.Sw * aircraft.b) * (aircraft.h_out/aircraft.h_out)**(0.5) * (aircraft.w_out/aircraft.w_out)**(1/3)   #(9-64)

    C_n_beta_prop = -0.053 * aircraft.CS_n_blades * (aircraft.X_LEMAC) * aircraft.CS_D_prop**2 / (aircraft.Sw * aircraft.b)

    C_n_beta_i = -0.017 #(9.67, high wing)

    C_n_beta_spec = 0.065  # This is the C_n_beta from the Cessna-152 and as we have to be at least as stable as them this is our requirements
    C_Y_v_alpha = 0.05*180/np.pi # This is an assumed value in [1/rad]. This depends on the chosen airfoil for the vertical tailplane. 
    l_v = aircraft.CS_l_h # It was assumed the horizontal tailplane was positioned at the same position as the vertical tailplane. 

    Sv_S = aircraft.b/l_v * (C_n_beta_spec - (C_n_beta_fus + C_n_beta_prop + C_n_beta_i)) / (C_Y_v_alpha * (aircraft.CS_Vv_V)**2)
    aircraft.CS_Sv_S = Sv_S


def ver_run(aircraft):
    lcg = aircraft.X_LEMAC + aircraft.X_cg_full * aircraft.MAC_length
    lateral_coefficients(aircraft, lcg)



