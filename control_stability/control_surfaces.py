import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects
from scipy.integrate import quad
from scipy.optimize import fsolve
from math import *

sys.path.append('..')

from parameters import atmosphere
atm      = atmosphere()

def aileron_design(aircraft):  # Using Aircraft design: A Systems Engineering Approach, Ch 12.4
    ''' '''
    '''Roll performance requirments'''
    phi_des = 60 * (pi/180)  # Roll requirement
    t_lim = 1.3  # Max time available to complete the roll requirement
    t = 3  # starting time

    '''Design parameters'''
    delta_a_max = 20 * (pi/180) # max aileron deflection
    tau = 0.48  # factor based on the fraction of chord that is the aileron
    Ixx = 1200  # mass moment of inertia x-axis [kg m^2] TODO Update value 
    ystart_a = 0.94 * aircraft.b/2 # staring location of aileron
    yend_a = 0.95 * aircraft.b/2  # end location of aileron
    y_d = 0.4 * aircraft.b/2  # average distance between the centre of gravity of the aircraft and the centre of drag [m] 
    Cdr = 0.9  # aircraft drag coefficient in rolling motion

    while t > t_lim:
        # v Roll moment calculation v
        C_l_delta_a = (2*aircraft.CLa_w_cruise* tau *aircraft.rootchord/(aircraft.Sw * aircraft.b)) * ((yend_a**2/2 + 2/3 * ((aircraft.taper - 1)/aircraft.b) * yend_a**3)- (ystart_a**2/2 + 2/3 * ((aircraft.taper - 1)/aircraft.b) * ystart_a**3))
        C_l_max = C_l_delta_a * delta_a_max
        L_a = 1/2 * atm.rho0 * (aircraft.V_s_min * 1.3)**2 * aircraft.Sw * aircraft.b * C_l_max  # Rolling moment generated by ailerons [Nm]
        
        # v Steady State Rollrate v
        Pss = sqrt(2*L_a/((aircraft.Sw + aircraft.AE_Sh_S * aircraft.Sw + aircraft.AE_Sv_S * aircraft.Sw) * Cdr * y_d**3))
    
        phi_1 = Ixx/(atm.rho0 * y_d**3* Cdr * (aircraft.Sw + aircraft.AE_Sh_S * aircraft.Sw + aircraft.AE_Sv_S * aircraft.Sw)) * np.log(Pss**2)

        P_roll_rate = Pss**2/(2*phi_1)

        t = sqrt(2*phi_des/P_roll_rate)

        ystart_a = ystart_a - 0.01*aircraft.b/2


    print("\n--------Aileron Design------------\n")
    print(f"The aileron will span from {round(ystart_a*100, 3)}% to 95% of the span of the wing")
    print(f"The time to roll {phi_des * 180/np.pi} degrees will be {t} s")

    aircraft.y_a_0 = ystart_a
    aircraft.y_a_1 = 0.95

    delta_L = L_a/2 /((aircraft.y_a_1 + aircraft.y_a_0)/2 * aircraft.b/2)

def elevator_design(aircraft):
    
    # C_L_c = 2 * aircraft.W_TO *atm.g0/(aircraft.rho * aircraft.V_cruise**2 * aircraft.Sw)
    # C_L_TO = C_L_c + aircraft.CS_dCLmax_TO
    # C_D_TO = aircraft.CD0 + 1/(aircraft.e * aircraft.A * np.pi)*C_L_TO**2
    # D_TO = 1/2 * atm.rho0 * C_D_TO * aircraft.Vs_min**2 * aircraft.Sw 
    # L_TO = 1/2 * atm.rho0 * C_L_TO * aircraft.Vs_min**2 * aircraft.Sw 
    # M_ac_wf = 1/2*atm.rho0 * aircraft.Cm_ac_w * aircraft.Vs_min**2 * aircraft.Sw * aircraft.MAC_length

    CL_h = -1
    L_h = 0.5 * 0.7 * CL_h * (aircraft.V_cruise)**2 * aircraft.AE_Sh_S * aircraft.Sw

def rudder_design(aircraft):
    
    V_w = 20 * (0.51444444444) # [m/s] side wind

    V_t = sqrt((1.3 * aircraft.V_s_min)**2 + V_w**2) # [m/s] total speed
    S_s = 1.02 * (aircraft.l_f * aircraft.h_out + aircraft.AE_Sv_S * aircraft.Sw)

    x_ca = aircraft.l_f * aircraft.h_out * aircraft.X_FCG + aircraft.AE_S_v * aircraft.AE_l_v
    d_ca = x_ca - aircraft.x_cg_position_aft

    C_d_y = 0.6 # Assumption

    F_w = 0.5 * atm.rho0 * C_d_y * V_w**2 * S_s

    angle_beta = np.arctan(V_w/(1.3 * aircraft.V_s_min))

    C_L_alpha_v = 3
    dsigma_dalpha = 0 
    eta_v = 0.96
    
    
    vertical_tail_volume = aircraft.AE_l_v*aircraft.AE_S_v/(aircraft.b*aircraft.Sw)
    print(f"vert_tail_vol:{vertical_tail_volume}")
    C_n_beta = 0.75 * C_L_alpha_v * (1 - dsigma_dalpha) * eta_v * vertical_tail_volume
    C_y_beta = -1.35 * C_L_alpha_v * (1 - dsigma_dalpha) * eta_v * aircraft.AE_Sv_S

    tau_r = 0.80
    C_Y_delta_r = C_L_alpha_v * eta_v * tau_r * 1 * aircraft.AE_Sv_S 
    C_n_delta_r = -C_L_alpha_v * vertical_tail_volume * eta_v * tau_r * 1

    # Define the equations
    def equations(x):
        delta_r, sigma = x
        equation1 = 0.5 * atm.rho0 * V_t**2 * aircraft.Sw * aircraft.b * (C_n_beta*(angle_beta - sigma) + C_n_delta_r * delta_r) + F_w * cos(sigma) * d_ca
        equation2 = 0.5 * atm.rho0 * V_w**2 * S_s * C_d_y - 0.5*atm.rho0 * V_t**2 * aircraft.Sw * (C_y_beta*(angle_beta-sigma)+C_Y_delta_r*delta_r)
        return [equation1, equation2]
    # Solve the system of equations
    initial_guess = [1, 1]
    solution = fsolve(equations, initial_guess)

    # Extract the values
    delta_r_value = solution[0]
    sigma_value = solution[1]
    
    if delta_r_value*180/pi > 30:
        print(f"\nRequired rudder deflection ({delta_r_value*180/pi} degrees) exceeds maximum rudder deflection (30 degrees)\n")
    else:
        print(f"delta_r_required:{delta_r_value*180/pi}")
        print(f"sigma:{sigma_value}")

def main_control_surface(aircraft):
    aileron_design(aircraft)
    rudder_design(aircraft)
    elevator_design(aircraft)