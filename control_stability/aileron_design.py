import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects
from scipy.integrate import quad
from sympy import Eq, symbols, solve
sys.path.append('..')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

from horizontal_tailplane import hor_run
hor_run(aircraft)

def aileron_design(aircraft):
    ''' '''
    '''Roll performance parameters'''
    phi_des = 60 * (np.pi/180)  # Roll requirement
    t_lim = 1.3  # Max time available to complete the roll requirement
    t = 3  # starting time

    '''Design parameters'''
    delta_a_max = 20 * (np.pi/180) # max aileron deflection
    tau = 0.48  # factor based on the fraction of chord that is the aileron
    Ixx = 1200  # mass moment of inertia x-axis [kg m^2] TODO Update value 
    p = 0.94
    y0 = 0.95 * aircraft.b/2    #distance symmetry line aircraft to start aileron

    while t > t_lim:
        y1 = p * aircraft.b/2       #distance symmetry line aircraft to end aileron

        C_l_delta_a = (2*aircraft.CLa_w_cruise* tau *aircraft.rootchord/(aircraft.Sw * aircraft.b)) * ((y0**2/2 + 2/3 * ((aircraft.taper - 1)/aircraft.b) * y0**3)- (y1**2/2 + 2/3 * ((aircraft.taper - 1)/aircraft.b) * y1**3))
        
        C_l = C_l_delta_a * delta_a_max

        L_a = 1/2 * atm.rho0 * (aircraft.V_s_min * 1.3)**2 * aircraft.Sw * aircraft.b * C_l
        
        y_d = 0.4 * aircraft.b/2        # [m] This is the average distance between the centre of gravity of the aircraft and the centre of drag
        Cdr = 0.9
        Pss = np.sqrt(2*L_a/((aircraft.Sw + aircraft.AE_Sh_S * aircraft.Sw + aircraft.AE_Sv_S * aircraft.Sw) * Cdr * y_d**3))
    
        phi_1 = Ixx/(atm.rho0 * y_d**3* Cdr * (aircraft.Sw + aircraft.AE_Sh_S * aircraft.Sw + aircraft.AE_Sv_S * aircraft.Sw)) * np.log(Pss**2)

        P_roll_rate = Pss**2/(2*phi_1)

        t = np.sqrt(2*phi_des/P_roll_rate)

        p = p - 0.01


    print("\n--------Aileron Design------------\n")
    print(f"The aileron will span from {round(p*100, 3)}% to 95% of the span of the wing")
    print(f"The time to roll {phi_des * 180/np.pi} degrees will be {t} s")

    aircraft.y_a_0 = p
    aircraft.y_a_1 = 0.95

    delta_L = L_a/2 /((aircraft.y_a_1 + aircraft.y_a_0)/2 * aircraft.b/2)
    print(f"delta_L:{delta_L}")

def elevator_design(aircraft):
    
    # C_L_c = 2 * aircraft.W_TO *atm.g0/(aircraft.rho * aircraft.V_cruise**2 * aircraft.Sw)
    # C_L_TO = C_L_c + aircraft.CS_dCLmax_TO
    # C_D_TO = aircraft.CD0 + 1/(aircraft.e * aircraft.A * np.pi)*C_L_TO**2
    # D_TO = 1/2 * atm.rho0 * C_D_TO * aircraft.Vs_min**2 * aircraft.Sw 
    # L_TO = 1/2 * atm.rho0 * C_L_TO * aircraft.Vs_min**2 * aircraft.Sw 
    # M_ac_wf = 1/2*atm.rho0 * aircraft.Cm_ac_w * aircraft.Vs_min**2 * aircraft.Sw * aircraft.MAC_length

    CL_h = -1
    L_h = 0.5 * 0.7 * CL_h * (aircraft.V_cruise)**2 * aircraft.AE_Sh_S * aircraft.Sw
    print(L_h)

elevator_design(aircraft)


def rudder_design(aircraft):
    
    V_w = 40 * (0.51444444444) # [m/s] side wind

    V_t = np.sqrt((1.3 * aircraft.V_s_min)**2 + V_w**2) # [m/s] total speed

    S_s = 1.02 * (aircraft.l_f * aircraft.h_out + aircraft.AE_Sv_S * aircraft.Sw)

    x_ca = aircraft.l_f * aircraft.h_out * aircraft.X_FCG + aircraft.CS_S_v * aircraft.l_v
    d_ca = x_ca - aircraft.x_cg_position_aft

    C_d_y = 0.6 # Assumption

    F_w = 0.5 * atm.rho0 * C_d_y * V_w**2 * S_s

    angle_beta = np.arctan(V_w/(1.3 * aircraft.V_s_min))
    
    C_L_alpha_v = 2.95
    dsigma_dalpha = 0 
    eta_v = 0.96
    
    vertical_tail_volume = l_v*aircraft.CS_S_v/(aircraft.b*aircraft.Sw)

    C_n_beta = 0.75 * C_L_alpha_v * (1-dsigma_dalpha) * eta_v * vertical_tail_volume
    C_y_beta = -1.35 * C_L_alpha_v * (1 - dsigma_dalpha) * eta_v * aircraft.AE_Sv_S

    tau_r = 0.51
    C_Y_delta_r = C_L_alpha_v * eta_v * tau_r * 1 * aircraft.AE_Sv_S 
    C_n_delta_r = -C_L_alpha_v * vertical_tail_volume * eta_v * tau_r * 1



    # Define the variables
    delta_r, sigma = symbols('delta_r sigma')

    # Define the equations
    equation1 = Eq(0.5 * atm.rho0 * V_t**2 * aircraft.Sw * aircraft.b * (C_n_beta*(angle_beta - sigma) + C_n_delta_r * delta_r) + F_w * np.cos(sigma) * d_ca, 0)
    equation2 = Eq(0.5 * atm.rho0 * V_w**2 * S_s * C_d_y - 0.5*atm.rho0 * V_t**2 * aircraft.Sw * (C_y_beta*(angle_beta-sigma)+C_Y_delta_r*delta_r))

    # Solve the system of equations
    solution = solve((equation1, equation2), (delta_r, sigma))

    # Extract the values
    delta_r_value = solution[delta_r]
    sigma_value = solution[sigma]

    print(f"delta_r_required:{delta_r_value}")
    print(f"sigma_required:{sigma_value}")

def main_control_surface(aircraft):
    aileron_design(aircraft)
    rudder_design(aircraft)
    elevator_design(aircraft)