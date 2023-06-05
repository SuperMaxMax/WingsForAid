import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('.')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

#################################################################################################################
"Determine Atmosphereic conditions"
#################################################################################################################

def Mach_calculation(atm, V, h):
    """Calculation of the Mach number based on the velocity and altitude"""
    atm.T =  atm.T0 - atm.lambd * h
    atm.speed_of_sound = (atm.gamma * atm.R * atm.T) ** 0.5
    Mach = V / atm.speed_of_sound
    return Mach

#################################################################################################################
"Determine Scissorplot Coefficients"
#################################################################################################################

def LiftRateCoefficient(aircraft, Mach, A, lambda_co2):  # lift rate coefficient tail/wing
    """Lift rate coefficient of a wing(/tail)
    Equation from SEAD Lecture 7, slide 41"""
    aircraft.CS_beta = np.sqrt(1 - Mach ** 2)
    CLa = 2 * np.pi * A / (2 + np.sqrt(4 + ((A * aircraft.CS_beta / aircraft.CS_eta)** 2) * (1 + np.tan(lambda_co2) ** 2  / aircraft.CS_beta ** 2)))
    aircraft.CLa = CLa
    return CLa

def TaillessLiftRateCoefficient(aircraft, CLa): 
    """Lift rate coefficient of aircraft without tail
    Equation from SEAD Lecture 7, slide 42"""
    CLa_Ah = aircraft.AE_cl_alpha * (1 + 2.15 * aircraft.w_out / aircraft.b) * aircraft.Sw / aircraft.Sw + (np.pi * aircraft.w_out ** 2) / (2 * aircraft.Sw)

    return CLa_Ah


def nacelle_influence(aircraft, CLa_Ah):
    """a.c. shift due to nacelles
    Equation from SEAD Lecture 7, slide 38"""
    #dx_ac_n = (-4) * (aircraft.w_out ** 2 * (aircraft.X_LEMAC + 0.25 * aircraft.MAC_length)) / (aircraft.Sw * aircraft.MAC_length * CLa_Ah)   
    l_p = aircraft.X_LEMAC + 0.25 * aircraft.MAC_length
    dx_ac_n = -0.05 * aircraft.CS_n_blades * aircraft.CS_D_prop**2 * l_p / (aircraft.Sw * aircraft.MAC_length * CLa_Ah)
    return dx_ac_n


def WingFuselageAC(aircraft, CLa_Ah):
    """Aerodynamic center of wingfuselage combination, x_ac_wf [ ] (already divided by MAC)
    Equation from SEAD Lecture 7, slide 37"""
    x_ac_wf = aircraft.CS_x_ac_w - (1.8 / CLa_Ah) * (aircraft.w_out *aircraft.h_out * (aircraft.X_LEMAC - aircraft.x_lemac)) / (aircraft.Sw * aircraft.MAC_length) + \
    (0.273 / (1 + aircraft.taper)) * ((aircraft.w_out * (aircraft.Sw / aircraft.b) * (aircraft.b - aircraft.w_out)) / ((aircraft.MAC_length ** 2) * (aircraft.b + 2.15 * aircraft.w_out))) * np.tan(aircraft.lambda_co4)
    return x_ac_wf

def Aerodynamic_centre_determination(aircraft):
    aircraft.Mcruise            = Mach_calculation(atm, aircraft.V_cruise, aircraft.h_cruise)
    aircraft.Mmin               = Mach_calculation(atm, aircraft.V_s_min, aircraft.h_TO)
    
    aircraft.CLa_h_cruise       = LiftRateCoefficient(aircraft, aircraft.Mcruise, aircraft.AE_A_h, aircraft.AE_lambda_co2_h)

    aircraft.CLa_w_cruise       = LiftRateCoefficient(aircraft, aircraft.Mcruise, aircraft.A, aircraft.lambda_co2)
    aircraft.CLa_w_approach     = LiftRateCoefficient(aircraft, aircraft.Mmin, aircraft.A, aircraft.lambda_co2)

    aircraft.CLa_Ah_cruise      = TaillessLiftRateCoefficient(aircraft, aircraft.CLa_w_cruise)
    aircraft.CLa_Ah_approach    = TaillessLiftRateCoefficient(aircraft, aircraft.CLa_w_approach)

    aircraft.dx_ac_n_cruise     = nacelle_influence(aircraft, aircraft.CLa_Ah_cruise)
    aircraft.dx_ac_n_approach   = nacelle_influence(aircraft, aircraft.CLa_Ah_approach)

    aircraft.x_ac_wf_cruise     = WingFuselageAC(aircraft, aircraft.CLa_Ah_cruise)
    aircraft.x_ac_wf_approach   = WingFuselageAC(aircraft, aircraft.CLa_Ah_approach)

    aircraft.x_ac_cruise        = aircraft.x_ac_wf_cruise + aircraft.dx_ac_n_cruise
    aircraft.x_ac_approach      = aircraft.x_ac_wf_approach + aircraft.dx_ac_n_approach

#################################################################################################################
"FIXME: Determine Control Surface Coefficients"
#################################################################################################################
# Using data from Torenbeek aroung page 533
flaptype = 'singleslotted'              # Can be 'singleslotted' or 'fowler'

if flaptype == 'singleslotted':
    CS_deltaf_TO                   = 20  # flap deflection angle at take-off [deg] (ADSEE-II, L3 S14)
    CS_deltaf_LD                   = 40  # flap deflection angle at landing [deg]
    CS_fc_c                        = 0.25  # flap chord length / wing chord length [-] NOTE: This influences aft spar position TODO: update value
    CS_dc_cf_TO                    = 0.2  # increase in chord length / flap chord length [-]
    CS_dc_cf_LD                    = 0.3  # increase in chord length / flap chord length [-]
    CS_c_prime                     = 1  # wing chord length when flaps are extended [m] TODO: update value
    CS_dClmax                      = 1.3  # Additional airfoil lift due to single slotted flap (ADSEE-II, L3 S36)

if flaptype == 'fowler':
    CS_deltaf_TO                   = 15  # flap deflection angle at take-off [deg] (ADSEE-II, L3 S15)
    CS_deltaf_LD                   = 40  # flap deflection angle at landing [deg]
    CS_fc_c                        = 0.35  # flap chord length / wing chord length [-] NOTE: This influences aft spar position TODO: update value
    CS_dc_cf_TO                    = 0.48  # increase in chord length / flap chord length [-]
    CS_dc_cf_LD                    = 0.625  # increase in chord length / flap chord length [-]
    CS_c_prime                     = 1  # wing chord length when flaps are extended [m] TODO: update value
    CS_dClmax                      = 1.3 # Additional airfoil lift due to fowler flap (ADSEE-II, L3 S36)

CS_Swf                          = 0.5 * aircraft.Sw # spanwise portion of wing influenced by flaps (ADSEE-II, L3 S31) TODO: update value
CS_lambda_hinge                 = 0.02 # hinge line sweep angle, likely parallel to aft spar [rad] TODO: update value
CS_dCLmax = 0.9 * CS_dClmax * (CS_Swf/aircraft.Sw) * np.cos(CS_lambda_hinge)  # Increase in CL due to flap extension
# NOTE: Above equation can be rewritten to calculate Swf for a given delta CLmax 

#################################################################################################################
"Controlability and Stability Curves"
################################################################################################################

def controlability_curve(aircraft, xcgRange): #TODO: change constants here 
    """For the controllability curve, the CL_h, CL_Ah and Cm_ac of the aircraft needs to be found. 
    For CL_h: SEAD L8 S17
    For CL_Ah: L = W @ approach
    For Cm_ac: SEAD L8 S19"""

    CL_h = -1
    CL_Ah = aircraft.W_TO * aircraft.g0 / (0.5 * aircraft.rho0 * (1.3 * aircraft.V_s_min)**2 * aircraft.Sw)

    #Calculation of Cm_ac starting here
    Cm_ac_w = aircraft.AE_cm0 * (aircraft.A * (np.cos(aircraft.lambda_co4))**2 / (aircraft.A + 2 * np.cos(aircraft.lambda_co4)))

    # FIXME: Everything beyond this point is not yet checked, therefore errors will be present
    dCm_ac_f = aircraft.CS_mu2 * ((-aircraft.CS_mu1) * CS_dClmax * 1.06426 - (CL_Ah + CS_dClmax * (1 - 42.47695 / aircraft.Sw)) * (1/8) * 1.06426 * (1.06426 - 1)) + 0.7 * (aircraft.A / (1 + 2 / aircraft.A)) * aircraft.CS_mu3 * CS_dClmax * np.tan(aircraft.lambda_co4)

    CL0_w = aircraft.AE_Cl0 * (np.cos(aircraft.lambda_co4)) ** 2 # CL0 of wing, ADSEE-II L1 slide 61
    CL0_tot = CL0_w + 1.067171 # TODO: Change flap constant - CL0 of wing with full flaps, constant is contribution of flaps

    dCm_ac_fus = (-1.8) * (1 - 2.5 * aircraft.w_out / aircraft.l_f) * ((np.pi * aircraft.w_out * aircraft.h_out * aircraft.l_f)/(4 * aircraft.Sw * aircraft.MAC_length)) * (CL0_tot / aircraft.CLa_Ah_approach)

    Cm_ac = Cm_ac_w + dCm_ac_f + dCm_ac_fus

    ControlFrac = 1 / ((CL_h / CL_Ah) * (aircraft.AE_l_h / aircraft.MAC_length) * aircraft.AE_Vh_V ** 2)
    ControlSh_S = ControlFrac * (xcgRange + Cm_ac / CL_Ah - aircraft.x_ac_approach)

    return ControlSh_S


def stability_curve(aircraft, xcgRange):
    
    # Making stability line

    StabilityFrac               = 1 / ((aircraft.CLa_h_cruise / aircraft.CLa_Ah_cruise) * (1 - aircraft.AE_dEpsilondA) * (aircraft.AE_l_h/aircraft.MAC_length) * aircraft.AE_Vh_V ** 2)
    StabilityMargin             = 0.05
    StabilitySh_S_margin        = StabilityFrac * xcgRange - StabilityFrac * (aircraft.x_ac_cruise - StabilityMargin)
    StabilitySh_S               = StabilityFrac * xcgRange - StabilityFrac * aircraft.x_ac_cruise

    return StabilitySh_S, StabilitySh_S_margin


##################################################################################################################
"Plotting and running"
##################################################################################################################

def plot_scissor_plot(aircraft):
    xcgRange                    = np.arange(-0.1005, 1.005, 0.005)
    ControlSh_s = controlability_curve(aircraft, xcgRange)
    StabilitySh_S, StabilitySh_S_margin, = stability_curve(aircraft, xcgRange)

    fig1, ax1 = plt.subplots(figsize=(15, 8))
    ax1.plot(xcgRange, StabilitySh_S_margin, color = 'black',
         path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])
    ax1.plot(xcgRange, StabilitySh_S, color = 'red',
         path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])

    ax1.plot(xcgRange, ControlSh_s, color = 'blue',
         path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])

    Sh_S_cont = controlability_curve(aircraft, aircraft.X_cg_fwd)
    Sh_S_stab = stability_curve(aircraft, aircraft.X_cg_aft)[1]

    if Sh_S_cont > Sh_S_stab:
        aircraft.AE_Sh_S = Sh_S_cont
    else:
        aircraft.AE_Sh_S = Sh_S_stab

    x_cg_limit = [aircraft.X_cg_fwd, aircraft.X_cg_aft]
    S_h_S_array = [aircraft.AE_Sh_S, aircraft.Sh_S]

    ax1.plot(x_cg_limit, S_h_S_array, color = 'orange')

    plt.xlim([0, 1])
    plt.ylim([0, 0.6])

    plt.xlabel("x/c [-]")
    plt.ylabel("S_h/S [-]")

    ax1.grid()

    # fig1.savefig("scissorplot")

    plt.show()
    

def run(aircraft):
    Aerodynamic_centre_determination(aircraft)
    plot_scissor_plot(aircraft)

run(aircraft)