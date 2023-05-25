import numpy as np
import sys
sys.path.append("..")

from parameters import UAV
aircraft = UAV('aircraft')

#################################################################################################################
"Determine coefficients"
#################################################################################################################
# Mach, CLa, CLa


def LiftRateCoefficient(aircraft, Mach, A, lamda_co5):  # lift rate coefficient tail/wing
    """Lift rate coefficient of a wing(/tail)
    Equation from Lecture 7, slide 41"""
    aircraft.CS_beta = np.sqrt(1 - Mach ** 2)
    CLa = 2 * np.pi * A / (2 + np.sqrt(4 + ((A * aircraft.CS_beta / aircraft.CS_eta)** 2) * (1 + np.tan(lamda_co5 * np.pi / 180) ** 2  / aircraft.CS_beta ** 2)))
    return CLa

def TaillessLiftRateCoefficient(aircraft, CLa): 
    """Lift rate coefficient of aircraft without tail
    Equation from Lecture 7, slide 42"""
    CLa_Ah = CLa * (1 + 2.15 * aircraft.w_out / aircraft.b) * aircraft.Sw / aircraft.Sw + (np.pi * aircraft.w_out ** 2) / (2 * aircraft.Sw)
    return CLa_Ah


def nacelle_influence(aircraft, CLa_Ah): # TODO: BRAM CONTINUE HERE
    """a.c. shift due to nacelles
    Equation from Lecture 7, slide 38"""
    dx_ac_n = (-4) * (aircraft.w_out ** 2 * (aircraft.X_LEMAC + 0.25 * aircraft.MAC_length)) / (aircraft.S * aircraft.c * CLa_Ah)   
    return dx_ac_n


def WingFuselageAC(aircraft, CLa_Ah, c_g):
    """Aerodynamic center of wingfuselage combination, x_ac_wf [ ] (already divided by MAC)
    Equation from Lecture 7, slide 37"""
    x_ac_wf = aircraft.x_ac_w - (1.8 / CLa_Ah) * (aircraft.b_f *aircraft. h_f * aircraft.l_fn) / (aircraft.Sw * aircraft.c) + \
    (0.273 / (1 + aircraft.taperratio)) * ((aircraft.b_f * c_g * (aircraft.b - aircraft.b_f)) / ((aircraft.c ** 2) * (aircraft.b + 2.15 * aircraft.b_f))) * np.tan(aircraft.Lambda_025c * np.pi / 180)
    return x_ac_wf

def Aerodynamic_centre_determination(aircraft):
    aircraft.CLa_h_cruise       = LiftRateCoefficient(aircraft, aircraft.Mcruise, aircraft.Ah, aircraft.lambda_co2_h)

    aircraft.CLa_w_cruise       = LiftRateCoefficient(aircraft, aircraft.Mcruise, aircraft.A, aircraft.lambda_co2)
    aircraft.CLa_w_approach     = LiftRateCoefficient(aircraft, aircraft.Mmin, aircraft.A, aircraft.lambda_co2)

    aircraft.CLa_Ah_cruise      = TaillessLiftRateCoefficient(aircraft, aircraft.CLa_w_cruise)
    aircraft.CLa_Ah_approach    = TaillessLiftRateCoefficient(aircraft, aircraft.CLa_w_approach)

    aircraft.dx_ac_n_cruise     = nacelle_influence(aircraft, aircraft.CLa_Ah_cruise)
    aircraft.dx_ac_n_approach   = nacelle_influence(aircraft, aircraft.CLa_Ah_approach)

    aircraft.x_ac_wf_cruise     = WingFuselageAC(aircraft, aircraft.CLa_Ah_cruise, aircraft.cg)
    aircraft.x_ac_wf_approach   = WingFuselageAC(aircraft, aircraft.CLa_Ah_approach, aircraft.cg)

    aircraft.x_ac_cruise        = aircraft.x_ac_wf_cruise + aircraft.dx_ac_n_cruise
    aircraft.x_ac_approach      = aircraft.x_ac_wf_approach + aircraft.dx_ac_n_approach


#################################################################################################################
"Controlability and Stability Curves"
################################################################################################################

def controlability_curve(aircraft): #TODO: change constants here 
    CL_h                        = -0.35 * aircraft.A_h ** (1/3)
    CL_Ah                       = 22350 * 9.80665 / (0.5 * 1.225 * (95 * 0.51444)**2 * aircraft.Sw)

    Cm_ac_w                     = -0.073 * (aircraft.A * (np.cos(aircraft.lambda_co4 * np.pi / 180))**2 / (aircraft.A + 2 * np.cos(aircraft.lambda_co4 * np.pi / 180)))


def stability_curve(aircraft):
    xcgRange                    = np.arange(0, 1.005, 0.005)
    # Making stability line
    StabilityFrac               = 1 / ((aircraft.CLa_h_cruise / aircraft.CLa_Ah_cruise) * (1 - aircraft.dEpsilondA) * (aircraft.l_h/aircraft.c) * aircraft.Vh_V ** 2)
    StabilityMargin             = 0.05
    StabilitySh_S_margin        = StabilityFrac * xcgRange - StabilityFrac * (aircraft.x_ac_cruise - StabilityMargin)
    StabilitySh_S               = StabilityFrac * xcgRange - StabilityFrac * aircraft.x_ac_cruise

    return StabilitySh_S, StabilitySh_S_margin


##################################################################################################################
"Plotting and running"
##################################################################################################################

def plot_scissor_plot(aircraft):
    controlability_curve(aircraft)
    StabilitySh_S, StabilitySh_S_margin, = stability_curve(aircraft)


def run(aircraft):
    Aerodynamic_centre_determination(aircraft)

    plot_scissor_plot(aircraft)
