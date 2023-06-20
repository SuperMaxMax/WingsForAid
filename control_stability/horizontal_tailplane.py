import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects
from scipy.integrate import quad

sys.path.append('..')

from parameters import atmosphere#, UAV
# aircraft = UAV('aircraft')
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

def LiftRateCoefficient(aircraft, Mach, A, sweep_co2):  # lift rate coefficient tail/wing
    """Lift rate coefficient of a wing(/tail)
    Equation from SEAD Lecture 7, slide 41"""
    aircraft.CS_beta = np.sqrt(1 - Mach ** 2)
    CLa = 2 * np.pi * A / (2 + np.sqrt(4 + ((A * aircraft.CS_beta / aircraft.CS_eta)** 2) * (1 + np.tan(sweep_co2) ** 2  / aircraft.CS_beta ** 2)))
    return CLa

def TaillessLiftRateCoefficient(aircraft, CLa): 
    """Lift rate coefficient of aircraft without tail
    Equation from SEAD Lecture 7, slide 42"""
    CLa_Ah = CLa * (1 + 2.15 * aircraft.w_out / aircraft.b) * aircraft.Sw / aircraft.Sw + (np.pi * aircraft.w_out ** 2) / (2 * aircraft.Sw)

    return CLa_Ah


def nacelle_influence(aircraft, CLa_Ah):
    """a.c. shift due to nacelles
    Equation from SEAD Lecture 7, slide 38"""
    #dx_ac_n = (-4) * (aircraft.w_out ** 2 * (aircraft.X_LEMAC + 0.25 * aircraft.MAC_length)) / (aircraft.Sw * aircraft.MAC_length * CLa_Ah)   
    l_p = aircraft.X_LEMAC + 0.25 * aircraft.MAC_length
    dx_ac_n = -0.05 * aircraft.CS_n_blades * (aircraft.prop_radius * 2)**2 * l_p / (aircraft.Sw * aircraft.MAC_length * CLa_Ah)
    return dx_ac_n


def WingFuselageAC(aircraft, CLa_Ah):
    """Aerodynamic center of wingfuselage combination, x_ac_wf [ ] (already divided by MAC)
    Equation from SEAD Lecture 7, slide 37"""
    x_ac_wf = aircraft.CS_x_ac_w - (1.8 / CLa_Ah) * (aircraft.w_out *aircraft.h_out * (aircraft.X_LEMAC - aircraft.x_lemac)) / (aircraft.Sw * aircraft.MAC_length) + \
    (0.273 / (1 + aircraft.taper)) * ((aircraft.w_out * (aircraft.Sw / aircraft.b) * (aircraft.b - aircraft.w_out)) / ((aircraft.MAC_length ** 2) * (aircraft.b + 2.15 * aircraft.w_out))) * np.tan(aircraft.sweep_co4)
    return x_ac_wf

def Aerodynamic_centre_determination(aircraft):
    aircraft.Mcruise            = Mach_calculation(atm, aircraft.V_cruise, aircraft.h_cruise)
    aircraft.Mmin               = Mach_calculation(atm, aircraft.V_s_min, aircraft.h_TO)
    
    aircraft.CLa_h_cruise       = LiftRateCoefficient(aircraft, aircraft.Mcruise, aircraft.AE_A_h, aircraft.sweep_co2_h)

    aircraft.CLa_w_cruise       = LiftRateCoefficient(aircraft, aircraft.Mcruise, aircraft.A, aircraft.sweep_co2)
    aircraft.CLa_w_approach     = LiftRateCoefficient(aircraft, aircraft.Mmin, aircraft.A, aircraft.sweep_co2)

    aircraft.CLa_Ah_cruise      = TaillessLiftRateCoefficient(aircraft, aircraft.CLa_w_cruise)
    aircraft.CLa_Ah_approach    = TaillessLiftRateCoefficient(aircraft, aircraft.CLa_w_approach)

    aircraft.dx_ac_n_cruise     = nacelle_influence(aircraft, aircraft.CLa_Ah_cruise)
    aircraft.dx_ac_n_approach   = nacelle_influence(aircraft, aircraft.CLa_Ah_approach)

    aircraft.x_ac_wf_cruise     = WingFuselageAC(aircraft, aircraft.CLa_Ah_cruise)
    aircraft.x_ac_wf_approach   = WingFuselageAC(aircraft, aircraft.CLa_Ah_approach)

    aircraft.x_ac_cruise        = aircraft.x_ac_wf_cruise + aircraft.dx_ac_n_cruise
    aircraft.x_ac_approach      = aircraft.x_ac_wf_approach + aircraft.dx_ac_n_approach

def C_m_alpha_calculation(aircraft, x_cg):
    aircraft.C_N_h_alpha = 2 * np.pi * aircraft.AE_A_h/(aircraft.AE_A_h + 2)
    aircraft.C_m_alpha = aircraft.CLa_w_cruise * (x_cg - aircraft.CS_x_ac_w) - aircraft.C_N_h_alpha * (1 - aircraft.AE_dEpsilondA) * (aircraft.AE_Vh_V)**2 * aircraft.Sh_S * aircraft.l_h/aircraft.MAC_length

#################################################################################################################
"FIXME: Determine Control Surface Coefficients"
#################################################################################################################


def Flaplength(aircraft, taperratio, rootchord, span, HLDroot, flappedsurface):
    """Calculate the end location of the flaps based on a starting point measured from the root chord

    taperratio: taper ratio of the wing [m]
    rootchord: chord length of the wing root [m]
    span: total span of the wing [m]
    HLDroot: starting point in span direction of the HLD measured from rootchord[m]
    flappedsurface: total area of wing that is flapped by HLDs [m^2]"""

    def Chordlength(y):
        return aircraft.rootchord * (1 - ((1-aircraft.taper) / (aircraft.b / 2)) * y)
    
    AvailableArea = aircraft.Sw - 2 * quad(Chordlength, 0, HLDroot)[0]  # accounting for overlap of fuselage
    
    if flappedsurface > (AvailableArea):
        print(f"\nUnfeasible flap size, flapped area: {flappedsurface} > wing area not occupied by fuselage: {round(AvailableArea,5)}")
    else:
        A = (-1) * (1 - taperratio) / span
        B = 1
        C = (-1) * ((flappedsurface/2) / rootchord + HLDroot - ((1-taperratio)/span) * (HLDroot ** 2))

        y1 = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        # y2 = (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        print(f"\nFlap ends at {round(y1, 4)} [m] and starts at {round(HLDroot, 4)} [m] --> flap length: {round(y1 - HLDroot, 4)} [m]")
        print(f"Total flapped area: {round(flappedsurface, 4)} [m^2]")

        return y1

def flaps(aircraft):
    flaptype = 'singleslotted'              # Can be 'plain', 'singleslotted' or 'fowler'
    CS_Swf = 0.7 * aircraft.Sw              # spanwise portion of wing influenced by flaps (ADSEE-II, L3 S31) NOTE: Used to calculate resulting dCLmax
    CS_lambda_hinge = 0.05468               # hinge line sweep angle, likely parallel to aft spar [rad] NOTE: For aft spar @ 0.8 chord

    # Using data from Torenbeek aroung page 533
    if flaptype == 'plain':  # FIXME values independent of flap chord, check if logical
        CS_deltaf_TO                   = 20  # flap deflection angle at take-off [deg] (ADSEE-II, L3 S14)
        CS_deltaf_LD                   = 60  # flap deflection angle at landing [deg]
        CS_fc_c                        = aircraft.xc_aft_spar  # flap chord length / wing chord length [-] (ADSEE-II, L3 S32) NOTE: This influences aft spar position
        CS_cprime_c_TO = 1  # wing chord length with take-off extended flaps / chord [-]
        CS_cprime_c_LD = 1  # wing chord length with landing extended flaps / chord [-]
        CS_dClmax_TO = 0.9 * 0.7  # Additional airfoil lift at take-off due to plain flap (ADSEE-II, L3 S36)
        CS_dClmax_LD = 0.9  # Additional airfoil lift at landing due to plain flap
        print("\nPlain flap")
        print(f"dClmax_TO: {CS_dClmax_TO}, dClmax_LD: {CS_dClmax_LD}")

    if flaptype == 'singleslotted':
        CS_deltaf_TO                   = 20  # flap deflection angle at take-off [deg] (ADSEE-II, L3 S14)
        CS_deltaf_LD                   = 40  # flap deflection angle at landing [deg]
        CS_fc_c                        = aircraft.xc_aft_spar  # flap chord length / wing chord length [-] (ADSEE-II, L3 S32) NOTE: This influences aft spar position
        CS_dc_cf_TO                    = 0.2  # increase in chord length / flap chord length [-] (ADSEE-II, L3 S37)
        CS_dc_cf_LD                    = 0.3  # increase in chord length / flap chord length [-]
        CS_cprime_c_TO = 1  # wing chord length with take-off extended flaps / chord [-]
        CS_cprime_c_LD = 1  # wing chord length with landing extended flaps / chord [-]
        CS_dClmax_TO = 1.3 * 0.7  # Additional airfoil lift at take-off due to single slotted flap (ADSEE-II, L3 S36)
        CS_dClmax_LD = 1.3  # Additional airfoil lift at landing due to single slotted flap
        print("\nSingle slotted flap")
        print(f"dClmax_TO: {CS_dClmax_TO}, dClmax_LD: {CS_dClmax_LD}")

    if flaptype == 'fowler':
        CS_deltaf_TO                   = 15  # flap deflection angle at take-off [deg] (ADSEE-II, L3 S15)
        CS_deltaf_LD                   = 40  # flap deflection angle at landing [deg]
        CS_fc_c                        = 0.35  # flap chord length / wing chord length [-] (ADSEE-II, L3 S32) NOTE: This influences aft spar position
        CS_dc_cf_TO                    = 0.48  # increase in chord length / flap chord length [-] (ADSEE-II, L3 S37)
        CS_dc_cf_LD                    = 0.625  # increase in chord length / flap chord length [-]
        CS_cprime_c_TO = 1 + CS_dc_cf_TO * CS_fc_c  # wing chord length with take-off extended flaps / chord [-] (Using notes of ADSEE-II, L3 S37)
        CS_cprime_c_LD = 1 + CS_dc_cf_LD * CS_fc_c  # wing chord length with landing extended flaps / chord [-]
        CS_dClmax_TO = 1.3 * 0.7 * CS_cprime_c_TO   # Additional airfoil lift at take-off due to fowler flap (ADSEE-II, L3 S36)
        CS_dClmax_LD = 1.3 * CS_cprime_c_LD  # Additional airfoil lift at landing due to fowler flap
        print('\nFowler flap')
        print(f"dClmax_TO: {CS_dClmax_TO}, dClmax_LD: {CS_dClmax_LD}")

    print(f"\nCLmax current wing in clean configuration: {aircraft.CL_max_clean}")

    calcCL = False # calculate dCLmax for a given Swf, or vice versa

    if calcCL:
        CS_dCLmax_TO = 0.9 * CS_dClmax_TO * (CS_Swf/aircraft.Sw) * np.cos(CS_lambda_hinge)  # Increase in CL due to flap extension at take-off [-] (ADSEE-II, L3 S35)
        CS_dCLmax_LD = 0.9 * CS_dClmax_LD * (CS_Swf/aircraft.Sw) * np.cos(CS_lambda_hinge)  # Increase in CL due to flap extension at landing [-]
        print(f"CS_dCLmax_TO: {CS_dCLmax_TO}, CS_dCLmax_LD: {CS_dCLmax_LD}")

    else:
        CS_dCLmax_TO = 0
        CS_dCLmax_LD = 0

        if (aircraft.FP_CL_max_TO - aircraft.CL_max_clean) > 0: CS_dCLmax_TO = (aircraft.FP_CL_max_TO - aircraft.CL_max_clean)
        if (aircraft.FP_CL_max_land - aircraft.CL_max_clean) > 0: CS_dCLmax_LD = (aircraft.FP_CL_max_land - aircraft.CL_max_clean)

        CS_Swf_TO = CS_dCLmax_TO / (0.9 * CS_dClmax_TO * np.cos(CS_lambda_hinge)) * aircraft.Sw
        CS_Swf_LD = CS_dCLmax_LD / (0.9 * CS_dClmax_LD * np.cos(CS_lambda_hinge)) * aircraft.Sw
        
        CS_Swf = max(CS_Swf_TO, CS_Swf_LD) # most constraining case becomes required flapped area
        
        print(f"\nNew CLmax during take-off: {aircraft.CL_max_clean + CS_dCLmax_TO}\nNew CLmax during landing: {aircraft.CL_max_clean + CS_dCLmax_LD}")
        print(f"\nSwf_TO: {round(CS_Swf_TO / aircraft.Sw * 100, 3)}% Sw") 
        print(f"Swf_LD: {round(CS_Swf_LD / aircraft.Sw * 100, 3)}% Sw")

    CS_yend_f = Flaplength(aircraft, aircraft.taper, aircraft.rootchord, aircraft.b, 0, CS_Swf)
    return CS_dClmax_LD, CS_cprime_c_LD, CS_Swf, CS_dCLmax_LD, CS_yend_f

       


#################################################################################################################
"Controlability and Stability Curves"
################################################################################################################

def controlability_curve(aircraft, xcgRange):
    """For the controllability curve, the CL_h, CL_Ah and Cm_ac of the aircraft needs to be found. 
    For CL_h: SEAD L8 S17
    For CL_Ah: L = W @ approach
    For Cm_ac: SEAD L8 S19"""
    
    CS_dClmax_LD, CS_cprime_c_LD, CS_Swf, CS_dCLmax_LD, CS_yend_f = flaps(aircraft)
    aircraft.yend_flap = CS_yend_f
    
    CL_Ah = aircraft.W_TO * aircraft.g0 / (0.5 * atm.rho0 * (1.3 * aircraft.V_s_min)**2 * aircraft.Sw)

    #Calculation of Cm_ac starting here
    Cm_ac_w = aircraft.af_cm0 * (aircraft.A * (np.cos(aircraft.sweep_co4))**2 / (aircraft.A + 2 * np.cos(aircraft.sweep_co4)))

    dCm_ac_f = aircraft.CS_mu2 * ((-aircraft.CS_mu1) * CS_dClmax_LD * CS_cprime_c_LD - (CL_Ah + CS_dClmax_LD * (1 - CS_Swf / aircraft.Sw)) * (1/8) * CS_cprime_c_LD * (CS_cprime_c_LD - 1)) \
        + 0.7 * (aircraft.A / (1 + 2 / aircraft.A)) * aircraft.CS_mu3 * CS_dClmax_LD * np.tan(aircraft.sweep_co4) - CL_Ah * (0.25 - aircraft.x_ac_approach)

    CL0_w = aircraft.af_Cl0 * (np.cos(aircraft.sweep_co4)) ** 2 # CL0 of wing, ADSEE-II L1 slide 61
    CL0_tot = CL0_w + CS_dCLmax_LD

    dCm_ac_fus = (-1.8) * (1 - 2.5 * aircraft.w_out / aircraft.l_f) * ((np.pi * aircraft.w_out * aircraft.h_out * aircraft.l_f)/(4 * aircraft.Sw * aircraft.MAC_length)) * (CL0_tot / aircraft.CLa_Ah_approach)

    Cm_ac = Cm_ac_w + dCm_ac_f + dCm_ac_fus

    ControlFrac = 1 / ((aircraft.CL_h / CL_Ah) * (aircraft.l_h / aircraft.MAC_length) * aircraft.AE_Vh_V ** 2)
    # print(f'Cm_ac:{Cm_ac}')
    # print(f"CL_Ah:{CL_Ah}")
    ControlSh_S = ControlFrac * (xcgRange + Cm_ac / CL_Ah - aircraft.x_ac_approach)

    return ControlSh_S


def stability_curve(aircraft, xcgRange):
    
    # Making stability line

    StabilityFrac = 1 / ((aircraft.CLa_h_cruise / aircraft.CLa_Ah_cruise) * (1 - aircraft.AE_dEpsilondA) * (aircraft.l_h/aircraft.MAC_length) * aircraft.AE_Vh_V ** 2)
    StabilityMargin = 0.05
    StabilitySh_S_margin = StabilityFrac * xcgRange - StabilityFrac * (aircraft.x_ac_cruise - StabilityMargin)
    StabilitySh_S = StabilityFrac * xcgRange - StabilityFrac * aircraft.x_ac_cruise

    return StabilitySh_S, StabilitySh_S_margin


##################################################################################################################
"Plotting and running"
##################################################################################################################

def plot_scissor_plot(aircraft, plot):
    xcgRange = np.arange(-0.105, 1.005, 0.0005)
    ControlSh_s = controlability_curve(aircraft, xcgRange)

    StabilitySh_S, StabilitySh_S_margin, = stability_curve(aircraft, xcgRange)

    # Finding smallest Sh/s for control and stability lines
    AbsXcgCont = abs(xcgRange - aircraft.X_cg_fwd)  # Array with distances between index and forward cg
    MinControlSh_s = ControlSh_s[AbsXcgCont.argmin()-1]  # Sh/S for most forward cg of control line, -1 to account for step inaccuracy
    AbsXcgStab = abs(xcgRange - aircraft.X_cg_aft)  # Array with distances between index and aft cg
    MinStabSh_s = StabilitySh_S_margin[AbsXcgStab.argmin()+1]  # Sh/S for most aft cg of control line, +1 to account for step inaccuracy

    if MinControlSh_s > MinStabSh_s:  # TODO FIXME DIT MOET VGM ALTIJD GEBEUREN
        aircraft.Sh_S = MinControlSh_s
    else:
        aircraft.Sh_S = MinStabSh_s

    aircraft.Sh = aircraft.Sh_S * aircraft.Sw
    if plot:
        fig1, ax1 = plt.subplots(figsize=(14, 7))
        # Stability lines
        ax1.plot(xcgRange, StabilitySh_S_margin, label = 'Stability Curve', color = 'black',
            path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])
        ax1.plot(xcgRange, StabilitySh_S, label = 'Neutral Stability Curve', color = 'red',
            path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])
        # Controllability line
        ax1.plot(xcgRange, ControlSh_s, label = 'Controllability Curve', color = 'blue',
            path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])
        # CG range
        ax1.plot([aircraft.X_cg_fwd, aircraft.X_cg_aft], [aircraft.Sh_S, aircraft.Sh_S], '-o', markersize=4, label='CG range loading diagram', color = 'orange')

        ax1.spines['bottom'].set_position('zero')
        ax1.set_ylim(bottom = 0., top=0.4)
        ax1.set_xlim(left = 0., right=1.)
        ax1.set_ylabel(r"$\dfrac{S_h}{S}$",rotation = 0, fontsize = 12)
        ax1.set_xlabel(r"$\dfrac{x_{cg}}{\bar{c}}$", fontsize = 12)
        ax1.yaxis.set_label_coords(-0.05, 0.95)
        ax1.xaxis.set_label_coords(1.03, -0.00)
        ax1.legend(loc='lower right')
        ax1.grid()

        plt.show()
        # fig1.savefig("scissorplot")
    

def hor_run(aircraft, plot):
    Aerodynamic_centre_determination(aircraft)
    C_m_alpha_calculation(aircraft, aircraft.X_cg_full)
    plot_scissor_plot(aircraft, plot)
    
