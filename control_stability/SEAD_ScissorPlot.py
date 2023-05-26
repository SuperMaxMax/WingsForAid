import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects

from SEAD_CG_calculation import loading_diagram
"""
Calculate scissor plot:
"""

#Include loading diagram
loading = True
# ~~~Constants~~~
# General
b_f = 2.865  # fuselage width [m]
h_f = 2.57  # fuselage height [m]
l_f = 27.165  #fuselage length [m]
l_fn = 11.382  # length from nose to LE wingroot [m]
dEpsilondA = 0  # Downwash gradient effect [ ]
l_h = 13.66  # Length from qc mac wing to qc mac tail (see Lec7 slide 37) [m]

# Wing properties
c = 2.45  # MAC length [m]
x_ac_w = 0.25  # Location of ac of wing, divided by MAC [ ] NOTE: for cruise velocity, found using L7 slide 34
c_t = 1.59  # tip chord length [m]
c_r = 2.57  # root chord length [m]
Lambda = c_t / c_r  # taper ratio [ ]
S = 61  # Wing surface area [m2]
S_net = S - c_r * b_f  # Wing surface area - fuselage overlap [m2]
b = 27.05  # wing span [m]
A = 12  # Wing Aspect ratio [ ]
c_g = S / b  # mean geometric chord [m]
Lambda_025c = 2.5  # quarter chord sweep angle [deg] 
Lambda_05c = 1  # half chord sweep angle [deg]

# Nacelle properties
b_n = 0.88  # nacelle/engine width [m]
l_n = 2.497  # length from front of nacelle to quarterchord sweep [m]

# Tail properties
Vh_V = 1  # Velocity ration Vh/V, set to 1 as T-tail configuration [ ]
Sh = 11.73  # Horizontal tail surface area [m2]
b_h = 7.34  # Horizontal tail span [m]
A_h = b_h ** 2 / Sh  # Horizontal tail aspect ratio [ ]
Lambda_05ch = 4  # Horizontal tail half chord sweep anlge [deg]
eta = 0.95  # Airfoil efficiency constant [ ] NOTE: Found in slide 41 of Lec7

# Mach calculations @ cruise and minimum velocity
# From Jane's: Vcruise = 509 km/h @ FL170 (97% MTOW)
# --> T = 254.46959999999999, p = 52715.4203698557, rho  = 0.7216720192950535
Mcruise = 509 / 3.6 / np.sqrt(1.4 * 287.053 * 254.4696)
Mmin = 95 * 0.51444 / np.sqrt(1.4 * 287.053 * 288.15)  # NOTE: See Excel for source used

#For loading diagram
reference_aircraft = {
    "Name" : "ATR72-600_reference", #Name of aircraft
    "OEW": 13178.4,                   #kg
    "MTOW": 22800.00,               #kg
    "Fuel_max": 5000,               #kg
    "max_payload": 7500,            #kg
    "pax": 72,                      #Number of passengers
    "pax_weight": 95,               #Weight of 1 passengers [kg]
    "front_volume": 4.6,            #Volume of front cargo bay [m^3]
    "aft_volume": 4.8,              #Volume of aft cargo bay [m^3]
    "MAC": 2.45,                    #Mean Aerodynamic Chord of aircraft
    "xlemac": 11.43,                #x position of the leading edge of the mean aerodynamic chord
    "xcg_OEW": 1.027711073,         #CG location OEW /MAC
    "seat_pitch": 29,               #Seat pitch [inch]
    "xcg_front_seat": 6.25,         #x position first seat from nose [m]
    "xcg_front_CB": 4.358164881,    #xcg front cargo bay [m]
    "xcg_aft_CB": 21.47884442,      #xcg aft cargo bay [m]
    "xcg_fuel": 12.655              #xcg of the fuel tanks [m]
}

modified_aircraft = {
    "Name" : "ATR72-600_modified", #Name of aircraft
    "OEW"  : 13837.13,             #kg
    "MTOW" : 22800.00,             #kg
    "Fuel_max" : 5000,             #kg #This needs to be reconsidered, less fuel
    "max_payload" : 7500 - 760,          #kg
    "pax" : 64,                    #Number of passengers
    "pax_weight" : 95,             #Weight of 1 passengers [kg]
    "front_volume" : 4.6,          #Volume of front cargo bay [m^3]
    "aft_volume" : 4.8,            #Volume of aft cargo bay [m^3]
    "MAC" : 2.45,                  #Mean Aerodynamic Chord of aircraft
    "xlemac" : 11.43,              #x position of the leading edge of the mean aerodynamic chord
    "xcg_OEW" : 1.300000062,       #CG location OEW /MAC
    "seat_pitch" : 29,             #Seat pitch [inch]
    "xcg_front_seat" : 6.25,       #x position first seat from nose [m]
    "xcg_front_CB" : 4.358164881,  #xcg front cargo bay [m]
    "xcg_aft_CB" : 21.47884442,    #xcg aft cargo bay [m]
    "xcg_fuel" : 12.655            #xcg of the fuel tanks [m]
}

def Scissorplot(AR, b_nac, l_nac):
    A = AR
    b_n = b_nac
    l_n = l_nac
    def WingFuselageAC(x_ac_w, CLa_Ah, b_f, h_f, l_fn, S, b, c, taperratio, Lambda_025c, c_g):
        """Aerodynamic center of wingfuselage combination, x_ac_wf [ ] (already divided by MAC)
        Equation from Lecture 7, slide 37"""
        x_ac_wf = x_ac_w - (1.8 / CLa_Ah) * (b_f * h_f * l_fn) / (S * c) + \
        (0.273 / (1 + taperratio)) * ((b_f * c_g * (b - b_f)) / ((c ** 2) * (b + 2.15 * b_f))) * np.tan(Lambda_025c * np.pi / 180)
        return x_ac_wf

    def LiftRateCoefficient(Mach, A, Lambda_05c):  # lift rate coefficient tail/wing
        """Lift rate coefficient of a wing(/tail)
        Equation from Lecture 7, slide 41"""
        beta = np.sqrt(1 - Mach ** 2)
        CLa = 2 * np.pi * A / (2 + np.sqrt(4 + ((A * beta / eta)** 2) * (1 + np.tan(Lambda_05c * np.pi / 180) ** 2  / beta ** 2)))
        return CLa

    def TaillessLiftRateCoefficient(CLa, b_f, b, S_net, S): 
        """Lift rate coefficient of aircraft without tail
        Equation from Lecture 7, slide 42"""
        CLa_Ah = CLa * (1 + 2.15 * b_f / b) * S_net / S + (np.pi * b_f ** 2) / (2 * S)
        return CLa_Ah

    # Lift rate coefficients Cruise
    CruiseCLa_h = LiftRateCoefficient(Mcruise, A_h, Lambda_05ch)  # Lift rate coefficient horizontal tail [1/rad]
    print("CruiseCLa_h", CruiseCLa_h)
    CruiseCLa_w = LiftRateCoefficient(Mcruise, A, Lambda_05c)  # Lift rate coefficient wing [1/rad]
    print("CruiseCLa_w:", CruiseCLa_w)

    CruiseCLa_Ah = TaillessLiftRateCoefficient(CruiseCLa_w, b_f, b, S_net, S)  # lift rate coefficient aircraft without tail [1/rad]
    print("CruiseCLa_Ah:", CruiseCLa_Ah)

    # Lift rate coefficients Approach
    ApproachCLa_w = LiftRateCoefficient(Mmin, A, Lambda_05c)
    print("ApproachCLa_w:", ApproachCLa_w)
    
    ApproachCLa_Ah = TaillessLiftRateCoefficient(ApproachCLa_w, b_f, b, S_net, S)
    print("ApproachCLa_Ah:", ApproachCLa_Ah)

    # ac shift due to nacelles - L7 slide 38
    # Cruise
    Cruisedx_ac_n = 2 * (-4) * (b_n ** 2 * l_n) / (S * c * CruiseCLa_Ah)
    print("Cruisedx_ac_n:", Cruisedx_ac_n)
    
    # Approach
    Approachdx_ac_n = 2 * (-4) * (b_n ** 2 * l_n) / (S * c * ApproachCLa_Ah)
    print("Approachdx_ac_n:", Approachdx_ac_n)
    
    # Total x_ac cruise
    Cruise_x_ac_wf = WingFuselageAC(x_ac_w, CruiseCLa_Ah, b_f, h_f, l_fn, S, b, c, Lambda, Lambda_025c, c_g) # ac of wing fuselage combo
    print("Cruise_x_ac_wf", Cruise_x_ac_wf)
    
    Cruise_x_ac = Cruise_x_ac_wf + Cruisedx_ac_n
    print("Cruise_x_ac:", Cruise_x_ac)
    
    # Total x_ac Approach
    Approach_x_ac_wf = WingFuselageAC(x_ac_w, ApproachCLa_Ah, b_f, h_f, l_fn, S, b, c, Lambda, Lambda_025c, c_g) # ac of wing fuselage combo
    print("Approach_x_ac_wf", Approach_x_ac_wf)
    
    Approach_x_ac = Approach_x_ac_wf + Approachdx_ac_n
    print("Approach_x_ac:", Approach_x_ac)

    xcgRange = np.arange(0, 1.005, 0.005)
    # Making stability line
    StabilityFrac = 1 / ((CruiseCLa_h / CruiseCLa_Ah) * (1 - dEpsilondA) * (l_h/c) * Vh_V ** 2)
    StabilityMargin = 0.05
    StabilitySh_S_margin = StabilityFrac * xcgRange - StabilityFrac * (Cruise_x_ac - StabilityMargin)
    StabilitySh_S = StabilityFrac * xcgRange - StabilityFrac * Cruise_x_ac

    # Making the controlability line
    CL_h = -0.35 * A_h ** (1/3)  # TODO: check if correct and if so move to coefficients (from L8 slide 17)
    CL_Ah = 22350 * 9.80665 / (0.5 * 1.225 * (95 * 0.51444)**2 * S) # NOTE: wing lift coefficient in landing configuration (max flaps), assumed that L = W_max_landing, using V = Vmin
    print("Cl_Ah", CL_Ah)

    # Moment contribution of wing
    Cm_ac_w = -0.073 * (A * (np.cos(Lambda_025c * np.pi / 180))**2 / (A + 2 * np.cos(Lambda_025c * np.pi / 180)))
    print("Cm_ac_w:", Cm_ac_w)

    # moment contribution of flaps
    mu1 = 0.24
    mu2 = 0.78
    mu3 = 0.525
    dClmax = 1.7028171
    dCm_ac_f = mu2 * ((-mu1) * dClmax * 1.06426 - (CL_Ah + dClmax * (1 - 42.47695 / S)) * (1/8) * 1.06426 * (1.06426 - 1)) + 0.7 * (A / (1 + 2 / A)) * mu3 * dClmax * np.tan(Lambda_025c * np.pi /180) 
    print("dCm_ac_f:", dCm_ac_f)

    # Moment contribution of fuselage
    CL0_w = 0.454 * (95 * 0.51444 * np.cos(Lambda_025c * np.pi / 180)) ** 2 / (95 * 0.51444) ** 2 # CL0 of wing, ADSEE L1 slide 61

    print("CL0_w", CL0_w)
    
    CL0_tot = CL0_w + 1.067171  # CL0 of wing and flaps, constant is contribution of flaps (see excel)
    print("CL0_tot (landing}:", CL0_tot)
    
    dCm_ac_fus = (-1.8) * (1 - 2.5 * b_f / l_f) * ((np.pi * b_f * h_f * l_f)/(4 * S * c)) * (CL0_tot / ApproachCLa_Ah)
    print("dCm_ac_fus", dCm_ac_fus)

    Cm_ac = Cm_ac_w + dCm_ac_f + dCm_ac_fus # NOTE, L8 slide 19
    print("Cm_ac:", Cm_ac)

    ControlFrac = 1 / ((CL_h / CL_Ah) * (l_h / c) * Vh_V ** 2)
    ControlSh_S = ControlFrac * (xcgRange + Cm_ac / CL_Ah - Approach_x_ac)
    return xcgRange, StabilitySh_S_margin, StabilitySh_S, ControlSh_S

print("\nReference Aircraft Data")
xcgRange, StabilitySh_S_margin_OG, StabilitySh_S_OG, ControlSh_S_OG = Scissorplot(A, b_n, l_n)
print("\nNew Aircraft Data")
xcgRange, StabilitySh_S_margin_NEW, StabilitySh_S_NEW, ControlSh_S_NEW = Scissorplot(A*1.2, b_n*1.25, l_n*1.15)



# Plotting Original Aircraft
fig1, ax1 = plt.subplots(figsize=(15, 8))
fig2, ax2 = plt.subplots(figsize=(15, 8))

ax1.plot(xcgRange, StabilitySh_S_margin_OG, label = 'Ref Stability Curve', color = 'black',
         path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])
ax1.plot(xcgRange, StabilitySh_S_OG, label = 'Ref Neutral Stability Curve', color = 'red')
ax1.plot(xcgRange, ControlSh_S_OG, label = 'Ref Controllability Curve', color = 'blue',
         path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])

ax2.plot(xcgRange, StabilitySh_S_margin_OG, label = 'Ref Stability Curve', linestyle = 'dashed', color = 'black')
ax2.plot(xcgRange, StabilitySh_S_OG, label = 'Ref Neutral Stability Curve', linestyle = 'dashed', color = 'red')
ax2.plot(xcgRange, ControlSh_S_OG, label = 'Ref Controllability Curve', linestyle = 'dashed', color = 'blue')
ax2.plot(xcgRange, StabilitySh_S_margin_NEW, label = 'New Stability Curve', color = 'black',
         path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])
ax2.plot(xcgRange, StabilitySh_S_NEW, label = 'New Neutral Stability Curve', color = 'red')
ax2.plot(xcgRange, ControlSh_S_NEW, label = 'New Controllability Curve', color = 'blue',
         path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])

if loading:
    xcg_aft, xcg_front = loading_diagram(reference_aircraft)
    xcg_loading = [xcg_front, xcg_aft]
    Shs_real = [Sh / S, Sh / S]
    ax1.plot(xcg_loading, Shs_real,'-o', label='Ref xcg range loading diagram', color='orange')
    xcg_aft_new, xcg_front_new = loading_diagram(modified_aircraft)
    xcg_loading_new = [xcg_front_new, xcg_aft_new]
    Shs_real_new = [Sh / S, Sh / S]
    ax2.plot(xcg_loading_new, Shs_real_new, '-o', label='New xcg range loading diagram', color='green')
    ax2.plot(xcg_loading, Shs_real,'|', label='Ref xcg range loading diagram', color='orange', linestyle = 'dashed', markersize=10)

ax1.legend(loc='lower right')
ax2.legend(loc='lower right')

# Move x-axis to x = 0
# ax1 = plt.gca()
ax1.spines['bottom'].set_position('zero')
ax1.set_ylim(bottom = 0., top=0.3)
ax1.set_xlim(left = 0., right=1.)
ax1.set_ylabel(r"$\dfrac{S_h}{S}$",rotation = 0, fontsize = 12)
ax1.set_xlabel(r"$\dfrac{x_{cg}}{\bar{x}}$", fontsize = 12)
ax1.yaxis.set_label_coords(-0.05, 0.95)
ax1.xaxis.set_label_coords(1.03, -0.00)
ax1.grid()

ax2.spines['bottom'].set_position('zero')
ax2.set_ylim(bottom = 0., top=0.3)
ax2.set_xlim(left = 0., right=1.)
ax2.set_ylabel(r"$\dfrac{S_h}{S}$",rotation = 0, fontsize = 12)
ax2.set_xlabel(r"$\dfrac{x_{cg}}{\bar{x}}$", fontsize = 12)
ax2.yaxis.set_label_coords(-0.05, 0.95)
ax2.xaxis.set_label_coords(1.03, -0.00)
ax2.grid()

# fig1.savefig("ref_scissorplot")
# fig2.savefig("new_scissorplot")

plt.show()