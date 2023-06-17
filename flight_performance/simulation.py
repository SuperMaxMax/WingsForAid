import sys
import os.path
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

# Start your import below this
# import packages
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import time
# Import other files
from parameters import UAV, airport, atmosphere

# aircraft = UAV("aircraft")
# airfield = airport("Sudan")
# atm      = atmosphere()
hp_to_watt = 745.699872

def takeoffweight(obj, W_F):
    ZFW = obj.W_OE + obj.n_boxes*obj.boxweight
    TOW = ZFW + W_F
    return TOW

def atm_parameters(atm_obj, h):
    T    = atm_obj.T0 + atm_obj.lambd * h
    rho  = atm_obj.rho0*np.power((T/atm_obj.T0), (-((atm_obj.g / (atm_obj.lambd * atm_obj.R))+1)))
    p    = atm_obj.p0*np.power((T/atm_obj.T0), (-(atm_obj.g / (atm_obj.lambd * atm_obj.R))))
    a    = np.sqrt(atm_obj.gamma*atm_obj.R*T)
    return p, T, rho, a

def dragpolar(ac_obj, CL):
    CD = ac_obj.CD0 + 1.15*CL**2/(np.pi*ac_obj.A*ac_obj.e)
    return CD

def prop_d(obj, power):
    d = 1.5 * (power ** (1/4)) * 0.3084
    obj.prop_radius = d/2
    return d

def climbrate(ac_obj, atm_obj, W_F, V, P_climb, plot=True): 
    start_time = time.time()                                                                # Counter to show computational effort
    V *= 0.5144                                                                             # Multiply V by 0.5144 to convert to m/s from kts
    CL_opt = np.sqrt(3*ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)                                  # Calculate the optimal climb CL
    print(f"The optimal CL for climb is {CL_opt}")                                          
    atm_parameters_vectorized = np.vectorize(lambda h: atm_parameters(atm_obj, h))          # Make a vectorized atm_parameters function
    alt_range = np.arange(0, ac_obj.ceiling, 0.5)                                           # Define an altitude range from 0 to the ceiling
    atm_obj.p, atm_obj.T, atm_obj.rho, atm_obj.a = atm_parameters_vectorized(alt_range)     # Calculate arrays of atmospheric parameters
    W       = takeoffweight(ac_obj, W_F) * atm_obj.g                                        # Calculate take-off weight
    ROC     = np.empty(0)                                                                   # Define empty arrays (used for plotting)
    Gamma   = np.empty(0)
    ROC_max = np.empty(0)
    Gammax  = np.empty(0)
    for i in range(len(alt_range)):
        p, T, rho, a = atm_obj.p[i], atm_obj.T[i], atm_obj.rho[i], atm_obj.a[i]             # Get the atmospheric parameters for altitude i
        CL = 2*W/(rho*ac_obj.Sw*V**2)                                                       # Calculate CL based on input climb speed
        CD = dragpolar(ac_obj, CL)                                                          # Calculate CD based on CL
        V_opt = np.sqrt(2*W/(rho*ac_obj.Sw*CL_opt))                                         # Calculate the optimum climb speed based on CL_opt
        CD_opt = dragpolar(ac_obj, CL_opt)                                                  # Calculate CD based on CL_opt
        D  = 1/2 * rho * V**2 * ac_obj.Sw * CD                                              # Calculate the drag
        Pr = D*V                                                                            # Calculate the power required based on drag and V
        Pr_min = 1/2 * rho * V_opt**3 * ac_obj.Sw * CD_opt                                  # Minimum power required
        Pa = ac_obj.power * P_climb * ac_obj.prop_eff * 745.699872 * (rho/atm_obj.rho0)**(3/4)  # Calculate power available
        roc= (Pa - Pr)/W                                                                    # Achievable rate of climb for altitude i at input speed
        roc_max = (Pa - Pr_min)/W                                                           # Maximum rate of climb for altitude at V_opt
        gamma = (roc/V)*(180/np.pi)                                                         # Climb angle for altitude i at input speed
        gamma_rocmax = (roc_max/V_opt)*(180/np.pi)                                          # Maximum climb angle for altitude i
        FMF= ac_obj.SFC * ac_obj.power * P_climb * 745.699872                               # Instantaneous fuel mass flow
        if i != 0:
            dt = (alt_range[i]-alt_range[i-1]) / (roc)
        else:
            dt = (alt_range[i+1]-alt_range[i]) / (roc) 
        dWF= FMF * dt
        W -= dWF
        ROC     = np.append(ROC, roc)
        Gamma   = np.append(Gamma, gamma)
        ROC_max = np.append(ROC_max, roc_max)
        Gammax  = np.append(Gammax, gamma_rocmax)

    end_time = time.time()
    if plot:
        plt.plot(alt_range, ROC, label=f"Climb power: {P_climb*100}%", color = 'green')
        plt.plot(alt_range, ROC_max, label=f"Max climb rate at power setting: {P_climb*100}%", color = 'red')
        plt.xlabel("Altitude [m]")
        plt.ylabel("Rate of Climb [m/s]")
        plt.legend()
        plt.show()
        plt.plot(alt_range, Gamma, label=f"Climb power: {P_climb*100}%", color= 'green')
        plt.plot(alt_range, Gammax, label=f"Climb angle at roc_max at power setting: {P_climb*100}%", color = 'red')
        plt.xlabel("Altitude [m]")
        plt.ylabel("Climb angle")
        plt.legend()
        plt.show()
    print("----------------------------------------------------------------------------")
    print(f"The max Rate of Climb is {ROC[ROC==np.max(ROC)]} [m/s] at altitude {alt_range[ROC==np.max(ROC)][0]} [m]")
    print(f"The max climb angle is {Gamma[Gamma==np.max(Gamma)]} [deg] at altitude {alt_range[Gamma==np.max(Gamma)][0]} [m]")
    print("----------------------------------------------------------------------------")
    print(f"This calculation took {end_time-start_time} seconds")
    return

def flightceiling(ac_obj, atm_obj, W_F, plot=True, result = False):
    W   = ac_obj.W_TO * atm_obj.g
    h   = 0.0
    Pa  = ac_obj.power * ac_obj.prop_eff * 735.49875
    CL_opt  = np.sqrt(3*ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    CD_opt  = dragpolar(ac_obj, CL_opt)
    V   = np.sqrt(2*W/(atm_obj.rho0*ac_obj.Sw*CL_opt))
    Pr  = 1/2 * atm_obj.rho0 * V**3 * ac_obj.Sw * CD_opt
    roc = (Pa-Pr)/W
    dt  = 1.0
    t   = 0.0
    Time    = np.empty(0)
    ROC     = np.empty(0)
    Weight  = np.empty(0)
    Height  = np.empty(0)
    while roc > 0.508:
        p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
        V   = np.sqrt(2*W/(rho*ac_obj.Sw*CL_opt))
        Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD_opt
        Pa  = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt
        roc = (Pa - Pr)/W
        if t%10 == 0:
            print(f"Altitude: {np.round(h, 2)} [m] | Power required: {np.round(Pr, 2)} [W] | Power available: {np.round(Pa, 2)} [W] | Rate of Climb: {np.round(roc, 2)} [m/s] | Velocity: {np.round(V, 2)} [m/s]")
        h   += roc * dt
        W   -= Pa * ac_obj.SFC * dt * atm_obj.g
        t   += dt
        Time= np.append(Time, t)
        ROC = np.append(ROC, roc)
        Height = np.append(Height, h)
        Weight = np.append(Weight, W)
    if result:
        print("----------------------------------------------------------------------------")
        print(f"The time to get to a cruising altitude of {ac_obj.h_cruise} [m] is {Time[np.abs(Height - ac_obj.h_cruise) <= 20.0][0]} [seconds]")
        print(f"The fuel used to get to a cruising altitude of {ac_obj.h_cruise} [m] is {np.round((Weight[0]-Weight[np.abs(Height - ac_obj.h_cruise) <= 20.0][0])/atm_obj.g)} [kg]")
        print(f"The maximum altitude is equal to {np.round(Height[-1])} [m]")
        print("----------------------------------------------------------------------------")
    if plot:
        # Create a figure and three subplots
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

        ax1.plot(Time, Height, color="red")
        ax1.set_ylabel('Altitude')
        ax1.set_title('Altitude vs Time')

        ax2.plot(Time, ROC, color='green')
        ax2.set_ylabel('Rate of Climb')
        ax2.set_title('Rate of Climb vs Time')

        ax3.plot(Time, Weight/atm_obj.g, color='red')
        ax3.set_xlabel('Time')
        ax3.set_ylabel('Weight')
        ax3.set_title('Weight vs Time')

        plt.tight_layout()

        plt.show()
    return Time, Height

# # ---------------- Assumptions for take-off equations of motion -----------------
# # Wind is included by take it into account in the speed: V_eff = V - V_wind
# # Runway slope is not zero
# # The present of rain is taken into account in the friction coefficient with the ground, mu
# # delta_rw = runway slope
# # D_g = force due to the ground friction, with
# # Thrust and lift are taken as average values


def TO_eom(obj, ap, atmos, max_runwayslope, max_hairport, max_headwind, max_tailwind, W_f, Plot=True):

    CL_all = []
    for i in range(0, 4):
        if i == 0:
            dic_constants = {'runway slope': np.arange(0, max_runwayslope),
                             'airport altitude': 0, 'wing surface area': obj.Sw,
                             'weight': takeoffweight(obj, W_f) * atmos.g,
                             'wind speed': 0, 'propeller power': obj.power * hp_to_watt,
                             'propeller efficiency': obj.prop_eff}
            p, T, rho, a = atm_parameters(atmos, dic_constants['airport altitude'])

        if i == 1:
            dic_constants['runway slope'] = 0
            dic_constants['wind speed'] = np.arange(0, max_headwind)

        if i == 2:
            dic_constants['wind speed'] = np.arange(0, max_tailwind, -1)

        if i == 3:
            dic_constants['wind speed'] = 0
            dic_constants['airport altitude'] = np.arange(0, max_hairport)
            p, T, rho, a = atm_parameters(atmos, dic_constants['airport altitude'])

        CL_max = obj.CL_max_clean
        CL = CL_max
        S = [600]
        while max(S) <= 750:
            S_old = S.copy()
            CL -= 0.01
            V_LOF = 1.05 * (np.sqrt(
                dic_constants['weight'] / dic_constants['wing surface area'] * 2 / rho * 1 / CL_max) - dic_constants[
                        'wind speed'])
            V_avg_sq = V_LOF ** 2 / 2
            CD = obj.CD0 + CL ** 2 / (np.pi * obj.A * obj.e)
            D = CD * V_avg_sq * rho / 2 * dic_constants['wing surface area']
            L = CL * V_avg_sq * rho / 2 * dic_constants['wing surface area']
            T = dic_constants['propeller power'] * dic_constants['propeller efficiency'] / np.sqrt(V_avg_sq)
            D_g = ap.mu_ground * (dic_constants['weight'] * np.cos(np.radians(dic_constants['runway slope'])) - L)
            a = atmos.g / dic_constants['weight'] * (
                        T - D - D_g - dic_constants['weight'] * np.sin(np.radians(dic_constants['runway slope'])))
            S = V_LOF ** 2 / (2 * a)
            if CL <= 0.75:
                break

        S = S_old
        CL += 0.01
        CL_all.append(CL)
        if len(S) == 1:
            print('One or more of the flight conditions; runway slope, headwind, tailwind, airport altitude, '
                  'are not achievable')
            Plot = False
            break
        else:
            CL_to = max(CL_all)
            print(CL_to, CL_max)
            obj.FP_CL_max_TO = CL_max
            obj.FP_CL_TO = CL_to

        if Plot:
            figure, axis = plt.subplots(2, 2)
            if i == 0:
                axis[0, 0].plot(dic_constants['runway slope'], S)
                axis[0, 0].set_title('runway slope vs C_L take-off')
                axis[0, 0].set_xlabel('runway slope[deg]')
                axis[0, 0].set_ylabel('C_L take-off [-]')

            elif i == 1:
                axis[1, 0].plot(dic_constants['wind speed'], S, color='red')
                axis[1, 0].set_title('headwind vs C_L take-off')
                axis[1, 0].set_xlabel('headwind speed [m/sec]')
                axis[1, 0].set_ylabel('C_L take-off [-]')

            elif i == 2:
                axis[1, 1].plot(dic_constants['wind speed'], S, color='green')
                axis[1, 1].set_title('tailwind vs C_L take-off')
                axis[1, 1].set_xlabel('tailwind speed [m/sec]')
                axis[1, 1].set_ylabel('C_L take-off [-]')

            else:
                axis[0, 1].plot(dic_constants['airport altitude'], S, color='black')
                axis[0, 1].set_title('airport altitude vs C_L take-off')
                axis[0, 1].set_xlabel('airport altitude [m]')
                axis[0, 1].set_ylabel('C_L take-off [-]')

            plt.subplots_adjust(hspace=0.6)
            plt.subplots_adjust(wspace=0.5)
            plt.suptitle('Take-Off')
            plt.show()

    return S


# TO_eom(aircraft, airfield, atm, 12, 4000, 12, -13.4, 65)

# Result:
# - If 12 boxes, then the slope limit is 11 degrees and max tailwind of 13 m/s = 25.3 kts
# - If we have a tailwind of 30knt we can't take more than 10 boxwes

# ------------------------------------------------------------------------

def turnperformance(ac_obj, atm_obj, W = None, h = None, V = None):
    if V == None:
        V = ac_obj.V_cruise
    else:
        V = V
    if W == None:
        W = ac_obj.W_TO * ac_obj.W1W_TO * ac_obj.W2W1 * ac_obj.W3W2 * ac_obj.W4W3
    else:
        W = W
    if h == None:
        h = ac_obj.h_cruise
    else:
        h = h
    rho = atm_parameters(atm_obj, h)[2]
    p   = atm_parameters(atm_obj, h)[0]
    # Standard turns
    standard_rates  = np.array([ac_obj.turnrate_half, ac_obj.turnrate_1, ac_obj.turnrate_2]) * (np.pi/180)  # Convert to rad
    turnradius_std  = V/standard_rates                                                                      # Radius of standard turns
    bankangle_std   = np.arctan(V**2/(atm_obj.g * turnradius_std))                                          # Required bank angles
    n_stdrates      = 1/np.cos(bankangle_std)                                                               # Associated load factors
    T_turn_std      = 360/(standard_rates * 180/np.pi)                                                      # Turning time
    print("-------------------------------------------------------------------------")
    print(f"Standard rate 1/2: V = {np.round(V, 2)} [m/s] | Turnradius = {np.round(turnradius_std[0], 2)} [m] | Bank angle = {np.round(bankangle_std[0]*180/np.pi, 2)} [deg] | Load factor = {np.round(n_stdrates[0], 2)} | Turning time = {np.round(T_turn_std[0], 2)} [s]")
    print(f"Standard rate 1  : V = {np.round(V, 2)} [m/s] | Turnradius = {np.round(turnradius_std[1], 2)} [m] | Bank angle = {np.round(bankangle_std[1]*180/np.pi, 2)} [deg] | Load factor = {np.round(n_stdrates[1], 2)} | Turning time = {np.round(T_turn_std[1], 2)} [s]")
    print(f"Standard rate 2  : V = {np.round(V, 2)} [m/s] | Turnradius = {np.round(turnradius_std[2], 2)} [m] | Bank angle = {np.round(bankangle_std[2]*180/np.pi, 2)} [deg] | Load factor = {np.round(n_stdrates[2], 2)} | Turning time = {np.round(T_turn_std[2], 2)} [s]")
    print("-------------------------------------------------------------------------")
    # Steepest turn (drag limited)
    Pa = ac_obj.power * ac_obj.prop_eff * p/atm_obj.p0 * 735.49875                             # Power available in Watts at altitude
    print(f"Power available: {Pa} [W]")
    bankangle_steep = 0
    n_steep  = 1/np.cos(bankangle_steep)
    CL_steep = 2*n_steep*W*atm_obj.g/(rho*V**2*ac_obj.Sw)
    CD_steep = dragpolar(ac_obj, CL_steep)
    Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD_steep
    print(f"Power required {Pr}")
    while Pa > Pr:
        bankangle_steep += np.pi/180                                  # Increase bankangle with 1 deg (expressed in rad) every calculation
        n_steep = 1/(np.cos(bankangle_steep))
        CL_steep = 2*n_steep*W*atm_obj.g/(rho*V**2*ac_obj.Sw)
        if CL_steep > ac_obj.CL_max_clean:
            bankangle_steep -= np.pi/180
            n_steep = 1/(np.cos(bankangle_steep))
            print(f"The maximum bank angle is stall limited")
            break
        CD_steep = dragpolar(ac_obj, CL_steep)
        Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD_steep
    print(f"The maximum bank angle at TAS = {np.round(V, 2)} [m/s], ALT = {np.round(h, 2)} [m] and Weight = {np.round(W, 2)} [kg] is {np.round(bankangle_steep*180/np.pi, 2)} [deg]")
    print("-------------------------------------------------------------------------")
    # Minimum turn radius
    R_min = V**2/(atm_obj.g*np.sqrt(n_steep**2 -1))
    print(f"The minimum turn radius at TAS = {np.round(V, 2)} [m/s], ALT = {np.round(h, 2)} [m] and Weight = {np.round(W, 2)} [kg] is {np.round(R_min, 2)} [m]")
    print("-------------------------------------------------------------------------")
    return
# turnperformance(aircraft, atm)

# -----------------------------------------------------------------

def cruiseperformance(ac_obj, atm_obj, n_boxes, W_F_TO, Range=None, V_cruise=None, h_cruise=None):
    if Range == None:
        R = ac_obj.R
    else:
        R = Range
    if V_cruise == None:
        V_cruise = ac_obj.V_cruise
    else:
        V_cruise = V_cruise
    if h_cruise == None:
        h_cruise = ac_obj.h_cruise
    else:
        h_cruise = h_cruise
    W_cr    = (ac_obj.W_OE + W_F_TO + n_boxes * ac_obj.boxweight) * ac_obj.W1W_TO * ac_obj.W2W1 * ac_obj.W3W2 * ac_obj.W4W3
    rho_cr  = atm_parameters(atm_obj, h_cruise)[2]
    p_cr    = atm_parameters(atm_obj, h_cruise)[0]
    r_it = 0.0
    t    = 0.0
    dt   = 0.1
    W    = W_cr * atm_obj.g
    D_max = 0
    while r_it < R:
        CL_cr   = 2*W/(rho_cr*V_cruise**2*ac_obj.Sw)
        CD_cr   = dragpolar(ac_obj, CL_cr)
        r_it    += (V_cruise * dt)
        t       += dt
        D       = 1/2 * rho_cr * V_cruise**2 * ac_obj.Sw * CD_cr
        if np.max(D) > D_max:
            D_max = np.max(D)
        Pa      = ac_obj.power * ac_obj.prop_eff * p_cr/(atm_obj.p0) * hp_to_watt
        P_req   = D*V_cruise/ac_obj.prop_eff
        if t%10 <= 0.01:
            print(f"Power available: {Pa} [W] | Power required: {P_req} [W]")
        F       = ac_obj.SFC * P_req
        W       -= (F*dt)
    W_F_used = W_cr - W
    kg_kgkm = W_F_used / (n_boxes * ac_obj.boxweight)
    print("---------------------------------------------------")
    print(f"Cruise performance - Range {R/1000} [km] - Cruise speed {V_cruise} [m/s] - Cruise height {h_cruise} [m]")
    print("---------------------------------------------------")
    print(f"The cruise time is {np.round(t, 2)} seconds ({np.round(t/3600, 2)} hours)")
    print(f"The fuel used during the cruise is {np.round(W_F_used)} [kilograms] ({np.round(((W_cr-W)/0.7429), 2)} [L] @ {ac_obj.fueldensity} [kg/m^3])")
    print("---------------------------------------------------")
    print("D", D_max)
    print('V', (ac_obj.power * hp_to_watt * ac_obj.eta_p)/D_max)
    return None
# cruiseperformance(aircraft, atm, 12, 65)

def payloadrange(ac_obj, atm_obj, V_cruise=None, h_cruise=None, plot=True):
    if V_cruise == None:
        V_cruise = ac_obj.V_cruise
    else:
        V_cruise = V_cruise
    if h_cruise == None:
        h_cruise = ac_obj.h_cruise
    else:
        h_cruise = h_cruise
    Fuel_loads  = np.arange(0, ac_obj.fuelcapacity+1.0, 1.0)
    Reserve     = ac_obj.M_res * ac_obj.fuelcapacity * ac_obj.fueldensity
    ZFW         = ac_obj.W_OE + Reserve + ac_obj.n_boxes * ac_obj.boxweight             # ac_ojb.n_boxes must be 12 here
    maxZFW_fuel = ac_obj.W_TO - ZFW
    ZFW_maxfuel = ac_obj.W_TO - ac_obj.fuelcapacity * ac_obj.fueldensity
    PL_maxfuel  = ZFW_maxfuel - Reserve - ac_obj.W_OE
    Ferryweight = ac_obj.W_OE + Reserve + ac_obj.fuelcapacity * ac_obj.fueldensity
    print("-------------------------------------------------------------------")
    print(f"Max ZFW: {np.round(ZFW, 2)} [kg]")
    print(f"Fuel @ max ZFW: {np.round(maxZFW_fuel, 2)} [kg] or {np.round(maxZFW_fuel/ac_obj.fueldensity, 2)} [L]")
    print(f"ZFW @ max fuel: {np.round(ZFW_maxfuel)} [kg]. The aircraft carries {np.round(PL_maxfuel)} [kg] of payload")
    print(f"The TOW @ ferry configuration is {np.round(Ferryweight)} [kg]")
    print("-------------------------------------------------------------------")
    rho_cr  = atm_parameters(atm_obj, h_cruise)[2]
    W = ZFW
    Range   = np.empty(0)
    n_boxes = ac_obj.n_boxes                                                            # 12
    W_PL    = np.empty(0)
    TOW     = np.empty(0)
    for i in range(len(Fuel_loads)):
        W_pl = n_boxes * ac_obj.boxweight
        W_f  = Fuel_loads[i] * ac_obj.fueldensity
        W = ac_obj.W_OE + Reserve + W_f + W_pl
        if W > ac_obj.W_TO and n_boxes > 0:
            n_boxes -= 2
            W_pl = n_boxes * ac_obj.boxweight
            W = ac_obj.W_OE + Reserve + W_f + W_pl
        TOW   = np.append(TOW, W)
        W_PL  = np.append(W_PL, W_pl)
        # R     = (ac_obj.prop_eff / ac_obj.SFC) * (CL_cr / CD_cr) * np.log(W/(W-(Fuel_loads[i]*ac_obj.fueldensity)))
        Mf_used = 0.0
        R   = 0.0
        t   = 0.0
        dt  = 1.0
        while Mf_used < W_f:
            CL_cr = 2*W*atm_obj.g/(rho_cr*V_cruise**2*ac_obj.Sw)
            CD_cr = dragpolar(ac_obj, CL_cr)
            P_br  = (1/ac_obj.prop_eff) * 1/2 * rho_cr * V_cruise**3 * ac_obj.Sw * CD_cr
            Mf_used += P_br * ac_obj.SFC
            t += dt
            R += V_cruise * dt
        print(f"The range with {int(n_boxes)} boxes and {np.round(W_f, 2)} [kg] ({np.round(Fuel_loads[i], 2)} [L]) of fuel is {np.round(R/1000)} [km] | Take-off weight: {np.round(W, 2)} [kg]")
        Range = np.append(Range, R)
    print("----------------------------------------------------------------")
    print(f"The number of boxes on the aircraft at max fuel is {n_boxes}")
    print("----------------------------------------------------------------")
    for j in range(1, int(n_boxes/2 + 1)):
        W_f   = ac_obj.fuelcapacity * ac_obj.fueldensity
        n_boxes -= 2
        W     = ac_obj.W_OE + Reserve + n_boxes * ac_obj.boxweight + W_f
        TOW   = np.append(TOW, W)
        W_PL  = np.append(W_PL, (n_boxes * ac_obj.boxweight))
        Mf_used = 0.0
        R   = 0.0
        t   = 0.0
        dt  = 1.0
        while Mf_used < W_f:
            CL_cr = 2*W*atm_obj.g/(rho_cr*V_cruise**2*ac_obj.Sw)
            CD_cr = dragpolar(ac_obj, CL_cr)
            P_br  = (1/ac_obj.prop_eff) * 1/2 * rho_cr * V_cruise**3 * ac_obj.Sw * CD_cr
            Mf_used += P_br * ac_obj.SFC
            t += dt
            R += V_cruise * dt
        # R     = (ac_obj.prop_eff / ac_obj.SFC) * (CL_cr / CD_cr) * np.log(W/(W-W_f))
        print(f"The range with {int(n_boxes)} boxes and {np.round(W_f, 2)} [kg] ({np.round(W_f/ac_obj.fueldensity, 2)} [L]) is {np.round(R/1000)} [km] | Take-off weight: {np.round(W, 2)} [kg]")
        Range = np.append(Range, R)
    Range /= 1000
    print("----------------------------------------------------------------")
    if plot:
        plt.plot(Range, W_PL, color = 'red', label = "Range-Payload")
        plt.plot(Range, TOW, color = 'blue', label = "Take-off weight - Payload")
        plt.xlabel("Range [km]")
        plt.ylabel("Weight [kg]")
        plt.legend()
        plt.show()


# -------------------------------- LANDING -----------------------------------
def LA_eom(obj, ap, atmos, max_runwayslope, max_hairport, max_headwind, max_tailwind, w_fuel, Plot=True):

    CL_all = []
    for i in range(0, 4):
        if i == 0:
            dic_constants = {'runway slope': np.arange(0, max_runwayslope, -1),
                             'airport altitude': 0, 'wing surface area': obj.Sw,
                             'weight': (obj.W_OE + w_fuel) * atmos.g,
                             'wind speed': 0, 'propeller power': obj.power * hp_to_watt,
                             'propeller efficiency': obj.eta_p}
            p, T, rho, a = atm_parameters(atmos, dic_constants['airport altitude'])

        if i == 1:
            dic_constants['runway slope'] = 0
            dic_constants['wind speed'] = np.arange(0, max_headwind)

        if i == 2:
            dic_constants['wind speed'] = np.arange(0, max_tailwind, -1)

        if i == 3:
            dic_constants['wind speed'] = 0
            dic_constants['airport altitude'] = np.arange(0, max_hairport)
            p, T, rho, a = atm_parameters(atmos, dic_constants['airport altitude'])

        S = [600.0]
        CL = 0.8
        CL_max = 2
        while max(S) <= 750:
            S_old = S.copy()
            CL += 0.01
            V_T = 1.2 * np.sqrt(
                dic_constants['weight'] / dic_constants['wing surface area'] * 2 / rho * 1 / CL_max) - dic_constants['wind speed']
            V_avg_sq = V_T ** 2 / 2
            CD = obj.CD0 + CL**2 / (np.pi * obj.A * obj.e)
            D = CD * V_avg_sq * rho / 2 * dic_constants['wing surface area']
            L = CL * V_avg_sq * rho / 2 * dic_constants['wing surface area']
            D_g = ap.mu_ground_break * (dic_constants['weight'] * np.cos(np.radians(dic_constants['runway slope'])) - L)
            # print(L/dic_constants['weight'] * np.cos(np.radians(dic_constants['runway slope'])))
            a = atmos.g / dic_constants['weight'] * (-D-D_g-dic_constants['weight']*np.sin(np.radians(dic_constants['runway slope'])))
            S = - V_T**2 / (2 * a)
            if abs(CL - (CL_max+0.01)) <= 0.0001:
                break

        S = S_old
        CL -= 0.01
        CL_all.append(CL)
        if len(S) == 1:
            print('One or more of the flight conditions; runway slope, headwind, tailwind, airport altitude, '
                  'are not achievable')
            Plot = False
            break
        else:
            CL_land = min(CL_all)
            print(CL_all, round(CL_land, 1))
            obj.FP_CL_max_land = CL_max
            obj.FP_CL_land = CL_land

    if Plot:
        figure, axis = plt.subplots(2, 2)
        if i == 0:
            axis[0, 0].plot(dic_constants['runway slope'], S)
            axis[0, 0].set_title('runway slope vs runway length')
            axis[0, 0].set_xlabel('runway slope[deg]')
            axis[0, 0].set_ylabel('runway length [m]')

        elif i == 1:
            axis[1, 0].plot(dic_constants['wind speed'], S, color='red')
            axis[1, 0].set_title('headwind vs runway length')
            axis[1, 0].set_xlabel('headwind speed [m/sec]')
            axis[1, 0].set_ylabel('runway length [m]')

        elif i == 2:
            axis[1, 1].plot(dic_constants['wind speed'], S, color='green')
            axis[1, 1].set_title('tailwind vs runway length')
            axis[1, 1].set_xlabel('tailwind speed [m/sec]')
            axis[1, 1].set_ylabel('runway length [m]')

        else:
            axis[0, 1].plot(dic_constants['airport altitude'], S, color='black')
            axis[0, 1].set_title('airport altitude vs C_L landing')
            axis[0, 1].set_xlabel('airport altitude [m]')
            axis[0, 1].set_ylabel('C_L landing [-]')

        plt.subplots_adjust(hspace=0.6)
        plt.subplots_adjust(wspace=0.5)
        plt.suptitle('Landing')
        plt.show()

    return CL_max


# LA_eom(aircraft, airfield, atm, -8, 4000, 14, -5.144, 5)

# ------------------------------------------------------------------------------
def descend(obj, atmos, V, W, P_br_max, h_descend, P_descend):
    # P_descend is the throttle setting while descending

    P_br_max *= hp_to_watt
    V *= 0.514444
    t = 0
    h = h_descend
    h_sc = 15.24  # m
    P_br_d = obj.power * P_descend * hp_to_watt

    # approach:

    throttle_setting = []
    altitude = []
    RC = np.empty(0)
    gamma_d = np.empty(0)
    while h > 0:
        p, T, rho, a = atm_parameters(atmos, h)
        # descending
        if h > h_sc:
            Pa = obj.power * P_descend * obj.eta_p * hp_to_watt * (rho / atmos.rho0) ** (3 / 4)
            C_L_descend = 2 * W / (rho * obj.Sw * V**2)
            C_D_descend = dragpolar(obj, C_L_descend)
            D_descend = C_D_descend * 1/2 * rho * V**2 * obj.Sw
            Pr = D_descend * V
            RC_current = (Pa - Pr) / W  # where P_a is less than P_r in descending
            gamma = - RC_current / V * (180/np.pi)
            RC = np.append(RC, RC_current)
            gamma_d = np.append(gamma_d, gamma)


        # Approach
        else:
            gamma = -3  # degrees
            C_L = W / obj.Sw * 2 / rho * 1 / (V ** 2)
            C_D = dragpolar(obj, C_L)
            D = C_D * 1 / 2 * rho * V ** 2 * obj.Sw

            T = W * np.sin(np.radians(gamma)) + D
            P_br_d = T * V / obj.eta_p
            P_throttle = P_br_d / P_br_max  # P_br_max is the max break power that the engine can deliver
            throttle_setting.append(P_throttle*100)

        # update constants
        dt = 0.1
        t += dt
        altitude.append(h)
        h += V*np.sin(np.radians(gamma)) * dt

        FMF = obj.SFC * P_br_d
        dWf = FMF * dt
        W -= dWf

    figure, axis = plt.subplots(2, 2)

    axis[0, 0].plot(altitude[0:len(RC)], RC)
    axis[0, 0].set_title('Altitude vs RC')
    axis[0, 0].set_xlabel('Altitude [m]')
    axis[0, 0].set_ylabel('Rate of Descend [m/s]')

    axis[0, 1].plot(altitude[0:len(RC)], gamma_d, color='purple')
    axis[0, 1].set_title('Altitude vs Descending Angle')
    axis[0, 1].set_xlabel('Altitude [m]')
    axis[0, 1].set_ylabel('Descending Angle [deg]')

    plt.subplots_adjust(wspace=0.5)
    plt.subplots_adjust(hspace=0.5)
    plt.show()

    plt.plot(throttle_setting, altitude[len(RC):], color='black')
    plt.title('Throttle Setting  in Approach vs Altitude')
    plt.xlabel('Throttle Setting [%]')
    plt.ylabel('Altitude [m]')

    plt.show()


# descend(aircraft, atm, 90, (aircraft.W_OE+10)*atm.g, 95, 500, 0.6)


def cruiseheight(distance, desired_alt):
    if 0.0 <= distance <= 10000.0:
        h_cruise = 1000 * 0.3048
    elif 10000.0 <= distance <= 50000.0:
        h_cruise = 1000 + (distance - 10000.0) * 0.1
    elif 50000.0 <= distance <= 100000.0:
        h_cruise = 5000 + (distance - 50000.0) * 0.05
    elif 100000.0 <= distance <= 200000.0:
        h_cruise = 7500 + (distance - 100000.0) * 0.025
    else:
        h_cruise = desired_alt
    h_cruise = round(h_cruise/500) * 500 * 0.3048
    return h_cruise

def climbmaneuver(ac_obj, atm_obj, h, h_cruise, W0):
    x   = 0.0
    t   = 0.0
    dt  = 1.0
    p, rho  = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
    CL_opt = np.sqrt(3*ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    W = W0
    W_F_used = 0.0
    while h < h_cruise:
        V   = np.sqrt(2*W*atm_obj.g / (rho*ac_obj.Sw*CL_opt))
        Pa  = ac_obj.power * ac_obj.prop_eff * p/atm_obj.p0 * hp_to_watt
        CD  = dragpolar(ac_obj, CL_opt) 
        Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD
        ROC = (Pa-Pr)/(W*atm_obj.g)
        gamma = ROC / V
        x   += V*np.cos(gamma)*dt
        h   += ROC * dt
        W   -= (Pa / ac_obj.prop_eff) * ac_obj.SFC * dt
        W_F_used += (Pa / ac_obj.prop_eff) * ac_obj.SFC * dt
        t   += dt
        p, rho  = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
    return t, x, W_F_used, h

def cruisecalc(ac_obj, atm_obj, h_cruise, distance, W0, V_cruise = None):
    n = 0
    cruiseNAT = False
    x   = 0.0
    t   = 0.0
    dt  = 1.0
    p, rho  = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
    W = W0
    W_F_used = 0.0
    if V_cruise == None:
            CL  = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
            CD  = dragpolar(ac_obj, CL)
    while x <= distance:
        if V_cruise == None:
            V   = np.sqrt(2*W*atm_obj.g/(rho*ac_obj.Sw*CL))
        else:
            CL  = 2*W*atm_obj.g/(rho*V_cruise**2*ac_obj.Sw)
            CD  = dragpolar(ac_obj, CL)
            V   = V_cruise
        x += V * dt
        Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD
        Pbr = Pr/ac_obj.prop_eff
        # if n%1000 == 0:
        #     print(f"Break horse power: {Pbr}")
        if Pbr > ac_obj.power * ac_obj.prop_eff * p/atm_obj.p0 * hp_to_watt:
            print("This cruise speed is not obtainable (Pr > Pa)")
            cruiseNAT = True
            break
        FF  = Pbr * ac_obj.SFC
        W_F_used += FF * dt
        W -= FF * dt
        t += dt
        n += 1
    return t, W_F_used, cruiseNAT

# def cruiseperf_varying(ac_obj, atm_obj):
#     W_F_used = np.empty(0)
#     h_cruise = np.arange(5000, 15000, 500) * 0.3048
#     V_cruise = np.arange(45, 70, 5)
#     x_cruise = np.arange(10000, 111000, 1000)
#     for j in range(3):
#         V = 50 + 5*j
#         for i in range(len(h_cruise)):
#             W_F_used_s = cruisecalc(ac_obj, atm_obj, h_cruise[i], 100000, 710, V)
#             W_F_used = np.append(W_F_used, W_F_used_s)
#         plt.plot(h_cruise, W_F_used[20*j:20*j + 20], label=f"Cruise speed: {V} [m/s]")
#     plt.xlabel("h_cruise [m]")
#     plt.ylabel("W_F_used [kg]")
#     plt.title("W_F_used versus h_cruise")
#     plt.legend()
#     plt.show()

# cruiseperf_varying(aircraft, atm)


def descentmaneuver(ac_obj, atm_obj, h_cruise, h_stop, W0):
    CL = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    CD = dragpolar(ac_obj, CL)
    LD_max = CL / CD
    # print(f"Maximum lift to drag ratio: {LD_max} [-]")
    gamma_LD_max = np.arctan(1/(LD_max))
    h = h_cruise
    p, rho  = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
    W = W0
    x = 0.0
    t = 0.0
    dt = 1.0
    FF_idle = 3 / 3600
    W_F_used = 0.0
    while h > h_stop:
        V   = np.sqrt(2*W*atm_obj.g/(rho*ac_obj.Sw*CL))
        D   = 1/2 * rho * V**2 * ac_obj.Sw * CD
        gamma = - D / (W*atm_obj.g)
        ROD = V*gamma
        h   += ROD * dt
        W   -= FF_idle * dt
        W_F_used += FF_idle * dt
        x   += V*np.cos(gamma)*dt
        p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
        t += dt
    return t, x, W_F_used, h

def loiter(ac_obj, atm_obj, h_loiter, t_loiter, W0, standardrate, goback = True, result = False):
    h_loiter *= 0.3048
    CL_opt  = np.sqrt(3*ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    CD_opt  = dragpolar(ac_obj, CL_opt)
    p, rho  = atm_parameters(atm_obj, h_loiter)[0], atm_parameters(atm_obj, h_loiter)[2]
    heading_t0  = 0.0
    heading     = 0.0
    standardrate *= (3 * np.pi / 180)
    x   = 0.0
    t   = 0.0
    dt  = 1.0
    W   = W0
    W_F_used = 0.0
    patterns = 0
    while t <= t_loiter:
        V   = np.sqrt(2*W*atm_obj.g/(rho*ac_obj.Sw*CL_opt))
        R_turn    = V / standardrate
        bankangle = np.arctan(V**2/(atm_obj.g * R_turn))
        n   = 1/np.cos(bankangle)
        Pr  = n * 1/2 * rho * V**3 * ac_obj.Sw * CD_opt
        Pbr = Pr/ac_obj.prop_eff
        W   -= Pbr * ac_obj.SFC * dt
        W_F_used += Pbr * ac_obj.SFC * dt
        heading += standardrate * dt
        if heading >= 2 * np.pi:
            heading -= 2 * np.pi
            patterns += 1
        x   += V * dt
        t   += dt
    heading -= 2 * np.pi
    print(f"heading at t = t_loiter: {round(heading*180/np.pi, 2)}")
    if goback == True:
        target_heading = heading_t0 + np.pi
    else:
        target_heading = heading_t0
    print(f"target heading: {round(target_heading*180/np.pi, 2)}")
    while heading < target_heading:
        V   = np.sqrt(2*W*atm_obj.g/(rho*ac_obj.Sw*CL_opt))
        R_turn    = V / standardrate
        bankangle = np.arctan(V**2/(atm_obj.g * R_turn))
        n   = 1/np.cos(bankangle)
        Pr  = n * 1/2 * rho * V**3 * ac_obj.Sw * CD_opt
        Pbr = Pr/ac_obj.prop_eff
        W   -= Pbr * ac_obj.SFC * dt
        W_F_used += Pbr * ac_obj.SFC * dt
        heading += standardrate * dt
        x   += V * dt
        t   += dt
    if not goback:
        patterns += 1
    else:
        patterns += 1/2
    if result:
        print(f"======================== Loiter Summary ========================")
        print(f"Loiter time: {t_loiter} [sec] | Standard turn rate: {np.round(standardrate*180/np.pi, 2)} [deg/s]")
        print(f"Fuel used during loiter: {np.round(W_F_used, 2)} | Number of patterns completed: {patterns}")
        print(f"================================================================")
    return t, W_F_used
# loiter(aircraft, atm, 7500, 1200, 700, 1, goback=True, result=True)

def fuelusesortie(ac_obj, atm_obj, n_boxes, n_drops, h_cruise, W_F, V_cruise = None, Range = None, dropregion = None, Summary = False, plot = False, savefig = False):
    cruiseNAT = False
    flight_profile = []
    # Define a number of arrays used for plotting
    h_array = np.empty(0)
    x_array = np.empty(0)
    t_array = np.empty(0)
    W_F_used_array = np.empty(0)
    # ===========================================================================
    starttime = time.time()                                                     #
    W_F_used = 0.0                                                              #
    x   = 0.0                                                                   #
    x_plot = 0.0                                                                #
    t   = 0.0                                                                   #
    dt  = 1.0                                                                   #
    h   = 15.0                                              #after screenheight #
    # Add start values to arrays                                                #
    x_array = np.append(x_array, x)                                             #
    h_array = np.append(h_array, h)                                             #
    t_array = np.append(t_array, t)                                             #
    W_F_used_array = np.append(W_F_used_array, W_F_used)                        #
    # ===========================================================================
    CL = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    CD = dragpolar(ac_obj, CL)
    LD_max = CL / CD
    if n_drops != 0:
        boxesperdrop = n_boxes / n_drops
        if boxesperdrop % 2 != 0:
            print(f"The amount of boxes ({n_boxes}) and the amount of drops ({n_drops}) are not compatible.")
            print(f"The boxes need to be dropped in multiples of two to avoid a lateral C.O.G. shift")
    # Enter h_cruise in feet
    if Range == None:
        Range = ac_obj.R / 2
    else:
        Range *= 1000
    if V_cruise == None:
        CL_cr_opt = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    else:
        V_cruise = V_cruise
    W_TO = ac_obj.W_OE + W_F + n_boxes * ac_obj.boxweight
    W_F_0 = W_TO
    W = W_TO
    if Summary:
        print("=====================================================")
        print(f"Take-off weight: {np.round(W_TO, 2)} [kg] | OEW: {np.round(ac_obj.W_OE, 2)} [kg] | Fuelweight: {np.round(W_F, 2)} [kg] / {np.round(W_F/ac_obj.fueldensity, 2)} [L] | Payload: {np.round(n_boxes, 2)} boxes / {np.round(n_boxes*ac_obj.boxweight, 2)} [kg]")
    W_a_TO = W_TO * ac_obj.W1W_TO * ac_obj.W2W1 * ac_obj.W3W2
    W = W_a_TO
    W_F_used = W_TO - W_a_TO
    print(f"W_F_used after take off: {round(W_F_used, 2)} | fraction: {round(1- W_F_used/W_F_0, 3)}")
    W_F -= W_F_used
    # Now comes climb to projected cruise
    if n_boxes == 0:                            # flight without payload, surveillance ??
        cruise_alt = cruiseheight(Range, h_cruise)
        t_climb, x_climb, W_F_used_climb, h = climbmaneuver(ac_obj, atm_obj, h, cruise_alt, W)
        t += t_climb
        x += x_climb
        x_array = np.append(x_array, x)
        h_array = np.append(h_array, h)
        t_array = np.append(h_array, t)
        W_F_used += W_F_used_climb
        W_F_used_array = np.append(W_F_used_array, W_F_used)
        W -= W_F_used_climb
        # Assume loiter altitude 5000 ft
        dtt_remaining = Range - x_climb
        descent_dist = (cruise_alt-(5000*0.3048)) * LD_max
        cruise_dist = dtt_remaining - descent_dist
        t_cruise, W_F_used_cruise, cruiseNAT = cruisecalc(ac_obj, atm_obj, cruise_alt, cruise_dist, W, V_cruise)
        t += t_cruise
        # print(f"Cruise time: {round(t_cruise, 2)} [sec]")
        x += cruise_dist
        x_array = np.append(x_array, x)
        h_array = np.append(h_array, h)
        t_array = np.append(h_array, t)
        x_between_drop = cruise_dist
        W_F_used += W_F_used_cruise
        W_F_used_array = np.append(W_F_used_array, W_F_used)
        # print(f"W_F_used after cruise {i}: {W_F_used_cruise}")
        W -= W_F_used_cruise
        t_descent, x_descent, W_F_used_descent, h = descentmaneuver(ac_obj, atm_obj, cruise_alt, 15.0, W)
        # print(f"Horizontal distance descent: {x_descent} [m]")
        t += t_descent
        x += x_descent
        x_between_drop += x_descent
        W_F_used += W_F_used_descent
        x_array = np.append(x_array, x)
        h_array = np.append(h_array, h)
        t_array = np.append(h_array, t)
        W_F_used_array = np.append(W_F_used_array, W_F_used)
        # print(f"W_F_used after descent {i}: {W_F_used_descent} | time descent {i}: {t_descent}")
        W -= W_F_used_descent
        t_loiter, W_F_used_loiter = loiter(ac_obj, atm_obj, 5000, 1200, W, 1, goback=True)
        t += t_loiter
        W_F_used += W_F_used_loiter
    elif dropregion == None:                      # either 1 drop or n evenly spaced drops
        target_dist = Range / n_drops
        cruise_alt  = cruiseheight(target_dist, h_cruise)
        print(f"The cruising altitude for this sortie is {round(cruise_alt/0.3048)} [ft] (in between all drops)")
        for i in range(n_drops):
            x_between_drop = 0.0
            # Climb part
            t_climb, x_climb, W_F_used_climb, h = climbmaneuver(ac_obj, atm_obj, h, cruise_alt, W)
            # print(f"Horizontal distance climb: {x_climb} [m]")
            t += t_climb
            x += x_climb
            x_array = np.append(x_array, x)
            h_array = np.append(h_array, h)
            t_array = np.append(h_array, t)
            x_between_drop += x_climb
            W_F_used += W_F_used_climb
            W_F_used_array = np.append(W_F_used_array, W_F_used)
            # print(f"W_F_used after climb {i}: {round(W_F_used_climb, 2)} [kg] in {t_climb} [sec] | fraction: {round(1 - W_F_used_climb/W_F_0, 3)}")
            W -= W_F_used_climb
            # print(f"Start weight cruise: {np.round(W, 2)} [kg]")
            # Cruise part - first calculate cruise distance
            dtt_remaining = target_dist - x_between_drop
            descent_dist = (cruise_alt-15.0) * LD_max
            cruise_dist = dtt_remaining - descent_dist
            # print(f"Climb distance: {x_climb} [m] | Cruise distance: {cruise_dist} [m] | Descent distance: {descent_dist} [m]")
            t_cruise, W_F_used_cruise, cruiseNAT = cruisecalc(ac_obj, atm_obj, cruise_alt, cruise_dist, W, V_cruise)
            # print(f"cruise time: {t_cruise}")
            t += t_cruise
            x += cruise_dist
            x_array = np.append(x_array, x)
            h_array = np.append(h_array, h)
            t_array = np.append(h_array, t)
            x_between_drop = cruise_dist
            W_F_used += W_F_used_cruise
            W_F_used_array = np.append(W_F_used_array, W_F_used)
            # print(f"Cruise time: {round(t_cruise, 2)} [sec]")
            # print(f"W_F_used after cruise {i}: {round(W_F_used_cruise, 2)} | fraction: {round(1 - W_F_used_cruise/W_F_0, 3)}")
            W -= W_F_used_cruise
            # Descent part
            t_descent, x_descent, W_F_used_descent, h = descentmaneuver(ac_obj, atm_obj, cruise_alt, 15.0, W)
            # print(f"Horizontal distance descent: {x_descent} [m]")
            t += t_descent
            x += x_descent
            x_between_drop += x_descent
            W_F_used += W_F_used_descent
            x_array = np.append(x_array, x)
            h_array = np.append(h_array, h)
            t_array = np.append(h_array, t)
            W_F_used_array = np.append(W_F_used_array, W_F_used)
            # print(f"W_F_used after descent {i}: {round(W_F_used_descent, 2)} | time descent {i}: {t_descent} | fraction: {round(1- W_F_used_descent/W_F_0, 3)}")
            W -= W_F_used_descent
            # After each descent, a certain amount of boxes are dropped
            W -= boxesperdrop * ac_obj.boxweight
            if i == 0:
                ttfd = t
                ac_obj.TTFD_s = t
                flight_profile.append(t)
    else:
        target_dist = Range - dropregion*1000                           # distance to first target
        cruise_alt_tft = cruiseheight(target_dist, h_cruise)            # tft = to first target
        # To first target
        # Climb part
        t_climb, x_climb, W_F_used_climb, h = climbmaneuver(ac_obj, atm_obj, h, cruise_alt_tft, W)
        t += t_climb
        x += x_climb
        W_F_used += W_F_used_climb
        x_array = np.append(x_array, x)
        h_array = np.append(h_array, h)
        t_array = np.append(h_array, t)
        W_F_used_array = np.append(W_F_used_array, W_F_used)
        W -= W_F_used_climb
        # Cruise part - first calculate cruise distance
        dtt_remaining = target_dist - x
        descent_dist = (cruise_alt_tft - 15.0) * LD_max
        cruise_dist = dtt_remaining - descent_dist
        # print(f"cruise dist to first target: {cruise_dist} [m]")
        t_cruise, W_F_used_cruise, cruiseNAT = cruisecalc(ac_obj, atm_obj, cruise_alt_tft, cruise_dist, W, V_cruise)
        t += t_cruise
        x += cruise_dist
        W_F_used += W_F_used_cruise
        W -= W_F_used_cruise
        x_array = np.append(x_array, x)
        h_array = np.append(h_array, h)
        t_array = np.append(h_array, t)
        W_F_used_array = np.append(W_F_used_array, W_F_used)
        # Descent part
        t_descent, x_descent, W_F_used_descent, h = descentmaneuver(ac_obj, atm_obj, cruise_alt_tft, 15.0, W)
        t += t_descent
        x += x_descent
        W_F_used += W_F_used_descent
        x_array = np.append(x_array, x)
        h_array = np.append(h_array, h)
        t_array = np.append(h_array, t)
        W_F_used_array = np.append(W_F_used_array, W_F_used)
        W -= W_F_used_descent
        W -= boxesperdrop * ac_obj.boxweight
        ttfd = t
        ac_obj.TTFD_s = t
        flight_profile.append(t)
        # ============================================================================================================================
        # Drops in dropregion
        interdropdist = (dropregion*1000) / (n_drops - 1)
        cruise_alt_int_drop = cruiseheight(interdropdist, h_cruise)
        for i in range(n_drops - 1):
            x_indropregion = 0.0
            t_cl_indrop, x_cl_indrop, W_F_used_cl_indrop, h = climbmaneuver(ac_obj, atm_obj, h, cruise_alt_int_drop, W)
            t += t_cl_indrop
            x += x_cl_indrop
            x_indropregion += x_cl_indrop
            # print(f"Climb distance in dropzone {i}: {x_cl_indrop} [m]")
            W_F_used += W_F_used_cl_indrop
            x_array = np.append(x_array, x)
            h_array = np.append(h_array, h)
            t_array = np.append(h_array, t)
            W_F_used_array = np.append(W_F_used_array, W_F_used)
            W -= W_F_used_cl_indrop
            # Cruise
            dtnd = interdropdist - x_indropregion
            d_des = LD_max * (cruise_alt_int_drop-15.0)
            cruise_dist_int_drop = dtnd - d_des
            print(f"Cruise dist in dropzone: {cruise_dist_int_drop} [m]")
            # print(f'Inter-drop cruise distance {i}: {cruise_dist_int_drop} [m]')
            t_cr_indrop, W_F_used_cr_indrop, cruiseNAT = cruisecalc(ac_obj, atm_obj, cruise_alt_int_drop, cruise_dist_int_drop, W, V_cruise)
            t += t_cr_indrop
            x += cruise_dist_int_drop
            x_indropregion += cruise_dist_int_drop
            W_F_used += W_F_used_cr_indrop
            x_array = np.append(x_array, x)
            h_array = np.append(h_array, h)
            t_array = np.append(h_array, t)
            W_F_used_array = np.append(W_F_used_array, W_F_used)
            W -= W_F_used_cr_indrop
            # Descent
            t_des_indrop, x_des_indrop, W_F_used_des_indrop, h = descentmaneuver(ac_obj, atm_obj, cruise_alt_int_drop, 15.0, W)
            t += t_des_indrop
            x += x_des_indrop
            # print(f"descent distance in dropzone {i}: {x_des_indrop} [m]")
            W_F_used += W_F_used_des_indrop
            W -= W_F_used_des_indrop
            x_array = np.append(x_array, x)
            h_array = np.append(h_array, h)
            t_array = np.append(h_array, t)
            W_F_used_array = np.append(W_F_used_array, W_F_used)
            # Drop box
            W -= boxesperdrop * ac_obj.boxweight
    x_plot = x
    # Return to base
    x_return = 0.0
    cruise_alt_RTB = cruiseheight(Range, 15000)
    t_cl_RTB, x_cl_RTB, W_F_used_cl_RTB, h = climbmaneuver(ac_obj, atm_obj, h, cruise_alt_RTB, W)
    t += t_cl_RTB
    x += x_cl_RTB
    x_return += x_cl_RTB
    x_plot -= x_cl_RTB
    W_F_used += W_F_used_cl_RTB
    x_array = np.append(x_array, x_plot)
    h_array = np.append(h_array, h)
    t_array = np.append(h_array, t)
    W_F_used_array = np.append(W_F_used_array, W_F_used)
    print(f"W_F_used for climb RTB: {round(W_F_used_cl_RTB, 2)} [kg] | fraction: {round(1 - W_F_used_cl_RTB/W_F_0, 2)}")
    W -= W_F_used_cl_RTB
    # Cruise return to base
    dtb_remaining = Range - x_return
    desc_dist_RTB = LD_max * cruise_alt_RTB
    cruise_dist_RTB = dtb_remaining - desc_dist_RTB
    # print(f"Climb distance: {x_cl_RTB} [m] | Cruise distance: {cruise_dist_RTB} [m] | Descent distance: {desc_dist_RTB} [m]")
    t_cr_RTB, W_F_used_cr_RTB, cruiseNAT = cruisecalc(ac_obj, atm_obj, cruise_alt_RTB, cruise_dist_RTB, W, V_cruise)
    t += t_cr_RTB
    x += cruise_dist_RTB
    x_plot -= cruise_dist_RTB
    W_F_used += W_F_used_cr_RTB
    x_array = np.append(x_array, x_plot)
    h_array = np.append(h_array, h)
    t_array = np.append(h_array, t)
    W_F_used_array = np.append(W_F_used_array, W_F_used)
    print(f"W_F_used for cruise RTB: {round(W_F_used_cr_RTB, 2)} [kg] | fraction: {round(1 - W_F_used_cr_RTB/W_F_0, 3)}")
    W -= W_F_used_cr_RTB
    # Descent
    t_des_RTB, x_des_RTB, W_F_used_des_RTB, h = descentmaneuver(ac_obj, atm_obj, cruise_alt_RTB, 15.0, W)
    t += t_des_RTB
    x += x_des_RTB
    x_plot -= x_des_RTB
    W_F_used += W_F_used_des_RTB
    x_array = np.append(x_array, x_plot)
    h_array = np.append(h_array, h)
    t_array = np.append(h_array, t)
    W_F_used_array = np.append(W_F_used_array, W_F_used)
    # print(f"W_F_used for descent RTB: {round(W_F_used_des_RTB, 2)} [kg] | fraction: {round(1 - W_F_used_des_RTB/W_F_0, 3)}")
    W -= W_F_used_des_RTB
    W_beforelanding = W
    # Landing, taxi, shutdown using fuel fractions
    W *= ac_obj.WfinalW10
    print(f"W_F_used landing, taxi, shutdown: {round(W_beforelanding - W, 2)} | fraction: {round(W/W_beforelanding, 3)}")
    W_F_used += (W_beforelanding - W)
    if W_F_used > W_F:
        print("The fuel used is higher than the fuel loaded, this is not possible")
    # Calculate the fuel use per kg payload per km range (that means half the distance flown)
    F_kg_km = np.round(W_F_used/(n_boxes*ac_obj.boxweight*Range/1000), 5)
    print(f"The fuel used to transport one kg of payload over one km is {F_kg_km} [1/km]")
    endtime = time.time()
    if Summary:
        if cruiseNAT == False:
            print(f"=================================== Sortie Summary ===================================")
            print(f"Payload: {n_boxes} boxes ({np.round(n_boxes*ac_obj.boxweight, 2)} [kg]) | {n_drops} dropping maneuvers")
            print(f"Range: {np.round(Range/1000, 2)} [km] | Distance flown: {np.round(x/1000, 2)} [km] | Flight time: {np.round(t/3600, 2)} [hrs]")
            if n_drops != 0:
                print(f"Time to first drop: {ttfd} [sec] / {np.round(t/3600, 2)} [hrs]")
            print(f"Fuel used: {np.round(W_F_used, 2)} [kg]")
            print(f"======================================================================================")
            print(f"This simulation took {np.round((endtime-starttime), 2)} [s]")
            print(f"======================================================================================")
    if plot:
        plt.plot(x_array, h_array)
        plt.xlabel("Horizontal distance [m]")
        plt.ylabel("Altitude [m]")
        if savefig:
            if dropregion != None:
                filepath = "C:\\Users\\ties\\Downloads\\flightprofile-"+ str(Range) + str(n_drops) + str(dropregion)+ ".png"
                plt.savefig(filepath)
            else:
                filepath = "C:\\Users\\ties\\Downloads\\flightprofile-"+ str(Range) + str(n_drops) + ".png"
                plt.savefig(filepath)
        plt.show()

    if len(flight_profile)==0: flight_profile.append(0)
    flight_profile.append(t)                # total sortie time
    flight_profile.append(W_F_used)         # fuel burnt

    ac_obj.T_s = t
    #ac_obj.Wf = W_F_used

    return flight_profile

# fuelusesortie(aircraft, atm, 12, 1, 10000, 45, 54.012, Range=250, Summary=True, plot=True, savefig=False)
