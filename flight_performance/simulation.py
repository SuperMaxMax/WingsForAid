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

aircraft = UAV("aircraft")
airfield = airport("Sudan")
atm      = atmosphere()
hp_to_watt = 745.699872

def takeoffweight(obj, W_F):
    ZFW = obj.W_OE + obj.n_boxes*obj.boxweight
    TOW = ZFW + W_F
    return TOW

def atm_parameters(obj, h):
    T    = atm.T0 + atm.lambd * h
    rho  = atm.rho0*np.power((T/obj.T0), (-((atm.g / (atm.lambd * atm.R))+1)))
    p    = atm.p0*np.power((T/obj.T0), (-(atm.g / (atm.lambd * atm.R))))
    a    = np.sqrt(atm.gamma*obj.R*T)
    return p, T, rho, a

def dragpolar(ac_obj, CL):
    CD = ac_obj.CD0 + CL**2/(np.pi*ac_obj.A*ac_obj.e)
    return CD

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

def flightceiling(ac_obj, atm_obj, W_F, plot=True):
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

# # ---------------- Assumptions for take-off equations of motion -----------------
# # Wind is included by take it into account in the speed: V_eff = V - V_wind
# # Runway slope is not zero
# # The present of rain is taken into account in the friction coefficient with the ground, mu
# # delta_rw = runway slope
# # D_g = force due to the ground friction, with
# # Thrust and lift are taken as average values


def TO_eom(obj, ap, atmos, Plot=True):

    figure, axis = plt.subplots(2, 2)

    for i in range(0, 4):
        if i == 0:
            dic_constants = {'runway slope': np.arange(0, 10),
                             'airport altitude': 0, 'wing surface area': 11,
                             'weight': takeoffweight(aircraft, 200) * atm.g,
                             'wind speed': 0, 'propeller power': aircraft.power * hp_to_watt,
                             'propeller efficiency': aircraft.eta_p}
            p, T, rho, a = atm_parameters(obj, dic_constants['airport altitude'])

        if i == 1:
            dic_constants['runway slope'] = 0
            dic_constants['wind speed'] = np.arange(0, 10)
            print(dic_constants, rho)

        if i == 2:
            dic_constants['wind speed'] = np.arange(0, -10, -1)

        if i == 3:
            dic_constants['wind speed'] = 0
            dic_constants['airport altitude'] = np.arange(0, 500)
            p, T, rho, a = atm_parameters(obj, dic_constants['airport altitude'])

        V_avg_sq = 0.55125 * (np.sqrt(dic_constants['weight']/dic_constants['wing surface area'] * 2/rho * 1/obj.CL_max_TO) -
                              dic_constants['wind speed']) ** 2

        A = - dic_constants['wing surface area'] / (np.pi * obj.A * obj.e) * V_avg_sq * rho/2 * atmos.g / dic_constants['weight']
        B = ap.mu_ground * dic_constants['wing surface area'] * V_avg_sq * rho/2 * atmos.g / dic_constants['weight']
        C = (dic_constants['propeller power'] * dic_constants['propeller efficiency'] / np.sqrt(V_avg_sq) - ap.mu_ground *
             dic_constants['weight'] *np.cos(np.radians(dic_constants['runway slope'])) - obj.CD0 * rho/2 * V_avg_sq
             * dic_constants['wing surface area'] - dic_constants['weight']*np.sin(np.radians(dic_constants['runway slope']))) \
             * atmos.g / dic_constants['weight'] - V_avg_sq/750

        sqrt = B**2 - 4*A*C
        CL_TO_1 = (-B + sqrt) / (2*A)
        CL_TO_2 = (-B - sqrt) / (2*A)

        if Plot:
            if i == 0:
                print(i)
                axis[0, 0].plot(dic_constants['runway slope'], CL_TO_1)
                axis[0, 0].plot(dic_constants['runway slope'], CL_TO_2)
                axis[0, 0].set_title('runway slope vs C_L take-off')
                axis[0, 0].set_xlabel('runway slope[deg]')
                axis[0, 0].set_ylabel('C_L take-off [-]')

            elif i == 1:
                print(i)
                axis[1, 0].plot(dic_constants['wind speed'], CL_TO_1, color='red')
                axis[1, 0].plot(dic_constants['wind speed'], CL_TO_2, color='red')
                axis[1, 0].set_title('headwind vs C_L take-off')
                axis[1, 0].set_xlabel('headwind speed [m/sec]')
                axis[1, 0].set_ylabel('C_L take-off [-]')

            elif i == 2:
                print(i)
                axis[1, 1].plot(dic_constants['wind speed'], CL_TO_1, color='green')
                axis[1, 1].plot(dic_constants['wind speed'], CL_TO_2, color='green')
                axis[1, 1].set_title('tailwind vs C_L take-off')
                axis[1, 1].set_xlabel('tailwind speed [m/sec]')
                axis[1, 1].set_ylabel('C_L take-off [-]')

            else:
                print(i)
                axis[0, 1].plot(dic_constants['airport altitude'], CL_TO_1, color='black')
                axis[0, 1].plot(dic_constants['airport altitude'], CL_TO_2, color='black')
                axis[0, 1].set_title('airport altitude vs C_L take-off')
                axis[0, 1].set_xlabel('airport altitude [m]')
                axis[0, 1].set_ylabel('C_L take-off [-]')

    plt.subplots_adjust(hspace=0.6)
    plt.subplots_adjust(wspace=0.5)
    plt.suptitle('Take-Off')
    plt.show()

    return CL_TO_1, CL_TO_2

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
turnperformance(aircraft, atm)

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
    while r_it < R:
        CL_cr   = 2*W/(rho_cr*V_cruise**2*ac_obj.Sw)
        CD_cr   = dragpolar(ac_obj, CL_cr)
        r_it    += (V_cruise * dt)
        t       += dt
        D       = 1/2 * rho_cr * V_cruise**2 * ac_obj.Sw * CD_cr
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
    return None
# cruiseperformance(aircraft, atm)

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
def LA_eom(obj, ap, atmos, Plot=True):

    figure, axis = plt.subplots(2, 2)

    for i in range(0, 4):
        if i == 0:
            dic_constants = {'runway slope': np.arange(0, 10),
                             'airport altitude': 0, 'wing surface area': 11,
                             'weight': aircraft.W_OE * atm.g,
                             'wind speed': 0, 'propeller power': aircraft.power * hp_to_watt,
                             'propeller efficiency': aircraft.eta_p}
            p, T, rho, a = atm_parameters(obj, dic_constants['airport altitude'])

        if i == 1:
            print('1')
            dic_constants['runway slope'] = 0
            dic_constants['wind speed'] = np.arange(0, 10)

        if i == 2:
            print('2')
            dic_constants['wind speed'] = np.arange(0, -10, -1)

        if i == 3:
            print('3')
            dic_constants['wind speed'] = 0
            dic_constants['airport altitude'] = np.arange(0, 500)
            p, T, rho, a = atm_parameters(obj, dic_constants['airport altitude'])

        V_avg_sq = 0.72 * (np.sqrt(dic_constants['weight'] / dic_constants['wing surface area'] * 2 / rho * 1 / obj.CL_max_land) -
                              dic_constants['wind speed']) ** 2
        A = - dic_constants['wing surface area'] / (np.pi * obj.A * obj.e) * V_avg_sq * rho / 2 * atmos.g / dic_constants['weight']
        B = ap.mu_ground * dic_constants['wing surface area'] * V_avg_sq * rho / 2 * atmos.g / dic_constants['weight']
        C = (800- ap.mu_ground * dic_constants['weight'] * np.cos(np.radians(dic_constants['runway slope'])) - obj.CD0 * rho / 2 * V_avg_sq
             * dic_constants['wing surface area'] - dic_constants['weight'] * np.sin(np.radians(dic_constants['runway slope']))) \
            * atmos.g / dic_constants['weight'] + V_avg_sq / 750

        sqrt = B ** 2 - 4 * A * C
        CL_LA_1 = (-B + sqrt) / (2 * A)
        CL_LA_2 = (-B - sqrt) / (2 * A)

        if Plot:
            if i == 0:
                axis[0, 0].plot(dic_constants['runway slope'], CL_LA_1)
                axis[0, 0].plot(dic_constants['runway slope'], CL_LA_2)
                axis[0, 0].set_title('runway slope vs C_L landing')
                axis[0, 0].set_xlabel('runway slope[deg]')
                axis[0, 0].set_ylabel('C_L landing [-]')

            elif i == 1:
                axis[1, 0].plot(dic_constants['wind speed'], CL_LA_1, color='red')
                axis[1, 0].plot(dic_constants['wind speed'], CL_LA_2, color='red')
                axis[1, 0].set_title('headwind vs C_L landing')
                axis[1, 0].set_xlabel('headwind speed [m/sec]')
                axis[1, 0].set_ylabel('C_L landing [-]')

            elif i == 2:
                axis[1, 1].plot(dic_constants['wind speed'], CL_LA_1, color='green')
                axis[1, 1].plot(dic_constants['wind speed'], CL_LA_2, color='green')
                axis[1, 1].set_title('tailwind vs C_L landing')
                axis[1, 1].set_xlabel('tailwind speed [m/sec]')
                axis[1, 1].set_ylabel('C_L landing [-]')

            else:
                axis[0, 1].plot(dic_constants['airport altitude'], CL_LA_1, color='black')
                axis[0, 1].plot(dic_constants['airport altitude'], CL_LA_2, color='black')
                axis[0, 1].set_title('airport altitude vs C_L landing')
                axis[0, 1].set_xlabel('airport altitude [m]')
                axis[0, 1].set_ylabel('C_L landing [-]')

    plt.subplots_adjust(hspace=0.6)
    plt.subplots_adjust(wspace=0.5)
    plt.suptitle('Landing')
    plt.show()

    return CL_LA_1, CL_LA_2


# ------------------------------------------------------------------------------
def descend(obj, atmos, V, W, P_br_max, h_descend, P_descend):
    # P_descend is the throttle setting while descending

    P_br_max *= hp_to_watt
    V *= 0.514444
    t = 0
    h = h_descend
    h_sc = 15.24  # m
    P_br_d = obj.power * P_descend * hp_to_watt

    # beginning of descend before approach: (same as climb but gamma is negative)
    # P
    # RC = (P_a - P_r) / W  # where P_a is less than P_r in descending


    # approach:

    throttle_setting = []
    altitude = []
    RC = np.empty(0)
    gamma_d = np.empty(0)
    while h > 0:
        p, T, rho, a = atm_parameters(obj, h)
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
        dt = 0.01
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

    axis[1, 0].plot(throttle_setting, altitude[len(RC):], color='black')
    axis[1, 0].set_title('Throttle Setting  in Approach vs Altitude')
    axis[1, 0].set_xlabel('Throttle Setting [%]')
    axis[1, 0].set_ylabel('Altitude [m]')

    plt.subplots_adjust(hspace=0.6)
    plt.subplots_adjust(wspace=0.5)
    plt.show()


# descend(aircraft, atm, 90, (aircraft.W_OE+100)*atm.g, 95, 500, 0.6)
