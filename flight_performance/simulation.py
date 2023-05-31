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

def takeoffweight(obj, W_F):
    ZFW = obj.W_OE + obj.n_boxes*obj.boxweight
    TOW = ZFW + W_F
    return TOW

def atm_parameters(obj, h):
    T    = (atm.T0 +15) + atm.lambd * h
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
    W = takeoffweight(ac_obj, W_F)*atm_obj.g
    atm_parameters_vectorized = np.vectorize(lambda h: atm_parameters(atm_obj, h))
    alt_range = np.arange(0, ac_obj.th_ceil, 0.5)
    atm_obj.p, atm_obj.T, atm_obj.rho, atm_obj.a = atm_parameters_vectorized(alt_range)
    stall_limit = np.empty(0)
    thr_lim_lo = np.empty(0)
    thr_lim_hi  = np.empty(0)
    for i in range(len(alt_range)):
        rho = atm_obj.rho[i]
        V_s = np.sqrt(2*W/(rho*ac_obj.Sw*ac_obj.CL_max_clean))
        P_a = ac_obj.power * ac_obj.prop_eff * 745.699872 * (rho/atm_obj.rho0)**(3/4)
        V_max = np.arange(110, 200, 0.5)*0.5144
        V_min = np.arange(40, 130, 0.5)*0.5144
        CL_hi = 2*W/(rho*ac_obj.Sw*V_max**2)
        CL_lo = 2*W/(rho*ac_obj.Sw*V_min**2)
        CD_hi = ac_obj.CD0 + CL_hi**2/(np.pi*ac_obj.A*ac_obj.e)
        CD_lo = ac_obj.CD0 + CL_lo**2/(np.pi*ac_obj.A*ac_obj.e)
        Pr_hi = 1/2 * rho * ac_obj.Sw * V_max**3 * CD_hi
        Pr_lo = 1/2 * rho * ac_obj.Sw * V_min**3 * CD_lo
        P_a_compare = np.full(len(V_max), P_a)
        PaPr_hi = P_a_compare - Pr_hi
        PaPr_lo = P_a_compare - Pr_lo
        V_max= V_max[np.argmin(np.abs(PaPr_hi))]
        V_min= V_min[np.argmin(np.abs(PaPr_lo))]
        stall_limit = np.append(stall_limit, V_s)
        thr_lim_hi = np.append(thr_lim_hi, V_max)
        thr_lim_lo = np.append(thr_lim_lo, V_min)
    if plot:
        plt.plot(thr_lim_hi, alt_range, color='green')
        plt.plot(stall_limit, alt_range, color='green')
        #plt.plot(thr_lim_lo, alt_range, color = 'green')
        plt.xlabel("Airspeed [m/s]")
        plt.ylabel("Altitude [m]")
        plt.show()
    return


# # ---------------- Assumptions for take-off equations of motion -----------------
# # Wind is included by take it into account in the speed: V_eff = V - V_wind
# # Runway slope is not zero
# # The present of rain is taken into account in the friction coefficient with the ground, mu
# # delta_rw = runway slope
# # D_g = force due to the ground friction, with
# # Thrust and lift are taken as average values


def TO_eom(obj, ap, atmos, constants):

    p, T, rho, a = atm_parameters(obj, constants['airport altitude'])

    V_avg_sq = 0.55125 * (np.sqrt(constants['weight']/constants['wing surface area'] * 2/rho * 1/obj.CL_max_TO) -
                          constants['wind speed']) ** 2
    
    A = - constants['wing surface area'] / (np.pi * obj.A * obj.e) * V_avg_sq * rho/2 * atmos.g / constants['weight']
    B = ap.mu_ground * constants['wing surface area'] * V_avg_sq * rho/2 * atmos.g / constants['weight']
    C = (constants['propeller power'] * constants['propeller efficiency'] / np.sqrt(V_avg_sq) - ap.mu_ground * 
         constants['weight'] *np.cos(np.radians(constants['runway slope'])) - obj.CD0 * rho/2 * V_avg_sq 
         * constants['wing surface area'] - constants['weight']*np.sin(np.radians(constants['runway slope']))) \
         * atmos.g / constants['weight'] - V_avg_sq/750
    
    sqrt = B**2 - 4*A*C
    C_L_TO_1 = (-B + sqrt) / (2*A)
    C_L_TO_2 = (-B - sqrt) / (2*A)

    return C_L_TO_1, C_L_TO_2


# ---------------- Run the plotting -----------------

# dictionary with constants:
hp_to_watt = 745.699872
# Plot for constant wind and different runway slopes, fixed runway slope with different wind speed with and against
dic_constants = {'runway slope': np.arange(0, 10),
    'airport altitude': 0, 'wing surface area': 11, 'weight': takeoffweight(aircraft, 200)*atm.g,
    'wind speed': 0, 'propeller power': aircraft.power*hp_to_watt, 'propeller efficiency': aircraft.eta_p}

plt_to = False
if plt_to:
    figure, axis = plt.subplots(2, 2)

    CL_TO_1 ,CL_TO_2 = TO_eom(aircraft, airfield, atm, dic_constants)
    axis[0, 0].plot(dic_constants['runway slope'], CL_TO_1)
    axis[0, 0].plot(dic_constants['runway slope'], CL_TO_2)
    axis[0, 0].set_title('runway slope vs C_L take-off')
    axis[0, 0].set_xlabel('runway slope[deg]')
    axis[0, 0].set_ylabel('C_L take-off [-]')

    dic_constants['runway slope'] = 0
    dic_constants['wind speed'] = np.arange(0, 10)
    CL_TO_1, CL_TO_2 = TO_eom(aircraft, airfield, atm, dic_constants)

    axis[1, 0].plot(dic_constants['wind speed'], CL_TO_1, color='red')
    axis[1, 0].plot(dic_constants['wind speed'], CL_TO_2, color='red')
    axis[1, 0].set_title('headwind vs C_L take-off')
    axis[1, 0].set_xlabel('headwind speed [m/sec]')
    axis[1, 0].set_ylabel('C_L take-off [-]')

    dic_constants['wind speed'] = np.arange(0, -10, -1)
    CL_TO_1, CL_TO_2 = TO_eom(aircraft, airfield, atm, dic_constants)

    axis[1, 1].plot(dic_constants['wind speed'], CL_TO_1, color='green')
    axis[1, 1].plot(dic_constants['wind speed'], CL_TO_2, color='green')
    axis[1, 1].set_title('tailwind vs C_L take-off')
    axis[1, 1].set_xlabel('tailwind speed [m/sec]')
    axis[1, 1].set_ylabel('C_L take-off [-]')

    dic_constants['wind speed'] = 0
    dic_constants['airport altitude'] = np.arange(0, 500)
    CL_TO_1, CL_TO_2 = TO_eom(aircraft, airfield, atm, dic_constants)

    axis[0, 1].plot(dic_constants['airport altitude'], CL_TO_1, color='black')
    axis[0, 1].plot(dic_constants['airport altitude'], CL_TO_2, color='black')
    axis[0, 1].set_title('airport altitude vs C_L take-off')
    axis[0, 1].set_xlabel('airport altitude [m]')
    axis[0, 1].set_ylabel('C_L take-off [-]')

    plt.subplots_adjust(hspace=0.6)
    plt.subplots_adjust(wspace=0.5)
    plt.suptitle('Take-off')
    plt.show()

# ------------------------------------------------------------------------

def turnperformance(ac_obj, atm_obj, phi, V, W, h, heading_change):
    rho = atm_parameters(atm_obj, h)[3] 
    phi *= np.pi/180
    R_turn = V**2/(atm_obj.g * np.tan(phi))
    n_turn = 1/np.cos(phi)
    CL = 2*W/(rho*V**2*ac_obj.Sw)
    CD = ac_obj.CD0 + CL**2/(np.pi*ac_obj.A*ac_obj.e)
    CL_CD = CL/CD
    Pr = n_turn*CL_CD*W*V
    turn_time = (2*np.pi*(heading_change/360)*R_turn)/V
    return R_turn, n_turn, Pr, turn_time

# -----------------------------------------------------------------

def cruiseperformance(ac_obj, atm_obj, Range=None, V_cruise=None, h_cruise=None, Payload_Range=False):
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
    W_cr    = ac_obj.W_TO * ac_obj.W1W_TO * ac_obj.W2W1 * ac_obj.W3W2 * ac_obj.W4W3 * atm_obj.g
    rho_cr  = atm_parameters(atm_obj, h_cruise)[2]
    r_it = 0.0
    t    = 0.0
    dt   = 0.1
    W    = W_cr
    while r_it < R:
        CL_cr   = 2*W/(rho_cr*V_cruise**2*ac_obj.Sw)
        CD_cr   = dragpolar(ac_obj, CL_cr)
        r_it    += (V_cruise * dt)
        t       += dt
        D       = 1/2 * rho_cr * V_cruise**2 * ac_obj.Sw * CD_cr
        P_req   = D*V_cruise/ac_obj.prop_eff
        F       = ac_obj.SFC * P_req
        W       -= (F*dt)
    print("---------------------------------------------------")
    print(f"Cruise performance - Range {R/1000} [km] - Cruise speed {V_cruise} [m/s] - Cruise height {h_cruise} [m]")
    print("---------------------------------------------------")
    print(f"The cruise time is {np.round(t, 2)} seconds ({np.round(t/3600, 2)} hours)")
    print(f"The fuel used during the cruise is {np.round(W_cr - W)} [kilograms] ({np.round(((W_cr-W)/0.7429), 2)} [L] @ {ac_obj.fueldensity} [kg/m^3])")
    print("---------------------------------------------------")
    return None


def payloadrange(ac_obj, atm_obj, V_cruise=None, h_cruise=None, plot=True):
    if V_cruise == None:
        V_cruise = ac_obj.V_cruise
    else:
        V_cruise = V_cruise
    if h_cruise == None:
        h_cruise = ac_obj.h_cruise
    else:
        h_cruise = h_cruise
    Fuel_loads  = np.arange(0, ac_obj.fuelcapacity, 1.0)
    Reserve     = ac_obj.M_res * ac_obj.fuelcapacity * ac_obj.fueldensity
    ZFW         = ac_obj.W_OE + Reserve + ac_obj.n_boxes * ac_obj.boxweight
    maxZFW_fuel = ac_obj.W_TO - ZFW
    ZFW_maxfuel = ac_obj.W_TO - ac_obj.fuelcapacity * ac_obj.fueldensity
    PL_maxfuel  = ZFW_maxfuel - Reserve - ac_obj.W_OE
    Ferryweight = ac_obj.W_OE + Reserve + ac_obj.fuelcapacity * ac_obj.fueldensity
    print("-------------------------------------------------------------------")
    print(f"Max ZFW: {np.round(ZFW, 2)} [kg]")
    print(f"Fuel @ max ZFW: {np.round(maxZFW_fuel, 2)} [kg] or {np.round(maxZFW_fuel / ac_obj.fueldensity, 2)} [L]")
    print(f"ZFW @ max fuel: {np.round(ZFW_maxfuel)} [kg]. The aircraft carries {np.round(PL_maxfuel)} [kg]")
    print(f"The TOW @ ferry configuration is {np.round(Ferryweight)} [kg]")
    print("-------------------------------------------------------------------")
    rho_cr = atm_parameters(atm_obj, h_cruise)
    W = ZFW
    for i in range(len(Fuel_loads)):
        CL_cr = 2 * W / (rho_cr * V_cruise ** 2 * ac_obj.Sw)
        CD_cr = dragpolar


# -------------------------------- LANDING -----------------------------------
def LA_eom(obj, ap, atmos, constants):

    p, T, rho, a = atm_parameters(obj, constants['airport altitude'])

    V_avg_sq = 0.72 * (np.sqrt(constants['weight'] / constants['wing surface area'] * 2 / rho * 1 / obj.CL_max_land) -
                          constants['wind speed']) ** 2

    A = - constants['wing surface area'] / (np.pi * obj.A * obj.e) * V_avg_sq * rho / 2 * atmos.g / constants['weight']
    B = ap.mu_ground * constants['wing surface area'] * V_avg_sq * rho / 2 * atmos.g / constants['weight']
    C = (800- ap.mu_ground * constants['weight'] * np.cos(np.radians(constants['runway slope'])) - obj.CD0 * rho / 2 * V_avg_sq
         * constants['wing surface area'] - constants['weight'] * np.sin(np.radians(constants['runway slope']))) \
        * atmos.g / constants['weight'] + V_avg_sq / 750

    sqrt = B ** 2 - 4 * A * C
    C_L_LA_1 = (-B + sqrt) / (2 * A)
    C_L_LA_2 = (-B - sqrt) / (2 * A)

    return C_L_LA_1, C_L_LA_2


# plot the results:
plt_land = False
if plt_land:

    figure, axis = plt.subplots(2, 2)

    dic_constants['weight'] = aircraft.W_OE * atm.g
    dic_constants['runway slope'] = np.arange(0, 10)
    dic_constants['airport altitude'] = 0
    CL_LA_1, CL_LA_2 = LA_eom(aircraft, airfield, atm, dic_constants)

    axis[0, 0].plot(dic_constants['runway slope'], CL_LA_1)
    axis[0, 0].plot(dic_constants['runway slope'], CL_LA_2)
    axis[0, 0].set_title('runway slope vs C_L take-off')
    axis[0, 0].set_xlabel('runway slope[deg]')
    axis[0, 0].set_ylabel('C_L take-off [-]')

    dic_constants['runway slope'] = 0
    dic_constants['wind speed'] = np.arange(0, 10)
    CL_LA_1, CL_LA_2 = LA_eom(aircraft, airfield, atm, dic_constants)

    axis[1, 0].plot(dic_constants['wind speed'], CL_LA_1, color='red')
    axis[1, 0].plot(dic_constants['wind speed'], CL_LA_2, color='red')
    axis[1, 0].set_title('headwind vs C_L take-off')
    axis[1, 0].set_xlabel('headwind speed [m/sec]')
    axis[1, 0].set_ylabel('C_L take-off [-]')

    dic_constants['wind speed'] = np.arange(0, -10, -1)
    CL_LA_1, CL_LA_2 = LA_eom(aircraft, airfield, atm, dic_constants)

    axis[1, 1].plot(dic_constants['wind speed'], CL_LA_1, color='green')
    axis[1, 1].plot(dic_constants['wind speed'], CL_LA_2, color='green')
    axis[1, 1].set_title('tailwind vs C_L take-off')
    axis[1, 1].set_xlabel('tailwind speed [m/sec]')
    axis[1, 1].set_ylabel('C_L take-off [-]')

    dic_constants['wind speed'] = 0
    dic_constants['airport altitude'] = np.arange(0,500)
    CL_LA_1, CL_LA_2 = LA_eom(aircraft, airfield, atm, dic_constants)

    axis[0, 1].plot(dic_constants['airport altitude'], CL_LA_1, color='black')
    axis[0, 1].plot(dic_constants['airport altitude'], CL_LA_2, color='black')
    axis[0, 1].set_title('airport altitude vs C_L take-off')
    axis[0, 1].set_xlabel('airport altitude [m]')
    axis[0, 1].set_ylabel('C_L take-off [-]')

    plt.subplots_adjust(hspace=0.6)
    plt.subplots_adjust(wspace=0.5)
    plt.suptitle('Landing')
    plt.show()

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


descend(aircraft, atm, 90, (aircraft.W_OE+100)*atm.g, 95, 500, 0.6)
