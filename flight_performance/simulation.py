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

def propthrust(obj, h, ):
    p, T, rho, a = atm_parameters(atm, h)

def climbrate(ac_obj, atm_obj, W_F, V, P_climb, plot=True):
    start_time = time.time()
    V *= 0.5144
    atm_parameters_vectorized = np.vectorize(lambda h: atm_parameters(atm_obj, h))
    alt_range = np.arange(0, ac_obj.ceiling, 0.5)
    atm_obj.p, atm_obj.T, atm_obj.rho, atm_obj.a = atm_parameters_vectorized(alt_range)
    W       = takeoffweight(ac_obj, W_F) * atm_obj.g
    ROC     = np.empty(0)
    Gamma   = np.empty(0)
    for i in range(len(alt_range)):
        p, T, rho, a = atm_obj.p[i], atm_obj.T[i], atm_obj.rho[i], atm_obj.a[i]
        # L = W gives the following
        CL = 2*W/(rho*ac_obj.Sw*V**2)
        CD = dragpolar(ac_obj, CL)
        D  = 1/2 * rho * V**2 * ac_obj.Sw * CD
        Pr = D*V
        Pa = ac_obj.power * P_climb * ac_obj.prop_eff * 745.699872 * (rho/atm_obj.rho0)**(3/4)      # Convert to Watts
        roc= (Pa - Pr)/W                                                                            # Climb angle in degrees
        gamma = (roc/V)*(180/np.pi)
        FMF= ac_obj.SFC * ac_obj.power * P_climb * 745.699872
        if i != 0:
            dt = (alt_range[i]-alt_range[i-1]) / (roc)
        else:
            dt = (alt_range[i+1]-alt_range[i]) / (roc)
        dWF= FMF * dt
        W -= dWF
        ROC= np.append(ROC, roc)
        Gamma = np.append(Gamma, gamma)
    end_time = time.time()
    if plot:
        plt.plot(alt_range, ROC, label=f"Climb power: {P_climb*100}%")
        plt.xlabel("Altitude [m]")
        plt.ylabel("Rate of Climb [m/s]")
        plt.legend()
        plt.show()
        plt.plot(alt_range, Gamma, label=f"Climb power: {P_climb*100}%")
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
    alt_range = (0, ac_obj.th_ceil)
    atm_parameters_vectorized = np.vectorize(lambda h: atm_parameters(atm_obj, h))
    alt_range = np.arange(0, ac_obj.ceiling, 0.5)
    atm_obj.p, atm_obj.T, atm_obj.rho, atm_obj.a = atm_parameters_vectorized(alt_range)
    stall_limit = np.empty(0)
    thr_lim_lo = np.empty(0)
    thr_lim_hi  = np.empty(0)
    for i in range(len(alt_range)):
        p, T, rho, a = atm_obj.p[i], atm_obj.T[i], atm_obj.rho[i], atm_obj.a[i]
        V_s = np.sqrt(2*W/(rho*ac_obj.Sw*ac_obj.CL_max_clean))
        P_a = ac_obj.power * ac_obj.prop_eff * 745.699872 * (rho/atm_obj.rho0)**(3/4)
        V_max = np.arange(150, 250, 0.5)*0.5144
        V_min = np.arange(50, 150, 0.5)*0.5144
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
        plt.plot(stall_limit[stall_limit >= thr_lim_lo], alt_range[stall_limit >= thr_lim_lo], color='green')
        plt.plot(thr_lim_lo[thr_lim_lo >= V_s], alt_range[thr_lim_lo >= V_s], color = 'green')
        plt.xlabel("Airspeed [m/s]")
        plt.ylabel("Altitude [m]")
        plt.show()
    return

# a = flightceiling(aircraft, atm, 60)



# ---------------- Assumptions for take-off equations of motion -----------------
# Wind is included by take it into account in the speed: V_eff = V - V_wind
# Runway slope is not zero
# The present of rain is taken into account in the friction coefficient with the ground, mu
# delta_rw = runway slope
# D_g = force due to the ground friction, with
# Thrust and lift are taken as average values

def TO_eom(obj, ap, atmos, constants):

    print(type(constants['runway slope']))
    p, T, rho, a = atm_parameters(obj, constants['runway altitude'])
    V_min = np.sqrt((constants['weight']*np.cos(np.radians(constants['runway slope'])))/constants['wing surface area'] * 2/rho * 1/obj.CL_max_TO) - constants['wind speed']
    V_LOF = 1.05 * V_min
    V_avg = V_LOF / np.sqrt(2)

    # Perpendicular to the runway:
    L_avg = obj.CL_TO * 0.5 * rho * ((V_avg)**2) * constants['wing surface area']
    N = constants['weight'] * np.cos(np.radians(constants['runway slope'])) - L_avg

    # Parallel to the runway:
    D_g = ap.mu_ground * N
    C_D = obj.CD0 + obj.CL_TO**2 / (np.pi * obj.A * obj.e)
    D = C_D * 0.5 * rho * V_avg**2 * constants['wing surface area']
    T_avg = constants['propeller power'] * constants['propeller efficiency'] / (V_avg)  # Form "Aircraft performance and design" page 457
    acc = atmos.g / constants['weight'] * (T_avg - D - D_g - constants['weight']*np.sin(np.radians(constants['runway slope'])))

    # lift off distance:
    s_LO = V_LOF**2 / (2 * acc)
    print(s_LO)
    # plot lift off distance to runway slope:
    return s_LO


# ---------------- Run the plotting -----------------
# dictionary with constants:
hp_to_watt = 745.699872
# Plot for constant wind and different runway slopes, fixed runway slope with different wind speed with and against
dic_constants = {'runway slope': 0, 'runway altitude': 0, 'wing surface area': 11, 'weight':
    takeoffweight(aircraft, 200)*atm.g, 'wind speed': np.arange(0,10), 'propeller power': aircraft.power*hp_to_watt # 80*hp_to_watt, much lower than 115000!,
    ,'propeller efficiency': aircraft.eta_p}

figure, axis = plt.subplots(2, 2)

s_LO = TO_eom(aircraft, airfield, atm, {'runway slope': np.array([0,1, 2,3,4,5,6,7,8,9,10]),
    'runway altitude': 0, 'wing surface area': 11, 'weight': takeoffweight(aircraft, 200)*atm.g,
    'wind speed': 0, 'propeller power': aircraft.power*hp_to_watt, 'propeller efficiency': aircraft.eta_p})
axis[0, 0].plot(np.array([0,1, 2,3,4,5,6,7,8,9,10]), s_LO)
axis[0, 0].set_title('runway slope vs runway length')
axis[0, 0].set_xlabel('runway slope[deg]')
axis[0, 0].set_ylabel('runway length [m]')

axis[1, 0].plot(np.arange(0, 10), TO_eom(aircraft, airfield, atm, {'runway slope': 0, 'runway altitude': 0,
    'wing surface area': 11, 'weight': takeoffweight(aircraft, 200)*atm.g, 'wind speed': np.arange(0, 10),
    'propeller power': aircraft.power*hp_to_watt, 'propeller efficiency': aircraft.eta_p}), color='red')
axis[1, 0].set_title('headwind vs runway length')
axis[1, 0].set_xlabel('headwind speed [m/sec]')
axis[1, 0].set_ylabel('runway length [m]')

axis[1, 1].plot(np.arange(-10, 0), TO_eom(aircraft, airfield, atm, {'runway slope': 0, 'runway altitude': 0,
    'wing surface area': 11, 'weight': takeoffweight(aircraft, 200)*atm.g, 'wind speed': np.arange(-10, 0),
    'propeller power': aircraft.power*hp_to_watt, 'propeller efficiency': aircraft.eta_p}), color='green')
axis[1, 1].set_title('tailwind vs runway length')
axis[1, 1].set_xlabel('tailwind speed [m/sec]')
axis[1, 1].set_ylabel('runway length [m]')

axis[0, 1].plot(np.arange(0,500), TO_eom(aircraft, airfield, atm, {'runway slope': 0,
    'runway altitude': np.arange(0,500), 'wing surface area': 11, 'weight': takeoffweight(aircraft, 200)*atm.g,
    'wind speed': 0, 'propeller power': aircraft.power*hp_to_watt, 'propeller efficiency': aircraft.eta_p}), color='black')
axis[0, 1].set_title('airport altitude vs runway length')
axis[0, 1].set_xlabel('airport altitude [m]')
axis[0, 1].set_ylabel('runway length [m]')

plt.show()

