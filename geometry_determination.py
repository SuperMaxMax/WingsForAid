import  numpy as np
from    parameters import *
import  matplotlib.pyplot as plt

def metertofeet(meter):
    feet = meter*3.2808399
    return feet


# W/S and W/P diagrams
# drag polar
def dragpolar(CL, CD0, e, A):
    CD = CD0 + (CL**2/(np.pi*A*e))
    return CD

def stallWS(V, rho, CL_max):
    WoS = 1/2 * rho * V**2 * CL_max
    return WoS

def TOP_calc(WS, sigma, CL_TO, BHP, W_TO):
    TOP = WS/(sigma*CL_TO*(BHP/W_TO))
    return TOP

def altitude_effects(h, Lambda, R, g, T0, rho0, BHP0):
    rho     = rho0*(1+ ((Lambda*h)/T0))**(-((g/(R*Lambda))+1))
    sigma   = rho/rho0
    BHP     = BHP0*(sigma)**(3/4)
    return rho, sigma, BHP

def WP_WSdiagrams(h, plot=True):
    # find atmospheric properties and power at altitude according to ISA
    rho, sigma, BHP = altitude_effects(h, Lambda, R, g0, T0, rho0, P_max)
    # define a W/S array to make graphs later on
    WS              = np.arange(100.0, 2201.0, 1.0)
    # vectorize function to apply to every CL_max_clean considered
    WS_v            = np.vectorize(stallWS)
    WS_stall_max    = WS_v(V_s_min, rho, CL_max_clean)
    # Take off parameter
    TOP_req = 250   #read from graph of ADSEE 1 slides, slide 29 lecture 3
    TW = np.empty(0)
    for i in range(len(CL_TO)):
        TW_TO   = (250/WS)*CL_TO[i]*sigma
        if i == 0:
            TW = np.append(TW, TW_TO)
        else:
            TW = np.vstack((TW, TW_TO))
    WS_ldg = (CL_LDG*rho*(LDG_dist/0.5915))/(2*f)
    #W/P cruise
    rho_cruise, sigma_cruise, BHP_cruise = altitude_effects(h_cruise, Lambda, R, g0, T0, rho0, P_max)
    WP_TO   = (power_setting/cruise_frac)*eta_p*sigma**(3/4)*(((CD0*1/2*rho_cruise*V_cruise**3)/(WS*cruise_frac))+(WS*(1/(np.pi*A*e*rho_cruise*V_cruise))))**(-1)
    #Climb performance
    climb_rate = (8.3/100)*V_climb
    WP_Climb = eta_p/(climb_rate+((np.sqrt(WS)*np.sqrt(2/rho_cruise))/(1.345*(A*e)**(3/4)/CD0**(1/4))))
    #plotting
    if plot:
        for i in range(len(TW)):
            plt.plot(WS, TW[i])
        WP_plotting = np.arange(0, 0.51, 0.01)
        for j in range(len(WS_ldg)):
            WS_ldg_plot = np.full(len(WP_plotting), WS_ldg[j])
            plt.plot(WS_ldg_plot, WP_plotting)
        for k in range(len(WS_stall_max)):
            WS_stall_plot = np.full(len(WP_plotting), WS_stall_max[k])
        plt.plot(WS, WP_TO)
        plt.plot(WS, WP_Climb)
        plt.xlabel("W/S [N/m^2]")
        plt.ylabel("W/P [N/W]")
        plt.ylim((0, 1.5))
        plt.show()
a = WP_WSdiagrams(0, plot=True)

    









