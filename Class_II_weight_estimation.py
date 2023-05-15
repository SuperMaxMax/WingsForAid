import numpy as np
from parameters import *


#########################################################################
"CLASS II WEIGHT ESTIMATION"
#########################################################################

def wing_weight(W_TO):   #Span, sweep, ultimate load factor, thickness over chord, cord length at root, wing loading, gross weight
    k_w = 4.9E-3
    b_s = b / np.cos(lambda_mid)
    b_ref = 1.905
    t_r = t_c * cwr

    W_w = k_w * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * n_ult**0.55 * ((b_s/t_r)/W_loading)**0.3 * W_TO

    #ADD 30% IF BRACED WING USED, 10% IF STRUT USED?
    return W_w

def tail_weight(): #ultimate load factor, tail surface area
    W_t = 0.64 * (n_ult * s_tail**2)**0.75
    return W_t
    #IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT

def gear_weight(W_TO):

    maingear_type = input("Is the main gear fixed or retractable [f/r]",)
    nosegear_type = input("Is the nose gear fixed or retractable [f/r]",)

    if (maingear_type!="f" and maingear_type!="r") and (nosegear_type!="f" and nosegear_type!="r"):
        print("Input f or r for the gear type")
        return

    if maingear_type=="f":
        A1 = 9.1
        B1 = 0.082
        C1 = 0.019
        D1 = 0
    else:
        A1 = 18.1
        B1 = 0.131
        C1 = 0.019
        D1 = 2.23 * 10**-5

    if nosegear_type=="f":
        A2 = 11.3
        B2 = 0
        C2 = 0.0024
        D2 = 0
    else:
        A2 = 9.1
        B2 = 0.082
        C2 = 0
        D2 = 2.97 * 10**-6

    #Main gear:
    W_uc1 = 1.08 * (A1 + B1 * W_TO**0.75 + C1 * W_TO + D1 * W_TO**1.5)

    #Nose gear:
    W_uc2 = 1.08 * (A2 + B2 * W_TO**0.75 + C2 * W_TO + D2 * W_TO**1.5)
    W_uc = W_uc1 + W_uc2
    return W_uc

def nacelle_weight(): #Take off power in hp
    W_n = 1.134 * P_TO**0.5
    return W_n

def equipment_weight(W_TO):
    W_eq = 0.008 * W_TO
    return W_eq
    #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF

def fuselage_weight():
    k_wf = 0.23
    W_fus = k_wf*(V_D*(l_t/(b_f+h_f)))**0.5*S_G**1.2
    return W_fus
    # Add 7% if the main landing gear is attached to the fuselage, but 4% can be subtracted from the fuselage weight if there is no attachment structure for the main landing gear
    # Add 10% for freighter aircraft
    # For booms l_t is defined as distance between local wing chord and horizontal tailplane

def control_surface_weight(W_TO):
    k_sc = 0.44*0.768
    W_sc = k_sc*W_TO**(2/3)
    return W_sc
    # Add 20% for LE flap or slat
    # Add 15% for lift dumper controls

def propulsion_weight():
    k_pg = 1.16 # tractor single propeller aircraft
    W_pg = k_pg*N_e*(W_e+0.109*P_TO)
    return W_pg
    # If number of cylinder and volume of cylinder are known use figure 4-12 Torenbeek

def weight_empty(W_pg, W_sc, W_fus, W_eq, W_n, W_t, W_w, W_uc):
    W_OEW = W_pg + W_sc + W_fus + W_eq + W_n + W_t + W_w  + W_uc

    # Print all weights
    print(f"W_pg:{round(W_pg,3)} [kg]")
    print(f"W_sc:{round(W_sc,3)} [kg]")
    print(f"W_fus:{round(W_fus)} [kg]")
    print(f"W_eq:{round(W_eq)}")
    print(f"W_n:{round(W_n)}")
    print(f"W_t:{round(W_t)} [kg]")
    print(f"W_w:{round(W_w)} [kg]")
    print(f"W_uc:{round(W_uc)} [kg]")

    # Operative empty weight
    print(f"W_OEW:{W_OEW}")

    return W_OEW

def cg_calc():
    wing_cg = 0
    if sweep_angle == 0:
        wing_cg = 0.4 * cwr + l_LE #40% of root chord plus Leading Edge location
    else:
        wing_cg = 0 # to be done later, depends on spar locations (table 8-15 Torenbeek)

    fus_cg = 0.335 * l_f #32-35% of fuselage length
    tail_cg = 0.42 * cwr + l_LE #42% of root chord plus Leading Edge location

    engine_cg = 0 # to be done later
    # c.g. or rotax 912is is at 32.7 cm from the front of the engine

    landing_gear_cg = 0 # to be done later
    # can be at airplane c.g. -> iteration needed, or use location main and nose landing gear

    cg = (wing_cg * W_w + fus_cg * W_f + tail_cg * W_t + engine_cg * W_pg + landing_gear_cg * W_uc) / (W_w + W_f + W_t + W_pg + W_uc)

if __name__ == "__main__":
    W_TO = 600
    W_w = wing_weight(W_TO)
    W_t = tail_weight()
    W_uc = gear_weight(W_TO)
    W_n = nacelle_weight()
    W_eq = equipment_weight(W_TO)
    W_fus = fuselage_weight()
    W_sc = control_surface_weight(W_TO)
    W_pg = propulsion_weight()

    weight_empty(W_pg, W_sc, W_fus, W_eq, W_n, W_t, W_w, W_uc)