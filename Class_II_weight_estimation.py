from math import cos

#########################################################################
"CLASS II WEIGHT ESTIMATION"
#########################################################################

def wing_weight(obj):   #Span, sweep, ultimate load factor, thickness over chord, cord length at root, wing loading, gross weight
    k_w = 4.9E-3
    b_s = obj.b / cos(obj.lambda_mid)
    b_ref = 1.905
    t_r = obj.t_c * obj.cwr

    W_w = k_w * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * obj.n_ult**0.55 * ((b_s/t_r)/(obj.W_TO/obj.S))**0.3 * obj.W_TO

    #ADD 30% IF BRACED WING USED, 10% IF STRUT USED?
    return W_w

def tail_weight(obj): #ultimate load factor, tail surface area
    W_t = 0.64 * (obj.n_ult * obj.s_tail**2)**0.75
    
    #IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT
    return W_t

def gear_weight(obj):
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
    W_uc1 = 1.08 * (A1 + B1 * obj.W_TO**0.75 + C1 * obj.W_TO + D1 * obj.W_TO**1.5)

    #Nose gear:
    W_uc2 = 1.08 * (A2 + B2 * obj.W_TO**0.75 + C2 * obj.W_TO + D2 * obj.W_TO**1.5)
    W_uc = W_uc1 + W_uc2

    return W_uc

def nacelle_weight(obj): #Take off power in hp
    W_n = 1.134 * obj.P_TO**0.5

    return W_n

def equipment_weight(obj):
    W_eq = 0.008 * obj.W_TO

    #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF
    return W_eq

def fuselage_weight(obj):
    k_wf = 0.23
    W_fus = k_wf*(obj.V_D*(obj.l_t/(obj.b_f+obj.h_f)))**0.5*obj.S_G**1.2

    # Add 7% if the main landing gear is attached to the fuselage, but 4% can be subtracted from the fuselage weight if there is no attachment structure for the main landing gear
    # Add 10% for freighter aircraft
    # For booms l_t is defined as distance between local wing chord and horizontal tailplane
    return W_fus


def control_surface_weight(obj):
    k_sc = 0.44*0.768
    W_sc = k_sc*obj.W_TO**(2/3)

    # Add 20% for LE flap or slat
    # Add 15% for lift dumper controls
    return W_sc

def propulsion_weight(obj):
    k_pg = 1.16 # tractor single propeller aircraft
    W_pg = k_pg*obj.N_e*(obj.W_e+0.109*obj.P_TO)

    # If number of cylinder and volume of cylinder are known use figure 4-12 Torenbeek
    return W_pg

def weight_empty(obj):
    W_w = wing_weight(obj)
    W_t = tail_weight(obj)
    W_uc = gear_weight(obj)
    W_n = nacelle_weight(obj)
    W_eq = equipment_weight(obj)
    W_fus = fuselage_weight(obj)
    W_sc = control_surface_weight(obj)
    W_pg = propulsion_weight(obj)

    obj.W_OE = W_pg + W_sc + W_fus + W_eq + W_n + W_t + W_w  + W_uc

    # # Print all weights
    # print(f"W_pg:{round(W_pg,3)} [kg]")
    # print(f"W_sc:{round(W_sc,3)} [kg]")
    # print(f"W_fus:{round(W_fus)} [kg]")
    # print(f"W_eq:{round(W_eq)}")
    # print(f"W_n:{round(W_n)}")
    # print(f"W_t:{round(W_t)} [kg]")
    # print(f"W_w:{round(W_w)} [kg]")
    # print(f"W_uc:{round(W_uc)} [kg]")

    # # Operative empty weight
    # print(f"W_OEW:{obj.W_OE}")

    # return obj.W_OE

# def cg_calc():
#     wing_cg = 0
#     if sweep_angle == 0:
#         wing_cg = 0.4 * cwr + l_LE #40% of root chord plus Leading Edge location
#     else:
#         wing_cg = 0 # to be done later, depends on spar locations (table 8-15 Torenbeek)

#     fus_cg = 0.335 * l_f #32-35% of fuselage length
#     tail_cg = 0.42 * cwr + l_LE #42% of root chord plus Leading Edge location

#     engine_cg = 0 # to be done later
#     # c.g. or rotax 912is is at 32.7 cm from the front of the engine

#     landing_gear_cg = 0 # to be done later
#     # can be at airplane c.g. -> iteration needed, or use location main and nose landing gear

#     cg = (wing_cg * W_w + fus_cg * W_f + tail_cg * W_t + engine_cg * W_pg + landing_gear_cg * W_uc) / (W_w + W_f + W_t + W_pg + W_uc)