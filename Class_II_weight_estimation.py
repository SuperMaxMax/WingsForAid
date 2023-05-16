from math import cos

#########################################################################
"CLASS II WEIGHT ESTIMATION"
#########################################################################

def wing_weight(obj):
    k_w = 4.9E-3
    b_s = obj.b / cos(obj.lambda_mid)
    b_ref = 1.905
    t_r = obj.t_c * obj.rootchord

    W_w = k_w * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * obj.n_ult**0.55 * ((b_s/t_r)/(obj.W_TO/obj.Sw))**0.3 * obj.W_TO

    if obj.braced_wing == True:
        W_w *= 0.7
    if obj.pos_main_carriage == 'fuselage':
        W_w *= 0.95

    return W_w

def tail_weight(obj): #ultimate load factor, tail surface area
    W_t = 0.64 * (obj.n_ult * obj.s_tail**2)**0.75
    
    # IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT
    return W_t

def gear_weight(obj):
    if obj.main_gear_type=="fixed":
        A1 = 9.1
        B1 = 0.082
        C1 = 0.019
        D1 = 0
    else:
        A1 = 18.1
        B1 = 0.131
        C1 = 0.019
        D1 = 2.23 * 10**-5

    if obj.nose_gear_type=="retractable":
        A2 = 11.3
        B2 = 0
        C2 = 0.0024
        D2 = 0
    else:
        A2 = 9.1
        B2 = 0.082
        C2 = 0
        D2 = 2.97 * 10**-6

    # Main gear:
    W_uc1 = 1.08 * (A1 + B1 * obj.W_TO**0.75 + C1 * obj.W_TO + D1 * obj.W_TO**1.5)

    # Nose gear:
    W_uc2 = 1.08 * (A2 + B2 * obj.W_TO**0.75 + C2 * obj.W_TO + D2 * obj.W_TO**1.5)
    W_uc = W_uc1 + W_uc2

    return W_uc

def nacelle_weight(obj):
    W_n = 1.134 * obj.P_TO**0.5 # Take off power in hp

    return W_n

def equipment_weight(obj):
    W_eq = 0.08 * obj.W_TO

    #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF
    return W_eq

def fuselage_weight(obj):
    k_wf = 0.23
    W_fus = k_wf*(obj.V_D*(obj.l_t/(obj.h_out+obj.w_out)))**0.5*obj.S_G**1.2

    if obj.engine_pos == 'pusher':
        W_fus *= 1.04

    if obj.pos_main_carriage == 'fuselage':
        W_fus *= 1.07
    else:
        W_fus *= 0.96

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
    obj.W_w = wing_weight(obj)
    obj.W_t = tail_weight(obj)
    obj.W_uc = gear_weight(obj)
    obj.W_n = nacelle_weight(obj)
    obj.W_eq = equipment_weight(obj)
    obj.W_fus = fuselage_weight(obj)
    obj.W_sc = control_surface_weight(obj)
    obj.W_pg = propulsion_weight(obj)

    obj.W_OE = obj.W_w + obj.W_t + obj.W_uc + obj.W_n + obj.W_eq + obj.W_fus + obj.W_sc + obj.W_pg

    obj.W_TO = obj.W_OE + obj.W_F + obj.W_PL    