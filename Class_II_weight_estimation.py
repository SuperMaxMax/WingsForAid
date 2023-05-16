from math import cos
import matplotlib.pyplot as plt

#########################################################################
"CLASS II WEIGHT ESTIMATION"
#########################################################################

def wing_weight(obj):   #Span, sweep, ultimate load factor, thickness over chord, cord length at root, wing loading, gross weight
    k_w = 4.9E-3
    b_s = obj.b / cos(obj.lambda_mid)
    b_ref = 1.905
    t_r = obj.t_c * obj.cwr

    W_w = k_w * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * obj.n_ult**0.55 * ((b_s/t_r)/(obj.W_TO/obj.Sw))**0.3 * obj.W_TO

    #ADD 30% IF BRACED WING USED, 10% IF STRUT USED?
    return W_w

def tail_weight(obj): #ultimate load factor, tail surface area
    W_t = 0.64 * (obj.n_ult * obj.s_tail**2)**0.75
    
    #IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT
    return W_t

def gear_weight(obj):
    maingear_type = 'f'#input("Is the main gear fixed or retractable [f/r]",)
    nosegear_type = 'F'#input("Is the nose gear fixed or retractable [f/r]",)

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
    W_eq = 0.08 * obj.W_TO

    #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF
    return W_eq

def fuselage_weight(obj):
    k_wf = 0.23
    W_fus = k_wf*(obj.V_D*(obj.l_t/(obj.b_f+obj.h_f)))**0.5*obj.S_G**1.2

    if obj.engine_pos == 'pusher':
        W_fus *= 1.04

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

def cg_calc(obj):
    # --- Wing group
    # Wing
    if obj.sweep_angle == 0:
        wing_cg = 0.4 * obj.cwr                 # 40% of root chord plus Leading Edge location
    else:
        wing_cg = 0.4 * obj.cwr                 # to be done later, depends on spar locations (table 8-15 Torenbeek)

    # Control surfaces
    control_surfaces_cg = obj.x_lemac + obj.MAC_length  # guess for now
    
    W_wing_gr = obj.W_w + obj.W_sc
    x_wcg = (wing_cg*obj.W_w + control_surfaces_cg*obj.W_sc)/(W_wing_gr)

    # --- Fuselage group # propellor to be done
    # Fuselage and engine
    if obj.engine_pos == 'tractor':
        fus_cg = 0.45 * obj.l_f                 # educated guess
        engine_cg = 0.327                       # based on Rotax 912is (.g. or rotax 912is is at 327 mm, total length is 665.1 mm)
    elif obj.engine_pos == 'pusher':
        fus_cg = 0.55 * obj.l_f                 # educated guess
        engine_cg = obj.l_f - (0.6651-0.327)    # based on Rotax 912is
    elif obj.engine_pos == 'fuselage':
        fus_cg = 0.53 * obj.l_f                 # educated guess
        engine_cg = 0.8 * obj.l_f               # educated guess

    # Tail
    if obj.boom == True:
        tail_cg = obj.l_f + 0.9*obj.l_f_boom    # educated guess
    else:
        tail_cg = 0.9*obj.l_f                   # for now at 0.9 of fuselage length

    # Equipment
    eq_cg = 0.5*obj.l_f                        # educated guess

    # Nacelle
    # at engine cg

    # Undercarriage
    # For now: cg assumed to be at aircraft cg -> not taken into account for X_FCG, but is part of OEW

    W_fus_gr = obj.W_fus + obj.W_pg + obj.W_t + obj.W_eq + obj.W_n + obj.W_uc
    X_FCG = (fus_cg*obj.W_fus + engine_cg*obj.W_pg + tail_cg*obj.W_t + eq_cg*obj.W_eq + engine_cg*obj.W_n)/(W_fus_gr - obj.W_uc)

    # X_LEMAC and xc_OEW
    xc_OEW = obj.xc_OEW_p*obj.MAC_length
    X_LEMAC = X_FCG + obj.MAC_length * ((x_wcg/obj.MAC_length)*(W_wing_gr/W_fus_gr)-(xc_OEW)*(1+W_wing_gr/W_fus_gr))

    # Final CG
    W_OEW = W_wing_gr+W_fus_gr
    X_OEW = X_LEMAC + xc_OEW
    print(f"W_OEW = {W_OEW} N, X_OEW = {X_OEW} m")

    # Fuel
    W_fuel_wi = obj.W_F
    X_fuel_wi = X_LEMAC + 0.5*obj.MAC_length
    print(f"W_fuel_wi = {W_fuel_wi} N, X_fuel_wi = {X_fuel_wi} m")

    # Payload
    if obj.engine_pos == 'tractor':
        dist_front = 0.6651 + 0.6  # [m]
    elif obj.engine_pos == 'pusher':
        dist_front = 0.6
    elif obj.engine_pos == 'fuselage':
        dist_front = 0.6
    
    # 2 boxes in front
    W_2box_f = 1/6*obj.W_PL
    X_2box_f = dist_front + 0.60/2
    # 4 boxes in the front
    W_4box_f = 1/3*obj.W_PL
    X_4box_f = dist_front + 2*0.60/2
    # 2 boxes in back
    W_2box_b = 1/6*obj.W_PL
    X_2box_b = dist_front + 5*0.6 + 0.6/2
    # 4 boxes in back
    W_4box_b = 1/3*obj.W_PL
    X_4box_b = dist_front + 4*0.6 + 2*0.6/2
    # all boxes
    W_allbox = obj.W_PL
    X_allbox = dist_front + 6*0.6/2

    # Calculate points to plot
    # OEW + fuel
    W_OEW_fuel_frac = (W_OEW + W_fuel_wi)/obj.W_TO
    X_OEW_fuel = (W_OEW*X_OEW + W_fuel_wi*X_fuel_wi)/(W_OEW + W_fuel_wi)
    print(f"W_OEW_fuel_frac = {W_OEW_fuel_frac}, X_OEW_fuel = {X_OEW_fuel}")
    # OEW + fuel + 2 boxes in front
    W_OEW_fuel_2box_f_frac = W_OEW_fuel_frac + W_2box_f/obj.W_TO
    X_OEW_fuel_2box_f = (W_OEW_fuel_frac*X_OEW_fuel + W_2box_f*X_2box_f)/(W_OEW_fuel_frac + W_2box_f)
    print(f"W_OEW_fuel_2box_f_frac = {W_OEW_fuel_2box_f_frac}, X_OEW_fuel_2box_f = {X_OEW_fuel_2box_f}")
    # OEW + fuel + 4 boxes in front
    W_OEW_fuel_4box_f_frac = W_OEW_fuel_frac + W_4box_f/obj.W_TO
    X_OEW_fuel_4box_f = (W_OEW_fuel_frac*X_OEW_fuel + W_4box_f*X_4box_f)/(W_OEW_fuel_frac + W_4box_f)
    print(f"W_OEW_fuel_4box_f_frac = {W_OEW_fuel_4box_f_frac}, X_OEW_fuel_4box_f = {X_OEW_fuel_4box_f}")
    # OEW + fuel + 2 boxes in back
    W_OEW_fuel_2box_b_frac = W_OEW_fuel_frac + W_2box_b/obj.W_TO
    X_OEW_fuel_2box_b = (W_OEW_fuel_frac*X_OEW_fuel + W_2box_b*X_2box_b)/(W_OEW_fuel_frac + W_2box_b)
    print(f"W_OEW_fuel_2box_b_frac = {W_OEW_fuel_2box_b_frac}, X_OEW_fuel_2box_b = {X_OEW_fuel_2box_b}")
    # OEW + fuel + 4 boxes in back
    W_OEW_fuel_4box_b_frac = W_OEW_fuel_frac + W_4box_b/obj.W_TO
    X_OEW_fuel_4box_b = (W_OEW_fuel_frac*X_OEW_fuel + W_4box_b*X_4box_b)/(W_OEW_fuel_frac + W_4box_b)
    print(f"W_OEW_fuel_4box_b_frac = {W_OEW_fuel_4box_b_frac}, X_OEW_fuel_4box_b = {X_OEW_fuel_4box_b}")
    # OEW + fuel + all boxes
    W_OEW_fuel_allbox_frac = W_OEW_fuel_frac + W_allbox/obj.W_TO
    X_OEW_fuel_allbox = (W_OEW_fuel_frac*X_OEW_fuel + W_allbox*X_allbox)/(W_OEW_fuel_frac + W_allbox)
    print(f"W_OEW_fuel_allbox_frac = {W_OEW_fuel_allbox_frac}, X_OEW_fuel_allbox = {X_OEW_fuel_allbox}")

    # Plot each point
    xs = [X_OEW, X_fuel_wi, X_2box_f, X_4box_f, X_2box_b, X_4box_b, X_allbox]
    w_fracs = [W_OEW_fuel_frac, W_OEW_fuel_frac, W_OEW_fuel_2box_f_frac, W_OEW_fuel_4box_f_frac, W_OEW_fuel_2box_b_frac, W_OEW_fuel_4box_b_frac, W_OEW_fuel_allbox_frac]
    labels = ['OEW', 'OEW + Fuel', 'OEW + Fuel + 2 boxes front', 'OEW + Fuel + 4 boxes front', 'OEW + Fuel + 2 boxes back', 'OEW + Fuel + 4 boxes back', 'OEW + Fuel + all boxes']
    # plot points with labels
    for x, w, label in zip(xs, w_fracs, labels):
        plt.scatter(x, w, label=label)
    plt.xlabel('X_cg [m]')
    plt.ylabel('Mass fraction [-]')
    plt.grid()
    plt.legend()
    plt.show()

    # Save most forward and most aft and fully loaded c.g. in object
    obj.X_cg_fwd = min(xs)
    obj.X_cg_aft = max(xs)
    obj.X_cg_range = obj.X_cg_aft - obj.X_cg_fwd
    obj.X_cg_full = X_OEW_fuel_allbox
    