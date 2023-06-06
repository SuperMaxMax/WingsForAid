import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects
from scipy.integrate import quad
import math


# Setting with which can be determined if first the aircraft is filled with fuel or with payload
fuel_first = False

def cg_calc(obj):
    # --- Wing group
    # Wing
    if obj.lambda_co4 == 0:
        wing_cg = 0.4 * obj.rootchord                 # 40% of root chord plus Leading Edge location
    else:
        wing_cg = 0.4 * obj.rootchord                 # to be done later, depends on spar locations (table 8-15 Torenbeek)

    # Control surfaces
    control_surfaces_cg = obj.x_lemac + 0.9*obj.MAC_length  # guess for now
    
    W_wing_gr = obj.W_w + obj.W_sc
    x_wcg = (wing_cg*obj.W_w + control_surfaces_cg*obj.W_sc)/(W_wing_gr)

    

    # --- Fuselage group # propellor to be done
    # Fuselage and engine
    prop_correction = 0.06                      # correction for propeller weight
    if obj.engine_pos == 'tractor':
        fus_cg = 0.48 * obj.l_f                 # educated guess
        engine_cg = obj.engine_cg - prop_correction     # based on Rotax 912is (.g. or rotax 912is is at 327 mm, total length is 665.1 mm)
    elif obj.engine_pos == 'pusher':
        fus_cg = 0.52 * obj.l_f                 # educated guess
        engine_cg = obj.l_f - (obj.engine_length-obj.engine_cg) + prop_correction # based on Rotax 912is
    elif obj.engine_pos == 'fuselage':
        fus_cg = 0.5 * obj.l_f                  # educated guess
        engine_cg = 0.8 * obj.l_f               # educated guess

    # Tail and boom
    if obj.boom == True:
        tail_cg = obj.l_f + 0.9*obj.l_f_boom    # educated guess
        boom_cg = obj.l_f + 0.5*obj.l_f_boom    # educated guess
    else:
        tail_cg = 0.9*obj.l_f                   # for now at 0.9 of fuselage length

    # Equipment
    eq_cg = 0.5*obj.l_f                         # educated guess

    # Nacelle
    nacelle_cg = engine_cg                      # nacelle cg assumed to be at engine cg

    # Undercarriage
    # For now: cg assumed to be at aircraft cg -> not taken into account for X_FCG, but is part of OEW

    if obj.boom == True:
        W_fus_gr = obj.W_fus + obj.W_pg + obj.W_t + obj.W_eq + obj.W_n + obj.W_uc + obj.W_boom
        X_FCG = (fus_cg*obj.W_fus + engine_cg*obj.W_pg + tail_cg*obj.W_t + eq_cg*obj.W_eq + nacelle_cg*obj.W_n + boom_cg*obj.W_boom)/(W_fus_gr - obj.W_uc)
    else:
        W_fus_gr = obj.W_fus + obj.W_pg + obj.W_t + obj.W_eq + obj.W_n + obj.W_uc
        X_FCG = (fus_cg*obj.W_fus + engine_cg*obj.W_pg + tail_cg*obj.W_t + eq_cg*obj.W_eq + nacelle_cg*obj.W_n)/(W_fus_gr - obj.W_uc)

    # X_LEMAC and xc_OEW
    xc_OEW = obj.xc_OEW_p*obj.MAC_length
    #X_LEMAC = X_FCG + obj.MAC_length * ((x_wcg/obj.MAC_length)*(W_wing_gr/W_fus_gr)-(xc_OEW)*(1+W_wing_gr/W_fus_gr))
    # X_LEMAC = obj.X_LEMAC
    X_LEMAC = 0.35 * (obj.l_f + obj.l_f_boom)
    obj.X_LEMAC = X_LEMAC
    aircraft.X_FCG = X_FCG
    print(aircraft.X_FCG)
    obj.X_FCG = X_FCG

    # Final CG
    W_OEW = W_wing_gr+W_fus_gr
    X_OEW = X_LEMAC + xc_OEW

    # Fuel
    W_fuel_wi = obj.W_F
    X_fuel_wi = X_LEMAC + 0.5*obj.MAC_length

    # Payload
    if obj.engine_pos == 'tractor':
        dist_front = obj.engine_length + 0.6  # [m]
    elif obj.engine_pos == 'pusher':
        dist_front = 0.4
    elif obj.engine_pos == 'fuselage':
        dist_front = 0.4

    if obj.n_boxes_abreast == 2:
        box_configs = [[2,0,0,0,0,0], [0,2,0,0,0,0], [0,0,2,0,0,0], [0,0,0,2,0,0], [0,0,0,0,2,0], [0,0,0,0,0,2], [2,2,0,0,0,0], [0,2,2,0,0,0], [0,0,2,2,0,0], [0,0,0,2,2,0], [0,0,0,0,2,2], [2,2,2,0,0,0], [0,2,2,2,0,0], [0,0,2,2,2,0], [0,0,0,2,2,2], [2,2,2,2,0,0], [0,2,2,2,2,0], [0,0,2,2,2,2], [2,2,2,2,2,0], [0,2,2,2,2,2], [2,2,2,2,2,2]]
        labels = ['200000', '020000', '002000', '000200', '000020', '000002', '220000', '022000', '002200', '000220', '000022', '222000', '022200', '002220', '000222', '222200', '022220', '002222', '222220', '022222', '222222']
        box_xs = [dist_front+0.2, dist_front+0.8, dist_front+1.4, dist_front+2.0, dist_front+2.6, dist_front+3.2]   
    elif obj.n_boxes_abreast == 3:
        box_configs = [[3,0,0,0], [0,3,0,0], [0,0,3,0], [0,0,0,3], [3,3,0,0], [0,3,3,0], [0,0,3,3], [3,3,3,0], [0,3,3,3], [3,3,3,3]]
        labels = ['3000', '0300', '0030', '0003', '3300', '0330', '0033', '3330', '0333', '3333']
        box_xs = [dist_front+0.2, dist_front+0.8, dist_front+1.4, dist_front+2.0]

    box_weights = [sum(i)*20 for i in box_configs]
    box_xcg_positions = [np.dot(i, box_xs)/sum(i) for i in box_configs]

    # OEW + fuel
    W_OEW_fuel_frac = (W_OEW + W_fuel_wi)/obj.W_TO
    X_OEW_fuel = (W_OEW*X_OEW + W_fuel_wi*X_fuel_wi)/(W_OEW + W_fuel_wi)

    # OEW + fuel + box configurations
    W_OEW_fuel_box_frac = [W_OEW_fuel_frac + i/obj.W_TO for i in box_weights]
    X_OEW_fuel_box = [((W_OEW+W_fuel_wi)*X_OEW_fuel + i*X_box)/(W_OEW+W_fuel_wi+i) for i, X_box in zip(box_weights, box_xcg_positions)]

    # OEW + box configurations
    W_OEW_box_frac = [W_OEW/obj.W_TO + i/obj.W_TO for i in box_weights]

    X_OEW_box = [(W_OEW*X_OEW + (i)*(X_box))/(W_OEW+i) for i, X_box in zip(box_weights, box_xcg_positions)]


    # Plot each point
    if fuel_first:
        Xs = [X_OEW, X_OEW_fuel] + X_OEW_fuel_box
        Xs = (np.array(Xs)-X_LEMAC)/obj.MAC_length
        w_fracs = [W_OEW/obj.W_TO, W_OEW_fuel_frac] + W_OEW_fuel_box_frac
        labels = ['OEW', 'OEW + Fuel'] + labels

    else:
        Xs = [X_OEW] + X_OEW_box + [X_OEW_fuel_box[-1]]
        Xs = np.array(Xs)
        Xs = (Xs-X_LEMAC)/obj.MAC_length
        w_fracs = [W_OEW/obj.W_TO] + W_OEW_box_frac + [W_OEW_fuel_box_frac[-1]]
        labels = ['OEW'] + labels + ['OEW + Box + Fuel']
        
    # plt.rcParams.update({'font.size': 14})
    for x, w, label, i in zip(Xs, w_fracs, labels, range(len(Xs))):
        plt.scatter(x, w)
        if i == 0 or i == len(Xs):
            rotation_t = 0
        else:
            rotation_t = 90
        plt.annotate(label, (x, w), textcoords="offset points", xytext=(0,10), ha='center', rotation=rotation_t)
            
    # Save most forward and most aft and fully loaded c.g. in object
    obj.X_cg_full = Xs[-1]
    obj.X_cg_range = max(Xs) - 0.20
    obj.X_cg_fwd = 0.20 - obj.X_cg_range * 0.05
    obj.X_cg_aft = max(Xs) + obj.X_cg_range * 0.05

    obj.AE_l_h = obj.l_f - (obj.X_LEMAC+ obj.X_cg_aft*obj.MAC_length) + obj.l_f_boom - 3/4 * obj.AE_rootchord_h
    print(obj.AE_l_h)
    print(obj.l_f + obj.l_f_boom)
    # Plot lines for forward and aft cg positions
    # plt.axvline(x=0, linestyle='--', color='red', label='0% MAC')
    # plt.axvline(x=1, linestyle=':', color='red', label='100% MAC')

    plt.axvline(x=obj.X_cg_fwd, linestyle='--', color='blue', label='most forward c.g. considered')
    plt.axvline(x=obj.X_cg_aft, linestyle='--', color='red', label='most aft c.g. considered')
    plt.xlim((0, 1))
    plt.xlabel('X_cg/MAC [-]', fontsize=12)
    plt.ylabel('Mass fraction [-]', fontsize=12)
    plt.grid()
    plt.legend()
    plt.title(f'Mass fraction vs X_cg/MAC for {obj.name}', loc='left')
    plt.show()

    return max(Xs), min(Xs), obj.X_cg_range
