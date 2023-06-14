from math import sqrt, pi, exp
from parameters import UAV
import numpy as np

#########################################################################
"CLASS I WEIGHT ESTIMATION"                                             #
#########################################################################

def weight_take_off(W_OE, W_F, W_PL):
    W_TO = W_OE + W_F + W_PL

    return W_TO

def weight_empty_operational(obj):
    # Operational empty weigth is a function of the empty weight and the trapped fuel weight, which are both a function of take_off weight
    W_E = obj.lin_par1 * obj.W_TO + obj.lin_par2                        # Relation from Roskam or own relation
    W_tfo = 0.002 * obj.W_TO
    W_OE = W_E + W_tfo

    return W_OE

def L_D_cruise(obj):
    # Lift over drag calculation
    L_D = sqrt(pi * obj.A * obj.e / (4 * obj.CD0))
    obj.L_D = L_D

    return L_D

def W5W4_calculation(obj):
    # Cruise fuel fraction calculation
    L_D = L_D_cruise(obj) # Lift over drag ratio from method above

    W5W4 = 1/exp(obj.R * obj.g0 * obj.c_p / (obj.prop_eff * L_D)) # Breguet's range equation

    return W5W4

def W7W5_calculation(obj):
    # Dropping manoeuvre
    W7W5 = obj.W4W3 * obj.W10W9

    return W7W5

def Mff_calculation(obj):
    W7W5 = W7W5_calculation(obj)
    W5W4 = W5W4_calculation(obj)

    Mff = obj.W1W_TO * obj.W2W1 * obj.W3W2 * obj.W4W3 * W5W4 * W7W5**obj.n_drops * obj.W10W9 * obj.WfinalW10
    obj.Mff = Mff
    return Mff

def weight_fuel_used(obj):
    Mff = Mff_calculation(obj)

    M_f_used = (1 - Mff) * obj.W_TO

    return M_f_used

def weight_fuel(obj):
    M_f_used = weight_fuel_used(obj)

    W_F = (1 + obj.M_res) * M_f_used
    return W_F

def profile(obj, n_boxes, n_drops, Range = None, dropregion = None, result = False):
    W0 = obj.W_TO
    Mf_used = 0.0
    h_cruise= obj.h_cruise
    L_D = L_D_cruise(obj)
    if Range == None:
        Range = obj.R / 2                           
    else:
        Range = Range*1000
    if dropregion == None:
        W = W0
        # ======== Take-off Mff ========
        Mff_untilTO = obj.W1W_TO * obj.W2W1 * obj.W3W2
        Mf_used_TO  = (1 - Mff_untilTO) * W
        W   -= Mf_used_TO
        Mf_used += Mf_used_TO
        print(f"Mf_used after Take-off: {Mf_used}")
        # == Climb - Cruise - Descent ==                          In between drops. If n_drops = 1 this is just climb cruise and descent once
        interdropdist = Range / n_drops                         # Distance between drops in meters
        if interdropdist <= 10000:
            interdrop_h_cruise = 500                            # [m]
        elif 10000 <= interdropdist <= 50000:
            interdrop_h_cruise = 1000                           # [m]
        elif 50000 <= interdropdist <= 100000:
            interdrop_h_cruise = 1000 + 0.01 * interdropdist
        elif 100000 <= interdropdist <= 150000:     
            interdrop_h_cruise = 2000 + 0.01 * (interdropdist - 100000) # [m]
        else:
            interdrop_h_cruise = h_cruise                       # 3048 [m], 10000 [ft]
        for i in range(n_drops):
            W4W3 = obj.W4W3**(interdrop_h_cruise / 2286)        # Typical cruising altitude single engine props 5000 ft - 10000 ft, change climb fraction accordingly
            print(f"Mf_used after Climb: {(1-W4W3)*W}")
            W5W4 = 1/exp(interdropdist * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
            print(f"Mf_used after Cruise: {(1-W5W4)*W}")
            print(f"Mf_used after Descent: {(1-obj.W10W9)*W}")
            Mff_cl_cr_de = W4W3*W5W4*obj.W10W9
            Mf_used_cl_cr_de = (1 - Mff_cl_cr_de) * W
            W = W - Mf_used_cl_cr_de - (n_boxes/n_drops) * obj.boxweight
            Mf_used += Mf_used_cl_cr_de
        # ======= Return to base =======
        climbfrac   = obj.W4W3**(h_cruise / 2286)
        print(f"Mf_used after Climb there: {(1-W4W3)*W}")
        cruisefrac  = 1/exp(Range * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
        print(f"Mf_used after Cruise back: {(1-W5W4)*W}")
        Mff_return  = climbfrac * cruisefrac * obj.W10W9 * obj.WfinalW10
        print(f"Mf_used after Descent: {(1-obj.W10W9)*W}")
        Mf_used_return  = (1-Mff_return) * W
        W -= Mf_used_return
        Mf_used += Mf_used_return
    else:
        W = W0
        # ======== Take-off Mff ========
        Mff_untilTO = obj.W1W_TO * obj.W2W1 * obj.W3W2
        Mf_used_TO  = (1 - Mff_untilTO) * W
        W   -= Mf_used_TO
        Mf_used += Mf_used_TO
        # ==== Range to first drop ====='
        R_firstdrop = Range - dropregion*1000
        if R_firstdrop <= 10000:
            h_cruise_to1stdrop = 500                            # [m]
        elif 10000 <= R_firstdrop <= 50000:
            h_cruise_to1stdrop = 1000                           # [m]
        elif 50000 <= R_firstdrop <= 100000:
            h_cruise_to1stdrop = 1000 + 0.01 * R_firstdrop
        elif 100000 <= R_firstdrop <= 150000:     
            h_cruise_to1stdrop = 2000 + 0.01 * (R_firstdrop - 100000) # [m]
        else:
            h_cruise_to1stdrop = h_cruise
        interdropdist = (dropregion * 1000) / n_drops           # Distance between drops in meters
        if interdropdist <= 10000:
            interdrop_h_cruise = 500                            # [m]
        elif 10000 <= interdropdist <= 50000:
            interdrop_h_cruise = 1000                           # [m]
        elif 50000 <= interdropdist <= 100000:
            interdrop_h_cruise = 1000 + 0.01 * interdropdist
        elif 100000 <= interdropdist <= 150000:     
            interdrop_h_cruise = 2000 + 0.01 * (interdropdist - 100000) # [m]
        else:
            interdrop_h_cruise = h_cruise                       # 3048 [m], 10000 [ft]
        # == S/U - TO - Climb - Cruise - Descent ==               To first drop
        Mff_tofirstdrop = obj.W1W_TO * obj.W2W1 * obj.W3W2 * (obj.W4W3**(h_cruise_to1stdrop / 2286)) * 1/exp(R_firstdrop * obj.g0 * obj.c_p / (obj.prop_eff * L_D)) * obj.W10W9
        Mf_used_to1stdrop = (1 - Mff_tofirstdrop) * W
        Mf_used += Mf_used_to1stdrop
        W -= Mf_used_to1stdrop
        # ======== Dropping in dropregion =========
        for i in range(n_drops):
            W4W3 = obj.W4W3**(interdrop_h_cruise / 2286)        # Typical cruising altitude single engine props 5000 ft - 10000 ft, change climb fraction accordingly
            W5W4 = 1/exp(interdropdist * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
            Mff_cl_cr_de = W4W3*W5W4*obj.W10W9
            Mf_used_cl_cr_de = (1 - Mff_cl_cr_de) * W
            W = W - Mf_used_cl_cr_de - (n_boxes/n_drops) * obj.boxweight
            Mf_used += Mf_used_cl_cr_de
        # ======= Return to base =======
        climbfrac   = obj.W4W3**(h_cruise / 2286)
        cruisefrac  = 1/exp(Range * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
        Mff_return  = climbfrac * cruisefrac * obj.W10W9 * obj.WfinalW10
        Mf_used_return  = (1-Mff_return) * W
        W -= Mf_used_return
        Mf_used += Mf_used_return
    if result:
        if n_drops == 1:
            print(f"The fuel used to drop {n_boxes} boxes at a range of {np.round(Range/1000,2)} [km] is {np.round(Mf_used, 2)} [kg]")
        if dropregion == None and n_drops > 1:
            print(f"The fuel used to drop {n_boxes} boxes in {n_drops} evenly spaced drops (furthest {np.round(Range/1000,2)} [km]) is {np.round(Mf_used, 2)} [kg]")
        if dropregion != None:
            print(f"The fuel used to drop {n_boxes} boxes in {n_drops} within a dropregion of {dropregion} [km] (furthest {np.round(Range/1000,2)} [km]) is {np.round(Mf_used, 2)} [kg]")
    # ==== Calculate total Mff =====
    obj.Mff = 1 - (Mf_used / W0)
    return Mf_used      

# def run(obj, n_boxes, n_drops, Range = None, dropregion = None):
#     obj.W_F  = profile(obj, n_boxes, n_drops, Range, dropregion)    # perform calculation of fuel weight
#     obj.W_OE = weight_empty_operational(obj)                        # perform calculation of operational empty weight
#     obj.W_TO = weight_take_off(obj.W_OE, obj.W_F, obj.W_PL)         # combines weights to find total weight
 