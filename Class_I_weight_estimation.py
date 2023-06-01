from math import sqrt, pi, exp
import numpy as np
from parameters import UAV

aircraft = UAV("Droppy", "tractor", boom=True, braced_wing=True)

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

def profile(obj, n_drops, n_boxes, furthest_drop = None, dropregion = None, Print=False):
    if furthest_drop == None:
        furthest_drop = obj.R / 2
    else:
        furthest_drop = furthest_drop*1000 #convert to meters
    if n_boxes%2 != 0:
        print("The amount of boxes you entered is not a multiple of two. Please restate the amount of boxes")
        n_boxes = int(input("Amount of boxes (multiple of 2): "))
    if int(n_boxes / n_drops) < 2:
        print("The boxes have to be dropped in (multiple) pairs, please restate the amount of drops")
        n_drops = int(input("Amount of drops: "))
    
    # n_drops is the amount of drops that are to be performed in the sortie
    # furthest drop is the dropzone located furthest away from the groundbase
    # dropregion is a circle in which the dropzones lie. This is centered around the furthest dropping point
    L_D = sqrt(pi * obj.A * obj.e / (4 * obj.CD0))
    W_TO= obj.W_TO
    if n_drops == 1:
        W4W3 = obj.W4W3**2
        R = 2 * furthest_drop                               # Fly to dropzone and back, hence 2 times furthest drop
        # On the way to the dropzone (payload as specified), Mff_there includes startup, taxi, take-off, climb, cruise, descent for drop
        Mff_cruise_there = 1/exp((R/2) * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
        Mff_there = obj.W1W_TO * obj.W2W1 * obj.W3W2 * W4W3 * Mff_cruise_there * obj.W10W9
        Mf_used_there = (1-Mff_there) * W_TO
        # On the way back to the airfield (payload is dropped), Mff_back includes climb from drop, cruise, descent, landing, taxi, shutdown
        W_back = W_TO - n_boxes * obj.boxweight - Mf_used_there
        Mff_cruise_back = 1/exp((R/2) * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
        Mff_back = W4W3 * Mff_cruise_back * obj.W10W9 * obj.WfinalW10
        Mf_used_back = (1-Mff_back) * W_back
        Mf_used = Mf_used_there + Mf_used_back
        obj.Mff = (W_TO - Mf_used)/W_TO
    else:
        if dropregion == None:
            Mf_used = 0.0
            interdropdist = furthest_drop/n_drops
            interdropcruiseheight = 11000 - 1000*n_drops
            climbfraction_power = ((interdropcruiseheight-5000)/5000) + 1
            W4W3 = obj.W4W3**climbfraction_power
            Mff_untilTO = obj.W1W_TO * obj.W2W1 * obj.W3W2
            Mf_used_TO  = (1-Mff_untilTO) * W_TO
            Mf_used += Mf_used_TO
            W = W_TO - Mf_used_TO  
            Mff_cruise_inter = 1/exp((interdropdist) * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
            Mff_drop = W4W3 * Mff_cruise_inter * obj.W10W9
            # Per drop: climb cruise descent, cruising at 5000 ft in between drops as climbing to 10k ft is not efficient for such a distance
            for i in range(n_drops):
                Mf_used_drop_i = (1-Mff_drop) * W
                Mf_used += Mf_used_drop_i
                W = W - Mf_used_drop_i - (n_boxes/n_drops)*obj.boxweight
            # Cruise back
            cruiseheight_back = 10000
            climbfraction_power = ((cruiseheight_back-5000)/5000) + 1
            Mff_cruise_back = 1/exp((furthest_drop) * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
            Mff_back = (W4W3**climbfraction_power) * Mff_cruise_back * obj.W10W9 * obj.WfinalW10
            Mf_used_back = (1-Mff_back)*W
            Mf_used += Mf_used_back
            obj.Mff = (W_TO - Mf_used)/W_TO
        else:
            dropregion *= 1000
            interdropdist   = dropregion/n_drops
            R_firstdrop     = furthest_drop - dropregion
            h_cruise_firstdrop = max((R_firstdrop/(250*1000))*10000, 1000)
            climbfraction_power = ((h_cruise_firstdrop-5000)/5000) + 1
            W4W3 = obj.W4W3**climbfraction_power
            W = W_TO
            Mf_used = 0.0
            Mff_untildrop1 = obj.W1W_TO*obj.W2W1*obj.W3W2*W4W3*1/exp((R_firstdrop) * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
            Mf_used += (1-Mff_untildrop1) * W_TO
            W -= Mf_used
            # Dropping in the overall dropregion
            interdropdist = dropregion/n_drops
            h_cruise_interdrop = max((interdropdist/(250*1000))*10000, 1000)
            climbfraction_power = ((h_cruise_interdrop -5000)/5000) + 1
            W4W3 = obj.W4W3**climbfraction_power
            Mff_cruise_inter = 1/exp((interdropdist) * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
            Mff_drop = W4W3 * Mff_cruise_inter * obj.W10W9
            for i in range(n_drops):
                Mf_used_drop_i = (1-Mff_drop) * W
                Mf_used += Mf_used_drop_i
                W = W - Mf_used_drop_i - (n_boxes/n_drops)*obj.boxweight
            cruiseheight_back = 10000
            climbfraction_power = ((cruiseheight_back-5000)/5000) + 1
            Mff_cruise_back = 1/exp((furthest_drop) * obj.g0 * obj.c_p / (obj.prop_eff * L_D))
            Mff_back = (W4W3**climbfraction_power) * Mff_cruise_back * obj.W10W9 * obj.WfinalW10
            Mf_used_back = (1-Mff_back)*W
            Mf_used += Mf_used_back
            obj.Mff = (W_TO - Mf_used)/W_TO
    if Print:
        print(f"Number of drops: {n_drops} | Number of boxes: {n_boxes} | Furthest drop: {furthest_drop/1000} [km] | Fuel used: {Mf_used} [kg] | Mff: {obj.Mff}")
    return Mf_used

def run(obj, n_drops, n_boxes, furthest_drop = None, dropregion = None):
    obj.W_OE = weight_empty_operational(obj)                    # perform calculation of operational empty weight
    obj.W_F  = profile(obj, n_drops, n_boxes, furthest_drop, dropregion)                          # perform calculation of fuel weight
    obj.W_TO = weight_take_off(obj.W_OE, obj.W_F, obj.W_PL)     # combines weights to find total weight
    print(f"Take-off Weight: {obj.W_TO} [kg] | OEW: {obj.W_OE} [kg] | Fuel used: {obj.W_F} [kg]")







