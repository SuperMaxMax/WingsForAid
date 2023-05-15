from math import sqrt, pi, exp

#########################################################################
"CLASS I WEIGHT ESTIMATION"
#########################################################################


def weight_take_off(W_OE, W_F, W_PL):
    # take_off weight is a function of operational empty weight, fuel weight and payload weight
    W_TO = W_OE + W_F + W_PL

    return W_TO

def weight_empty_operational(obj):
    # Operational empty weigth is a function of the empty weight and the trapped fuel weight, which are both a function of take_off weight
    W_E = obj.lin_par1 * obj.W_TO + obj.lin_par2  # relation from Roskam
    W_tfo = 0.002 * obj.W_TO
    W_OE = W_E + W_tfo

    return W_OE

def L_D_cruise(obj):
    # Lift over drag calculation
    L_D = sqrt(pi * obj.A * obj.e / (4 * obj.CD0))

    return L_D

def W5W4_calculation(obj):
    # Cruise fuel fraction calculation
    L_D = L_D_cruise(obj) # Lift over drag ratio from method above
    W5W4 = 1/exp(obj.R * obj.g * obj.c_p / (obj.prop_eff * L_D)) #Breguet's range equation

    return W5W4

def W7W5_calculation(obj):
    # Dropping manoeuvre
    W7W5 = obj.W4W3 * obj.W10W9

    return W7W5

def Mff_calculation(obj):
    W7W5 = W7W5_calculation(obj)
    W5W4 = W5W4_calculation(obj)

    Mff = obj.W1W_TO * obj.W2W1 * obj.W3W2 * obj.W4W3 * W5W4 * W7W5**obj.n_drops * obj.W10W9 * obj.WfinalW10

    return Mff

def weight_fuel_used(obj):
    Mff = Mff_calculation(obj)

    M_f_used = (1 - Mff) * obj.W_TO

    return M_f_used

def weight_fuel(obj):
    M_f_used = weight_fuel_used(obj)

    W_F = (1 + obj.M_res) * M_f_used

    return W_F

def iteration(obj):
    it = True
    while it:
        W_F = weight_fuel(obj)                        # perform calculation of fuel weight
        W_OE = weight_empty_operational(obj)           # perform calculation of operational empty weight
        W_TO_new = weight_take_off(W_OE, W_F, obj.W_PL)         # combines weights to find total weight
        change = (W_TO_new - obj.W_TO)/obj.W_TO
        if abs(change) < 0.001:
            obj.W_TO = W_TO_new    # change between iteration is smaller than 0.1 percent

            it = False
        else:
            obj.W_TO = W_TO_new
            obj.W_OE = W_OE
            obj.W_F = W_F

    print(f"W_OE:{obj.W_OE}")
    print(f"W_F:{obj.W_F}")
    print(f"W_PL:{obj.W_PL}")
    print(f"W_TO:{obj.W_TO}")