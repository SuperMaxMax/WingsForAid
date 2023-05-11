import numpy as np



def wing_weight(b, lambda_mid, n_ult, t_c, cwr, W_loading, W_G):   #Span, sweep, ultimate load factor, thickness over chord, cord length at root, wing loading, gross weight
    k_w = 4.9*10**(-3)
    b_s = b / np.cos(lambda_mid)
    b_ref = 1.905
    t_r = t_c * cwr

    W_w = kw * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * n_ult**0.55 * ((b_s/t_r)/W_loading)**0.3 * W_G
    
    #ADD 30% IF BRACED WING USED, 10% IF STRUT USED?
    
    return W_w
    
def tail_weight(n_ult, s_tail): #ultimate load factor, tail surface area
    W_t = 0.64 * (n_ult * s_tail**2)**0.75

    #IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT
    return W_t

def gear_weight(W_TO):
    A1 = input("Give A for main gear (Table 8-6 Torenbeek)", )
    B1 = input("Give B for main gear (Table 8-6 Torenbeek)", )
    C1 = input("Give C for main gear (Table 8-6 Torenbeek)", )
    D1 = input("Give D for main gear (Table 8-6 Torenbeek)", )
    A2 = input("Give A for nose gear (Table 8-6 Torenbeek)", )
    B2 = input("Give B for nose gear (Table 8-6 Torenbeek)", )
    C2 = input("Give C for nose gear (Table 8-6 Torenbeek)", )
    D2 = input("Give D for nose gear (Table 8-6 Torenbeek)", )

    #Main gear:
    W_uc1 = 1.08 * (A + B * W_TO**0.75 + C * W_TO + D * W_TO**1.5)

    #Nose gear:
    W_uc2 = 1.08 * (A + B * W_TO**0.75 + C * W_TO + D * W_TO**1.5)

    return W_uc1 + W_uc2

def nacelle_weight(P_TO): #Take off power in hp
    W_n = 1.134 * P_TO**0.5


def equipment_weight(W_TO):
    return 0.08 * W_TO #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF












