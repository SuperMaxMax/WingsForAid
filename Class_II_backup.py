def wing_weight(b, lambda_mid, n_ult, t_c, cwr, W_loading, W_G):   #Span, sweep, ultimate load factor, thickness over chord, cord length at root, wing loading, gross weight
    k_w = 4.9*10**(-3)
    b_s = b / np.cos(lambda_mid)
    b_ref = 1.905
    t_r = t_c * cwr

    W_w = k_w * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * n_ult**0.55 * ((b_s/t_r)/W_loading)**0.3 * W_G
    
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
    W_uc1 = 1.08 * (A1 + B1 * W_TO**0.75 + C1 * W_TO + D1 * W_TO**1.5)

    #Nose gear:
    W_uc2 = 1.08 * (A2 + B2 * W_TO**0.75 + C2 * W_TO + D2 * W_TO**1.5)

    return W_uc1 + W_uc2

def nacelle_weight(P_TO): #Take off power in hp
    W_n = 1.134 * P_TO**0.5


def equipment_weight(W_TO):
    return 0.08 * W_TO #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF

def fuselage_weight(V_D, l_t, b_f, h_f, S_G):
    k_wf = 0.23
    W_f = k_wf*(V_D*(l_t/(b_f+h_f)))**0.5*S_G**1.2

    # Add 7% if the main landing gear is attached to the fuselage, but 4% can be subtracted from the fuselage weight if there is no attachment structure for the main landing gear
    # Add 10% for freighter aircraft
    # For booms l_t is defined as distance between local wing chord and horizontal tailplane

    return W_f

def control_surface_weight(W_to):
    k_sc = 0.44*0.768
    W_sc = k_sc*W_to**(2/3)

    # Add 20% for LE flap or slat
    # Add 15% for lift dumper controls

    return W_sc

def propulsion_weight(N_e, W_e, P_to):
    k_pg = 1.16 # tractor single propeller aircraft
    W_pg = k_pg*N_e*(W_e+0.109*P_to)

    # If number of cylinder and volume of cylinder are known use figure 4-12 Torenbeek

    return W_pg

# def oil_and_residual_fuel_weight(W_to, W_f):
#     W_orf = 0.008*W_to+0.045*W_f
