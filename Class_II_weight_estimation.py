import numpy as np



def wing_weight(b, lambda_mid, n_ult, t_c, cwr, W_loading, W_G):
    k_w = 4.9*10**(-3)
    b_s = b / np.cos(lambda_mid)
    b_ref = 1.905
    t_r = t_c * cwr

    W_w = kw * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * n_ult**0.55 * ((b_s/t_r)/W_loading)**0.3 * W_G
    
    #ADD 30% IF BRACED WING USED, 10% IF STRUT USED?
    
    return W_w
    
def tail_weight(k_wt, n_ult, s_tail):
    W_t = 0.64 * (n_ult * s_tail**2)**0.75

    #IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT
    return W_t





