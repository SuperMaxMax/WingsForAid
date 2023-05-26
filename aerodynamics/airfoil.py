# clcd_max
# clcd32_max
# clcd12_max
# clmax
# alpha_s
# cd0
# cl0
# cl_alpha
# cm_alpha
# cm0

import pandas as pd

df = pd.DataFrame([[5,4,2,8],
                  [2,1,8,2],
                  [2,1,2,5]])


def tau_df(df):    # Reshape dataframe to workable values for tau,
                        # combining parameters that need to be combined and 
                        # inverting parameters where a maximum is desired.
                        # Assign weights if necessary
    return

def eta_df(df):     # Weights based on mission profile
    cruise_weight   =
    loiter_weight   =
    ascend_weight   =
    descend_weight  =
    cruise_par      =
    loiter_par      =
    ascend_par      =
    descend_par     =

    eta = cruise_weight * cruise_par + loiter_weight * loiter_par + ascend_weight * ascend_par + descend_weight * descend_par

    return eta




def tau(df):        #Input must 
    max_df = df.max()
    taulist = []
    for i in range(len(df)):
        row = df.iloc[[i]] #Row of an airfoil
        norm_list = []
        for j in range(len(row.columns)):
            norm = row[j]/max_df[j]
            norm_list.append(norm)
        tau = sum(norm_list)
        taulist.append(tau)
    return taulist

print(tau(df))


def eta(df):

    return