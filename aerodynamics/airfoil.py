"clcd_max"
"clcd32_max"
"clcd12_max"
"clmax"
"alpha_s"
"cd0"
"cl0"
"cl_alpha"
"cm_alpha"
"cm0"

import pandas as pd

df = pd.DataFrame([[5,4,2,8],
                    [2,1,8,2],
                    [2,1,2,5]], columns = ["A", "B", "C", "D"])

# print(df)
# def invert(df, parameter):
#     df[parameter] = df[parameter].apply(lambda x: 1/x)
#     return 

# def combine(df, parameter1, parameter2, newname):
#     new = df[parameter1] * df[parameter2]
#     df[newname] = new
#     return 



def remove(df, parameter):
    del df[parameter]
    return 


def tau_df(df):    # Reshape dataframe to workable values for tau,
                              # combining parameters that need to be combined and 
                              # inverting parameters where a maximum is desired.
                              # Assign weights if necessary
    return


def tau(df):        #Input must 
    max_df = df.max()
    taulist = []

    

    for i in range(len(df)):
        row = df.iloc[[i]] #Row of an airfoil 
        norm_list = []
        for j in range(len(row.columns)): 
            norm = row.iloc[:, j]/max_df[j]
            norm_list.append(norm)
        tau = sum(norm_list)
        taulist.append(tau)
    return taulist



def eta(df):     # Weights based on mission profile
    cruise_weight   = 0.475
    loiter_weight   = 0.038
    ascend_weight   = 0.174
    descend_weight  = 0.228
    cruise_par      = "clcd12_max"
    loiter_par      = "clcd32_max"
    ascend_par      = "clcd_max"
    descend_par     = "clcd_max"


    eta = cruise_weight * df[cruise_par] + loiter_weight * [loiter_par] + ascend_weight * [ascend_par] + descend_weight * [descend_par]

    return eta




