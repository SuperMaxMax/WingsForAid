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

df = pd.read_csv("airfoil_data_XFLR5.csv")

#df = pd.DataFrame([[5,4,2,8, 2, 3, 2],
                    #[2,1,8,2, 2, 3, 2],
                    #[2,1,2,5, 2, 3, 2]], columns = ["cd0", "cm0", "C", "D", "cl0", "cl_alpha", "cm_alpha" ])

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


def tau_init(df):        #Input must 
    df.drop("cd0", axis = 1)
    df.drop("cl0", axis = 1)
    df.drop("cl_alpha", axis = 1)
    df.drop("cm_alpha", axis = 1)
    df.drop("cm0", axis = 1)


    max_df = df.max()

    # max_df["cd0"] = df["cd0"].min()
    # max_df["cm0"] = df["cm0"].min()

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


def tau_add(df):
    cd0_min = df["cd0"].min()
    cd0_norm = cd0_min / df["cd0"]

    cl0_max = df["cl0"].max()
    cl0_norm = df["cl0"] / cl0_max

    cla_max = df["cl_alpha"].max()
    cla_norm = df["cl_alpha"] / cla_max

    cma_max = df["cm_alpha"].max()
    cma_norm = df["cm_alpha"] / cma_max

    cm0_min = abs(df["cm0"]).min()
    cm0_norm = abs(cm0_min / df["cm0"])

    add = cd0_norm + cl0_norm * (cla_norm + cma_norm + cm0_norm)

    return add

import numpy as np

def tau(df):
    f = np.array(tau_init(df))
    f = np.transpose(f)[0]
    g = np.array(tau_add(df))
    for i in range(len(f)):
        f[i] = f[i] + g[i]
    print(g)  #  f.add(g, fill_value=0)
    return f



def eta(df):     # Weights based on mission profile
    cruise_weight   = 0.475
    loiter_weight   = 0.038
    ascend_weight   = 0.174
    descend_weight  = 0.228
    cruise_par      = "clcd12_max"
    loiter_par      = "clcd32_max"
    ascend_par      = "clcd_max"
    descend_par     = "clcd_max"


    eta = cruise_weight * df[cruise_par] + loiter_weight * df[loiter_par] + ascend_weight * df[ascend_par] + descend_weight * df[descend_par]

    return eta

print("tau = ", tau(df))
print("eta = ", eta(df))
