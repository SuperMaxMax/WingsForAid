from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import Class_II_cg_estimation as c2cg
import geometry_determination as geo
import V_n_diagrams as Vn
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from parameters import UAV
import aerodynamics.main_aero as ae
import 



aircraft = UAV("aircraft")

# start
plot = False
remove_duplicates = False

# create dataframe with members and values, to save all concepts in
df = pd.DataFrame()

for concept in [aircraft]:
    # --- iteration
    n = 1
    it = True
    W_TO_c2_old = 750
    

    while it:
        
        ae.wp.main_wing_planform(aircraft)
        ae.htd.horizontal_tail_planform(aircraft)
        ae.vtd.horizontal_tail_planform(aircraft)
       

        # check if change is small enough
        change = (W_TO_c2 - W_TO_c2_old)/W_TO_c2_old

        if np.abs(change) < 0.001:
            it = False
        else:
            n += 1
            W_TO_c2_old = W_TO_c2
        
    if plot == True:
        plt.show()
    
    # --- saving
    # save all attributes of object to csv file
    members = [attr for attr in dir(concept) if not callable(getattr(concept, attr)) and not attr.startswith("__")]
    values = [getattr(concept, member) for member in members]

    # remove brackets and round values
    values = [value[0] if isinstance(value, np.ndarray) else value for value in values]
    values = [round(value, 4) if isinstance(value, float) else value for value in values]

    # add to dataframe
    df[concept.name] = values

# set index of dataframe
df.index = members

# export dataframe of current design to csv file
df['DET_CON_2_braced'].to_csv('DET_CON_2_braced.csv', sep=';')

# remove row in dataframe if all values in that row are the same
if remove_duplicates == True:
    for i in df.index:
        if all(element == df.loc[i].values[0] for element in df.loc[i].values):
            df.drop(i, inplace=True)
        
# save dataframe to csv file
df.to_csv('concept_comparison.csv', sep=';')

cs.





