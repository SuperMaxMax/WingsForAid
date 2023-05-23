from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import Class_II_cg_estimation as c2cg
import geometry_determination as geo
import V_n_diagrams as Vn
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

DET_CON_1 = UAV('DET_CON_1', 'tractor', boom=False, braced_wing=False)
DET_CON_1_braced = UAV('DET_CON_1_braced', 'tractor', boom=False, braced_wing=True)
DET_CON_2 = UAV('DET_CON_2', 'tractor', boom=True, braced_wing=False)
DET_CON_2_braced = UAV('DET_CON_2_braced', 'tractor', boom=True, braced_wing=True)
DET_CON_3 = UAV('DET_CON_3', 'pusher', boom=False, braced_wing=False)
DET_CON_3_braced = UAV('DET_CON_3_braced', 'pusher', boom=False, braced_wing=True)
DET_CON_4 = UAV('DET_CON_4', 'pusher', boom=True, braced_wing=False)
DET_CON_4_braced = UAV('DET_CON_4_braced', 'pusher', boom=True, braced_wing=True)
DET_CON_5 = UAV('DET_CON_5', 'fuselage', boom=False, braced_wing=False)
DET_CON_5_braced = UAV('DET_CON_5_braced', 'fuselage', boom=False, braced_wing=True)
DET_CON_6 = UAV('DET_CON_6', 'fuselage', boom=True, braced_wing=False)
DET_CON_6_braced = UAV('DET_CON_6_braced', 'fuselage', boom=True, braced_wing=True)

# start
plot = False
remove_duplicates = False

# create dataframe with members and values, to save all concepts in
df = pd.DataFrame()

for concept in [DET_CON_1, DET_CON_1_braced, DET_CON_2, DET_CON_2_braced, DET_CON_3, DET_CON_3_braced, DET_CON_4, DET_CON_4_braced, DET_CON_5, DET_CON_5_braced, DET_CON_6, DET_CON_6_braced]:
    # --- iteration
    n = 1
    it = True
    W_TO_c2_old = 750

    while it:
        # class 1
        c1.run(concept)
        
        # geometry determination
        geo.geometry_determination(concept)
       # concept.WS = 70.805

        # class 2
        c2.weight_empty(concept)
        W_TO_c2 = concept.W_TO

        # update load factor
        concept.n_ult = Vn.max_n(concept)*1.5

        # check if change is small enough
        change = (W_TO_c2 - W_TO_c2_old)/W_TO_c2_old

        if abs(change) < 0.00001:
            it = False
        else:
            W_TO_c2_old = W_TO_c2
            n += 1

    # --- plotting of concept
    print(f"{concept.name} done in {n} iterations \n")
    # cg calculation
    plt.figure(1)
    plt.subplot(121)
    c2cg.cg_calc(concept)

    # V-n diagram
    plt.subplot(122)
    Vn.plot_all(concept)
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
df['DET_CON_1_braced'].to_csv('output.csv', sep=';')

# remove row in dataframe if all values in that row are the same
if remove_duplicates == True:
    for i in df.index:
        if all(element == df.loc[i].values[0] for element in df.loc[i].values):
            df.drop(i, inplace=True)
        
# save dataframe to csv file
df.to_csv('concept_comparison.csv', sep=';')