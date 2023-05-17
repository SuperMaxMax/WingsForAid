from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import Class_II_cg_estimation as c2cg
import geometry_determination as geo
import V_n_diagrams as Vn
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

CON_1 = UAV('CON_1', 'tractor', boom=True, braced_wing=False)
CON_1_braced = UAV('CON_1_braced', 'tractor', boom=True, braced_wing=True)
CON_2 = UAV('CON_2', 'tractor', boom=False, braced_wing=False)
CON_2_braced = UAV('CON_2_braced', 'tractor', boom=False, braced_wing=True)
CON_3 = UAV('CON_3', 'pusher', boom=False, braced_wing=False)
CON_3_braced = UAV('CON_3_braced', 'pusher', boom=False, braced_wing=True)
CON_4 = UAV('CON_4', 'pusher', boom=False, braced_wing=False)
CON_4_braced = UAV('CON_4_braced', 'pusher', boom=False, braced_wing=True)
CON_5 = UAV('CON_5', 'fuselage', boom=False, braced_wing=False)
CON_5_braced = UAV('CON_5_braced', 'fuselage', boom=False, braced_wing=True)

# start
plot = True
remove_duplicates = False

# create dataframe with members and values, to save all concepts in
df = pd.DataFrame()

for concept in [CON_1, CON_1_braced, CON_2, CON_2_braced, CON_3, CON_3_braced, CON_4, CON_4_braced, CON_5, CON_5_braced]:
    # --- iteration
    n = 1
    it = True
    W_TO_c2_old = 750

    while it:
        # class 1
        print(f"- Iteration number: {n}, concept: {concept.name} - \n")
        c1.run(concept)
        W_TO_c1 = concept.W_TO
        
        # geometry determination
        geo.geometry_determination(concept)

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
    # cg calculation
    plt.figure(1)
    plt.subplot(121)
    c2cg.cg_calc(concept, plot=False)

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

# remove row in dataframe if all values in that row are the same
if remove_duplicates == True:
    for i in df.index:
        if all(element == df.loc[i].values[0] for element in df.loc[i].values):
            df.drop(i, inplace=True)
        
# save dataframe to csv file
df.to_csv('concept_comparison.csv', sep=';')