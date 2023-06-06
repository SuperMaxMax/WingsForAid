from parameters import UAV, Cessna_172
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import Class_II_cg_estimation as c2cg
import geometry_determination as geo
# import V_n_diagrams as Vn
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

aircraft = UAV("droppy", 'tractor', braced_wing=True, boom=True)

# start
plot = False
remove_duplicates = False

# create dataframe with members and values, to save all concepts in
df = pd.DataFrame()

for concept in [aircraft]:
    # --- iteration
    n_boxes = 12
    n = 1
    it = True
    W_TO_c2_old = 750

    while it:
        # class 1
        W_OE = c1.weight_empty_operational(concept)
        W_F  = c1.profile(concept, 12, 1)
        W_PL = n_boxes*concept.boxweight
        concept.W_TO = W_OE + W_F + W_PL
        print(f"Class I estimation {n}: MTOW [kg] = {np.round(concept.W_TO)}")
        concept.W_OE = W_OE
        concept.W_F  = W_F
        
        # geometry determination
        geo.geometry_determination(concept)
        #concept.WS = 70.805

        # class 2
        c2.weight_empty(concept)
        concept.W_TO = concept.W_OE + concept.W_F + n_boxes*concept.boxweight
        W_TO_c2 = concept.W_OE + concept.W_F + n_boxes*concept.boxweight
        print(f"Class II estimation {n}: MTOW = {np.round(concept.W_TO, 2)} [kg]")
        print(f"Fuel weight: {np.round(concept.W_F, 2)} [kg]")
        print(f"Operative empty weight: {np.round(concept.W_OE, 2)} [kg]")

        # update load factor
        # concept.n_ult = Vn.max_n(concept)*1.5

        # check if change is small enough
        change = (W_TO_c2 - W_TO_c2_old)/W_TO_c2_old

        if abs(change) < 0.00001:
            it = False
        else:
            W_TO_c2_old = W_TO_c2
            n += 1

    # --- plotting of concept
    print(f"{concept.name} done in {n} iterations \n")
    geo.geometry_determination(concept, plot = True)
    # cg calculation
    plt.figure(1)
    plt.subplot(121)
    c2cg.cg_calc(concept)

    

    # V-n diagram
    # plt.subplot(122)
    # Vn.plot_all(concept)
    # if plot == True:
    #     plt.show()
    
    # --- saving
    # save all attributes of object to csv file
    members = [attr for attr in dir(concept) if not callable(getattr(concept, attr)) and not attr.startswith("__")]
    values = [getattr(concept, member) for member in members]

    # remove brackets and round values
    # values = [value[0] if isinstance(value, np.ndarray) else value for value in values]
    # values = [round(value, 4) if isinstance(value, float) else value for value in values]

    # add to dataframe
    df[concept.name] = values

# set index of dataframe
df.index = members

# export dataframe of current design to csv file
df['droppy'].to_csv('output.csv', sep=';')

# remove row in dataframe if all values in that row are the same
if remove_duplicates == True:
    for i in df.index:
        if all(element == df.loc[i].values[0] for element in df.loc[i].values):
            df.drop(i, inplace=True)
        
# save dataframe to csv file
df.to_csv('concept_comparison.csv', sep=';')