from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import Class_II_cg_estimation as c2cg
import geometry_determination as geo
import V_n_diagrams as Vn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import copy

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

# create dataframe with members and values, to save all concepts in
df = pd.DataFrame()

out_folder_path = "output_folder"
os.makedirs(out_folder_path, exist_ok=True)

with open("plotting_data.csv", 'r') as data:
    rows = data.readlines()

for concept in [CON_1, CON_1_braced, CON_2, CON_2_braced, CON_3, CON_3_braced, CON_4, CON_4_braced, CON_5, CON_5_braced]:
    print(concept.name)
    # iteration
    concept_dir = os.path.join(out_folder_path, f"{concept.name}")
    os.makedirs(concept_dir, exist_ok=True)
    W_TO_c2_old = 1500
    m = 0
    for row in rows:
        row_values = row.strip().split(";")
        row_array = np.array(row_values)

        # var_array = np.arange(starting_value, ending value, step) * scale
        var_array = np.arange(float(row_array[1]), float(row_array[2]), float(row_array[3])) * float(row_array[4])
        save = var_array
        if m >= 10:
            var_array = [np.array([i]) for i in var_array]
            print(var_array)
            print(row_array)
        W_TO_array = []
        W_F_array = []
        W_OE_array = []
        m += 1
        n = 1
        for i in var_array:
            concept_analysis = copy.deepcopy(concept)
            var = row_array[0]
            setattr(concept_analysis, var, i)

            it = True
            while it:
                # class 1
                c1.run(concept_analysis)
                
                # geometry determination
                geo.geometry_determination(concept_analysis)

                # class 2
                c2.weight_empty(concept_analysis)
                W_TO_c2 = concept_analysis.W_TO

                # update load factor
                concept_analysis.n_ult = Vn.max_n(concept_analysis)*1.5

                # check if change is small enough
                change = (W_TO_c2 - W_TO_c2_old)/W_TO_c2_old

                if abs(change) < 0.001:
                    W_TO_array.append(float(concept_analysis.W_TO))
                    W_F_array.append(float(concept_analysis.W_F))
                    W_OE_array.append(float(concept_analysis.W_OE))
                    it = False
                else:
                    W_TO_c2_old = W_TO_c2
                    n += 1
                
            del concept_analysis

        print(f"{n} total iterations for {var}")

        # convert to numpy array
        W_TO_array = np.array(W_TO_array)
        W_F_array = np.array(W_F_array)
        W_OE_array = np.array(W_OE_array)

        var_array = save
        # save indices for highlighted part
        subset_indices = np.where((var_array >= float(row_array[5])) & (var_array <= float(row_array[6])))
        var_array_subset = var_array[subset_indices[0]]
        W_TO_range_subset = W_TO_array[subset_indices[0]]
        W_F_range_subset = W_F_array[subset_indices[0]]
        W_OE_range_subset = W_OE_array[subset_indices[0]]

        # plot
        plt.figure(figsize=(20, 5))

        # W_TO plot
        plt.subplot(1, 3, 1)
        plt.plot(var_array, W_TO_array, color="cornflowerblue")
        plt.plot(var_array_subset, W_TO_range_subset, color="red", lw="3")
        plt.xlabel(row_array[7], fontsize=13)
        plt.ylabel("W_TO [kg]", fontsize=13)
        plt.grid()

        # W_F plot
        plt.subplot(1, 3, 2)
        plt.plot(var_array, W_F_array, color="cornflowerblue")
        plt.plot(var_array_subset, W_F_range_subset, color="red", lw='3')
        plt.xlabel(row_array[7], fontsize=13)
        plt.ylabel("W_F [kg]", fontsize=13)
        plt.grid()

        # W_OE plot
        plt.subplot(1, 3, 3)
        plt.plot(var_array, W_OE_array, color="cornflowerblue")
        plt.plot(var_array_subset, W_OE_range_subset, color="red", lw= '3')
        plt.xlabel(row_array[7], fontsize=13)
        plt.ylabel("W_OE [kg]", fontsize=13)
        plt.subplots_adjust(wspace=0.2)
        plt.grid()

        file_name = f"figure_{row_array[0]}_{concept.name}"
        file_path = os.path.join(concept_dir, file_name)
        plt.savefig(file_path, bbox_inches='tight')

        plt.close()
    print("\n")