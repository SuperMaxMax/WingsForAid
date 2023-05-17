from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import Class_II_cg_estimation as c2cg
import geometry_determination as geo
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

concept_1 = UAV('concept_1', 'tractor', boom=True)
concept_2 = UAV('concept_2', 'tractor', boom=False)
concept_3 = UAV('concept_3', 'pusher', boom=False)
concept_4 = UAV('concept_4', 'pusher', boom=False)
concept_5 = UAV('concept_5', 'fuselage', boom=False)

# create dataframe with members and values, to save all concepts in
df = pd.DataFrame()

out_folder_path = "output_folder"
os.makedirs(out_folder_path, exist_ok=True)

for concept in [concept_1, concept_2, concept_3, concept_4, concept_5]:
    # iteration
    concept_dir = os.path.join(out_folder_path, f"{concept.name}")
    os.makedirs(concept_dir, exist_ok=True)
    n = 1
    W_TO_c2_old = 750
    with open("plotting_data.csv", 'r') as data:
        rows = data.readlines()

        for row in rows:
            row_values = row.strip().split(";")

            # Convert the row values to an array
            row_array = np.array(row_values)
            #var_array = np.arange(starting_value, ending value, step) * scale
            var_array = np.arange(float(row_array[1]), float(row_array[2]), float(row_array[3])) * float(row_array[4])

        # var_array = np.arange(4, 20, 0.5)
        # var_array = np.arange(0.5, 1, 0.02)
        # var_array = np.arange(0.01, 0.06, 0.001)
        # var_array = np.arange(0.5, 1, 0.02)
        # var_array = np.arange(1, 7, 1)
        # var_array = np.arange(50, 240, 1) * (1.852/3.6)
        # var_array = np.arange(1.3, 1.9, 0.05)
        # var_array = np.arange(5000, 18000, 50) * 0.3048
        # var_array = [np.array([i]) for i in var_array]
        # var_array = np.arange(60, 90, 0.02) * (1.852/3.6)
        # var_array = np.arange(1, 7, 0.02)
        # var_array = np.arange(50000, 600000, 1000)
            W_TO_array = []
            W_F_array = []
            W_OE_array = []
            for i in var_array:
                concept_analysis = concept

                var = row_array[0]

                setattr(concept_analysis, var, i)
                # concept.A = i

                # concept.e = i
                # concept.CD0 = i
                # concept.prop_eff = i
                # concept.n_drops = i
                # concept.V_cruise = i
                # concept.CL_LDG = i
                # concept.h_cruise = i
                # concept.V_climb = i
                # concept.n_ult = i
                # concept.R = i
                it = True
                while it:
                    # class 1
                    print(f"- Iteration number: {n}, concept: {concept.name} - \n")
                    c1.run(concept_analysis)
                    print("After class 1 iteration:")
                    print(f"Mff: {concept_analysis.Mff}, L/D: {concept_analysis.L_D}, W_TO: {concept_analysis.W_TO}, W_OE: {concept_analysis.W_OE}, W_F: {concept_analysis.W_F}")
                    W_TO_c1 = concept_analysis.W_TO

                    # geometry determination
                    geo.geometry_determination(concept_analysis)
                    print("After geometry determination:")
                    print(f"b: {concept_analysis.b} Sw: {concept_analysis.Sw} S_G: {concept_analysis.S_G}")

                    # class 2
                    c2.weight_empty(concept_analysis)
                    print("After class 2 iteration:")
                    print(f"W_OE: {concept_analysis.W_OE}, W_TO: {concept_analysis.W_TO}")
                    W_TO_c2 = concept_analysis.W_TO

                    print("")
                    print("-------------------------------------------------")
                    print("")

                    # iterate between class 1 and class 2
                    change = (W_TO_c2 - W_TO_c2_old)/W_TO_c2_old

                    if abs(change) < 0.00001:
                        W_TO_array.append(float(concept_analysis.W_TO))
                        W_F_array.append(float(concept_analysis.W_F))
                        W_OE_array.append(float(concept_analysis.W_OE))
                        it = False
                    else:
                        W_TO_c2_old = W_TO_c2
                        n += 1

            W_TO_array = np.array(W_TO_array)
            W_F_array = np.array(W_F_array)
            W_OE_array = np.array(W_OE_array)

            subset_indices = np.where((var_array >= float(row_array[5])) & (var_array <= float(row_array[6])))

            var_array_subset = var_array[subset_indices[0]]

            W_TO_range_subset = W_TO_array[subset_indices[0]]
            W_F_range_subset = W_F_array[subset_indices[0]]
            W_OE_range_subset = W_OE_array[subset_indices[0]]

            "W_TO PLOT"
            plt.figure(figsize=(20, 5))

            plt.subplot(1, 3, 1)
            plt.plot(var_array, W_TO_array, color="cornflowerblue")

            plt.plot(var_array_subset, W_TO_range_subset, color="red")

            plt.xlabel(row_array[7])
            plt.ylabel("W_TO [kg]")
            plt.grid()


            "W_F PLOT"

            plt.subplot(1, 3, 2)
            plt.plot(var_array, W_F_array, color="cornflowerblue")

            plt.plot(var_array_subset, W_F_range_subset, color="red")

            plt.xlabel(row_array[7])
            plt.ylabel("W_F [kg]")
            plt.grid()


            "W_OE PLOT"
            plt.subplot(1, 3, 3)
            plt.plot(var_array, W_OE_array, color="cornflowerblue")

            plt.plot(var_array_subset, W_OE_range_subset, color="red")

            plt.xlabel(row_array[7])
            plt.ylabel("W_OE [kg]")

            plt.subplots_adjust(wspace=0.2)
            plt.grid()

            #plt.show()

            file_name = f"figure_{row_array[0]}_{concept.name}"
            file_path = os.path.join(concept_dir, file_name)
            plt.savefig(file_path)

            plt.close()


#     # cg calculation
#     c2cg.cg_calc(concept, plot=False)
#
#     # save all attributes of object to csv file
#     members = [attr for attr in dir(concept) if not callable(getattr(concept, attr)) and not attr.startswith("__")]
#     values = [getattr(concept, member) for member in members]
#
#     # round values
#     values = [round(value, 3) if isinstance(value, float) else value for value in values]
#
#     # add to dataframe
#     df[concept.name] = values
#
# # set index of dataframe
# df.index = members
#
# # save dataframe to csv file
# df.to_csv('output.csv', sep=';')