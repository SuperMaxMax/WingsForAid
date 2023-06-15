from parameters import UAV
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import aerodynamics.main_aero as ae
import control_stability.main_stab_cont as cs
import structures.wingbox_full as wb

aircraft = UAV("aircraft")

# start
plot = False
remove_duplicates = False

ae.run_aero(aircraft)

cs.main_stab_control(aircraft, True, False)                     # FIXME: Tomorrow ask Theo about W_eq and calculate W_sc and W_tail

# print all attributes of object
print(aircraft.__dict__)

# wb.all_wingbox(aircraft)

# # create dataframe with members and values, to save all aircrafts in
# df = pd.DataFrame()

# # --- iteration
# n = 1
# it = True
# W_TO_c2_old = 750

# while it:
#     ae.wp.main_wing_planform(aircraft)
#     ae.htd.horizontal_tail_planform(aircraft)
#     ae.vtd.horizontal_tail_planform(aircraft)

#     # check if change is small enough
#     change = (W_TO_c2 - W_TO_c2_old)/W_TO_c2_old

#     if np.abs(change) < 0.001:
#         it = False
#     else:
#         n += 1
#         W_TO_c2_old = W_TO_c2
    
# if plot == True:
#     plt.show()

# # --- saving
# # save all attributes of object to csv file
# members = [attr for attr in dir(aircraft) if not callable(getattr(aircraft, attr)) and not attr.startswith("__")]
# values = [getattr(aircraft, member) for member in members]

# # remove brackets and round values
# values = [value[0] if isinstance(value, np.ndarray) else value for value in values]
# values = [round(value, 4) if isinstance(value, float) else value for value in values]

# # add to dataframe
# df[aircraft.name] = values

# # set index of dataframe
# df.index = members

# # export dataframe of current design to csv file
# df['DET_CON_2_braced'].to_csv('DET_CON_2_braced.csv', sep=';')

# # remove row in dataframe if all values in that row are the same
# if remove_duplicates == True:
#     for i in df.index:
#         if all(element == df.loc[i].values[0] for element in df.loc[i].values):
#             df.drop(i, inplace=True)
        
# # save dataframe to csv file
# df.to_csv('aircraft_comparison.csv', sep=';')