from parameters import UAV, atmosphere, airport
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import aerodynamics.main_aero as ae
import control_stability.main_stab_cont as cs
import structures.wingbox_full as wb
import flight_performance.simulation as fp


aircraft = UAV("aircraft")
atm      = atmosphere()
airfield = airport("Sudan")

# start
plot = False
remove_duplicates = False
jan = False
jarno = False

n_iteration = 0
something = True

while something:
    print('Iteration: ', n_iteration)
    print(f'MTOW: {aircraft.W_TO:.2f} kg, OEW: {aircraft.W_OE:.2f}')
    W_TO_old = aircraft.W_TO

    ae.run_aero(aircraft)
    cs.main_stab_control(aircraft, False, False) # FIXME: Tomorrow ask Theo about W_eq and calculate W_sc and W_tail
    print(aircraft.n_stringers, aircraft.n_ribs, aircraft.W_w)

    if not jarno:
        if n_iteration % 1  == 0 and n_iteration != 0:
            wb.all_wingbox(aircraft, True)
        else:
            wb.all_wingbox(aircraft, False)
    
    aircraft.W_F = fp.fuelusesortie(aircraft, atm, 12, 1, 10000, aircraft.W_F, 54.012, Summary = True)[0] + 5
    fp.TO_eom(aircraft, airfield, atm, 11, 4000, 12.86, -7.716, aircraft.W_F, Plot = False)
    fp.LA_eom(aircraft, airfield, atm, -8, 4000, 12.86, -5.14, 5, Plot = False)

    if np.abs((aircraft.W_TO - W_TO_old)/W_TO_old) < 0.001:
        something = False

    print(aircraft.__dict__)

    print("=================================")
    print("TAKE-OFF WEIGHT:", aircraft.W_TO)
    print("=================================")
    n_iteration += 1
    
if jan:
    import aerodynamics.avl as avl
    avl

    
# print all attributes of object
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