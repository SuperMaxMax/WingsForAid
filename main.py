from parameters import UAV, atmosphere, airport
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import aerodynamics.main_aero as ae
import control_stability.main_stab_cont as cs
import structures.wingbox_full as wb
import flight_performance.simulation as fp
import operations.sortie as op


aircraft = UAV("aircraft")
atm      = atmosphere()
airfield = airport("Sudan")

# start
plot = False
jan = False
jarno = False

n_iteration = 0
something = True
full_wingbox_loop = True

#TODO landing distance
#TODO max ferry range
#TODO ceiling altitude
#TODO max climb rate
#TODO max endurance

while something:
    print('Iteration: ', n_iteration)
    print(f'MTOW: {aircraft.W_TO:.2f} kg, OEW: {aircraft.W_OE:.2f}')
    W_check= aircraft.W_OE + aircraft.W_F
    # aircraft.Sw = aircraft.W_TO/aircraft.WS
    ae.run_aero(aircraft)
    cs.main_stab_control(aircraft, True, False) # FIXME: Tomorrow ask Theo about W_eq and calculate W_sc and W_tail
    print(aircraft.n_stringers, aircraft.n_ribs, aircraft.W_w)

    if not jarno: #Jarno can't run wingbox
        if full_wingbox_loop: # and n_iteration != 0:
            wb.all_wingbox(aircraft, True)
        else:
            wb.all_wingbox(aircraft, False)
    

    aircraft.W_F = fp.fuelusesortie(aircraft, atm, 12, 1, 10000, aircraft.W_F, 54.012, Summary = True)[0] + 5
    fp.TO_eom(aircraft, airfield, atm, 11, 4000, 12.86, -7.716, aircraft.W_F, Plot = False)
    fp.LA_eom(aircraft, airfield, atm, -8, 4000, 12.86, -5.14, 5, Plot = False)

    # op.operations_eval(aircraft)

    if np.abs((aircraft.W_OE + aircraft.W_F - W_check)/W_check) < 0.001:
        something = False

    print(aircraft.__dict__)
    aircraft.__dict__

    print("=================================")
    print("TAKE-OFF WEIGHT:", aircraft.W_TO)
    print("=================================")
    n_iteration += 1

if jan: #Jan's path is linked in avl so otherwise code breaks
    import aerodynamics.avl as avl
    avl.export(aircraft)


# # --- saving
df = pd.DataFrame()
# save all attributes of object to csv file
members = [attr for attr in dir(aircraft) if not callable(getattr(aircraft, attr)) and not attr.startswith("__")]
values = [getattr(aircraft, member) for member in members]

# remove brackets and round values
values = [value[0] if isinstance(value, np.ndarray) else value for value in values]
values = [round(value, 4) if isinstance(value, float) else value for value in values]

# add to dataframe
df[aircraft.name] = values

# set index of dataframe
df.index = members

# export dataframe of current design to csv file
df['aircraft'].to_csv('finaldesign.csv', sep=';')