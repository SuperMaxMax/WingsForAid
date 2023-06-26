from parameters import UAV_final, atmosphere, airport
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import aerodynamics.main_aero as ae
import control_stability.main_stab_cont as cs
import structures.wingbox_full as wb
import structures.fuselage_fairing_buckling as ffb
import flight_performance.simulation as fp
import operations.sortie as op

# aircraft = UAV_final()
aircraft = UAV_final()
atm      = atmosphere()
airfield = airport("Sudan")

index_df = ["MTOW", "OEW", "Fuel weight", "Tail weight", "Wing weight", "Fuselage weight"]
df_iterations = pd.DataFrame(index = index_df)

# start
plot = True
jan = False
theo = False

n_iteration = 1
running = True

# for wingbox
jarno = False               # Jarno can't run wingbox
full_wingbox_loop = False    # if true the full optimization will run, else it will calculate the weight for 36 stringers, 17 ribs
ele_span = 500              # number of elements in spanwise direction (smaller value is faster, but less accurate)

#TODO landing distance
#TODO max ferry range
#TODO ceiling altitude
#TODO max climb rate
#TODO max endurance

while running:
    W_check = aircraft.W_OE + aircraft.W_F

    print(f"=================== AERO-{n_iteration} =====================")
    ae.run_aero(aircraft)
    print("================================================\n")

    print(f"=================== CS-{n_iteration} =======================")
    cs.main_stab_control(aircraft, plot, True, False)
    print("================================================\n")

    print(f"=================== WB-{n_iteration} =======================")
    if not jarno: #Jarno can't run wingbox
        if full_wingbox_loop and n_iteration == 1:
            wb.all_wingbox(aircraft, ele_span, True)
        elif full_wingbox_loop and n_iteration > 1:
            wb.all_wingbox(aircraft, ele_span, False, True)
        else:
            wb.all_wingbox(aircraft, ele_span, False)
    print("================================================\n")

    print(f"=================== FFB-{n_iteration} =======================")
    ffb.fuselage_fairing(aircraft)
    print("================================================\n")
    
    print(f"=================== FP-{n_iteration} =======================")
    aircraft.W_F = fp.fuelusesortie(aircraft, atm, 12, 1, 10000, aircraft.W_F, 54.012, Summary = True)[0] + 5
    fp.TO_eom(aircraft, airfield, atm, 11, 4000, 12.86, -7.716, aircraft.W_F, Plot = False)
    fp.LA_eom(aircraft, airfield, atm, -8, 4000, 12.86, -5.14, 5, Plot = False)
    print("=================================================\n")

    print(f"=================== MANUAL UPDATES-{n_iteration} =======================")
    if jan: #Jan's path is linked in avl so otherwise code breaks
        import aerodynamics.avl as avl
        avl.export(aircraft)

    # aircraft.CL_max_clean = float(input("Wing CL_max_clean: "))
    print(f"=================================================\n")

    aircraft.W_OE = aircraft.W_eq + aircraft.W_n + aircraft.W_pg + aircraft.W_sc + aircraft.W_t + aircraft.W_strut + aircraft.ST_W_fus + aircraft.ST_W_boom + aircraft.ST_W_uc + aircraft.W_w
    aircraft.W_TO = aircraft.W_F + aircraft.W_OE + aircraft.W_PL

    if np.abs((aircraft.W_OE + aircraft.W_F - W_check)/W_check) < 0.001:
        running = False

    aircraft.Sw = ((aircraft.W_OE + aircraft.W_F + aircraft.n_boxes*aircraft.boxweight)*atm.g)/aircraft.WS
    aircraft.b = np.sqrt(aircraft.A*aircraft.Sw)

    df_iterations[f"Iteration {n_iteration}"] = [aircraft.W_TO, aircraft.W_OE, aircraft.W_F, aircraft.W_t, aircraft.W_w, aircraft.ST_W_fus]

    print(f"================= GENERAL-INFO ==================")
    print(df_iterations)
    print("===================================================\n")

    n_iteration += 1

if theo:
    fp.fuelusesortie(aircraft, atm, aircraft.n_boxes, 1, aircraft.h_cruise / 0.3048, aircraft.W_F, 54.012,
                        aircraft.OP_Range, Summary=True, plot = True)            
    op.operations_eval(aircraft)


# # --- saving
df = pd.DataFrame()
# save all attributes of object to csv file
members = [attr for attr in dir(aircraft) if not callable(getattr(aircraft, attr)) and not attr.startswith("__")]
values = [getattr(aircraft, member) for member in members]

# remove brackets and round values
values = [list(value) if isinstance(value, np.ndarray) else value for value in values]
# values = [round(value, 4) if isinstance(value, float) else value for value in values]

# add to dataframe
df[aircraft.name] = values

# set index of dataframe
df.index = members

# export dataframe of current design to csv file
df['aircraft'].to_csv('finaldesign.csv', sep=';')