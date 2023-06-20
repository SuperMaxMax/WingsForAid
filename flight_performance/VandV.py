# ========== Set working directory ==========
import sys
import os.path
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

# ============ Import packages ==============
import numpy as np
import matplotlib.pyplot as plt

# ========== Import other files =============
import flight_performance.simulation as fp
from parameters import UAV, atmosphere, airport, UAV_final, Cessna_172

# ============= Define objects ==============
aircraft = UAV_final()
atm      = atmosphere()
Cessna   = Cessna_172("tractor", braced_wing=True, boom=False)

# ========= V&V sortie simulation ===========
def VV_sortiesim_expected(ac_obj, atm_obj, plot = False):
    # Expected value tests for sortie simulation: varying range, n_drops, dropregion, cruise speed, cruise height, n_boxes
    ranges  = np.arange(50.0, 250.0, 25.0)                      # Default value: 250 km
    drops   = np.arange(1, 6, 1)                                # Default value: 1
    dropregions = np.arange(10.0, 100.0, 10.0)                  # Default value: None
    cruise_V    = np.arange(40.0, 55.0, 2.5)                    # Default value: 54.012 [m/s] / 55 [kts]
    cruise_h    = np.arange(5000.0, 15000.0, 1000.0) * 0.3048   # Default value: 10000 [ft]
    n_boxes     = np.arange(2, 12, 2)                           # Default value: 12
    vary_param  = [ranges, drops, dropregions, cruise_V, cruise_h, n_boxes]
    # Define empty arrays to store results in
    fuel_used = []
    sort_time = []
    for i in range(len(vary_param)):
        fuel_used_i = np.empty(0)
        sort_time_i = np.empty(0)
        for j in range(len(vary_param[i])):
            if i == 0:
                sortie = fp.fuelusesortie(ac_obj, atm_obj, 12, 1, 10000, 50, 54.012, Range = vary_param[i][j])
                fuel_used_i = np.append(fuel_used_i, sortie[0])
                sort_time_i = np.append(sort_time_i, sortie[2])
            if i == 1:
                sortie = fp.fuelusesortie(ac_obj, atm_obj, 12, vary_param[i][j], 10000, 50, 54.012)
                fuel_used_i = np.append(fuel_used_i, sortie[0])
                sort_time_i = np.append(sort_time_i, sortie[2])
            if i == 2:
                sortie = fp.fuelusesortie(ac_obj, atm_obj, 12, 1, 10000, 50, 54.012, dropregion=vary_param[i][j])
                fuel_used_i = np.append(fuel_used_i, sortie[0])
                sort_time_i = np.append(sort_time_i, sortie[2])
            if i == 3:
                sortie = fp.fuelusesortie(ac_obj, atm_obj, 12, 1, 10000, 50, vary_param[i][j])
                fuel_used_i = np.append(fuel_used_i, sortie[0])
                sort_time_i = np.append(sort_time_i, sortie[2])
            if i == 4:
                sortie = fp.fuelusesortie(ac_obj, atm_obj, 12, 1, vary_param[i][j], 50, 54.012)
                fuel_used_i = np.append(fuel_used_i, sortie[0])
                sort_time_i = np.append(sort_time_i, sortie[2])
            if i == 5:
                sortie = fp.fuelusesortie(ac_obj, atm_obj, vary_param[i][j], 1, 10000, 50, 54.012)
                fuel_used_i = np.append(fuel_used_i, sortie[0])
                sort_time_i = np.append(sort_time_i, sortie[2])
        fuel_used.append(fuel_used_i)
        sort_time.append(sort_time_i)
    if plot:
        fig1, axes1 = plt.subplots(nrows = 1, ncols = 2)
        axes1[0].plot(ranges, fuel_used[0])
        axes1[0].set_xlabel("Range [km]")
        axes1[0].set_ylabel("Fuel used [kg]")
        axes1[1].plot(ranges, sort_time[0])
        axes1[1].set_xlabel("Range [km]")
        axes1[1].set_ylabel("Sortie time [sec]")
        fig1.suptitle("Fuel used and Sortie time versus Range")
        fig1.tight_layout()
        plt.savefig("C:\\Users\\ties\\Downloads\\Sensitivity_range")
        plt.show()
        fig2, axes2 = plt.subplots(nrows = 1, ncols = 2)
        axes2[0].plot(drops, fuel_used[1])
        axes2[0].set_xlabel("Drops")
        axes2[0].set_ylabel("Fuel used [kg]")
        axes2[1].plot(drops, sort_time[1])
        axes2[1].set_xlabel("Drops")
        axes2[1].set_ylabel("Sortie time [sec]")
        fig2.suptitle("Fuel used and Sortie time versus number of drops")
        fig2.tight_layout()
        plt.savefig("C:\\Users\\ties\\Downloads\\Sensitivity_n_drops")
        plt.show()
        fig3, axes3 = plt.subplots(nrows = 1, ncols = 2)
        axes3[0].plot(dropregions, fuel_used[2])
        axes3[0].set_xlabel("Dropregion size [km]")
        axes3[0].set_ylabel("Fuel used [kg]")
        axes3[1].plot(dropregions, sort_time[2])
        axes3[1].set_xlabel("Dropregion size [km]")
        axes3[1].set_ylabel("Sortie time [sec]")
        fig3.suptitle("Fuel used and Sortie time versus dropregion size")
        fig3.tight_layout()
        plt.savefig("C:\\Users\\ties\\Downloads\\Sensitivity_dropregion_size")
        plt.show()
        fig4, axes4 = plt.subplots(nrows = 1, ncols = 2)
        axes4[0].plot(cruise_V, fuel_used[3])
        axes4[0].set_xlabel("Cruise speed [m/s]")
        axes4[0].set_ylabel("Fuel used [kg]")
        axes4[1].plot(cruise_V, sort_time[3])
        axes4[1].set_xlabel("Cruise speed [m/s]")
        axes4[1].set_ylabel("Sortie time [sec]")
        fig4.suptitle("Fuel used and Sortie time versus cruise speed")
        fig4.tight_layout()
        plt.savefig("C:\\Users\\ties\\Downloads\\Sensitivity_cruisespeed")
        plt.show()
        fig5, axes5 = plt.subplots(nrows = 1, ncols = 2)
        axes5[0].plot(cruise_h, fuel_used[4])
        axes5[0].set_xlabel("Cruise height [m]")
        axes5[0].set_ylabel("Fuel used [kg]")
        axes5[1].plot(cruise_h, sort_time[4])
        axes5[1].set_xlabel("Cruise height [m]")
        axes5[1].set_ylabel("Sortie time [sec]")
        fig5.suptitle("Fuel used and Sortie time versus cruise height")
        fig5.tight_layout()
        plt.savefig("C:\\Users\\ties\\Downloads\\Sensitivity_cruiseheight")
        plt.show()
        fig6, axes6 = plt.subplots(nrows = 1, ncols = 2)
        axes6[0].plot(n_boxes, fuel_used[5])
        axes6[0].set_xlabel("Number of boxes")
        axes6[0].set_ylabel("Fuel used [kg]")
        axes6[1].plot(n_boxes, sort_time[5])
        axes6[1].set_xlabel("Number of boxes")
        axes6[1].set_ylabel("Sortie time [sec]")
        fig5.suptitle("Fuel used and Sortie time versus cruise height")
        fig5.tight_layout()
        plt.savefig("C:\\Users\\ties\\Downloads\\Sensitivity_n_boxes")
        plt.show()
    return 

def ferryrange(ac_obj, atm_obj, aux_fuel_boxes, h0, summary = False, plot = False):
    # =============== Starting parameters =================
    W_F = ac_obj.fuelcapacity * ac_obj.fueldensity          # Max fuel
    W_F += aux_fuel_boxes * 20                              # aux_fuel_boxes represents the number of boxes swapped for fuel (stowed in boxlike objects)
    W_TO = ac_obj.W_OE + W_F                                # Ferry range means no payload
    W = W_TO                                                # Set startweight equal to take-off weight
    W_F_used = 0.0                                          # Fuel used at start of flight is 0 kg
    t = 0.0                                                 # Start time is 0 sec
    h = 0.0                                                 # Start altitude is 0 m
    x = 0.0                                                 # Start distance flown is 0 m
    # ========== Empty arrays used for plotting ===========
    x_arr = np.empty(0)
    h_arr = np.empty(0)
    # =========== Take-off ===========
    W *= ac_obj.W1W_TO * ac_obj.W2W1 * ac_obj.W3W2
    W_F_used += W_TO - W
    W_F -= W_F_used
    h += 15.0                                               # Take off maneuver ends at screenheight = 15.0 m
    x_arr = np.append(x_arr, x)
    h_arr = np.append(h_arr, h)
    # ============ Climb =============
    t_climb, x_climb, W_F_used_climb, h = fp.climbmaneuver(ac_obj, atm_obj, h, h0, W)
    t += t_climb
    x += x_climb
    W_F_used += W_F_used_climb
    W   -= W_F_used_climb
    W_F -= W_F_used_climb
    x_arr = np.append(x_arr, x)
    h_arr = np.append(h_arr, h)
    # =========== Cruise =============                      # For cruise, assume that the descent is started with 10 kg of fuel left
    CL  = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    CD  = fp.dragpolar(ac_obj, CL)
    dt  = 1.0
    while W_F > 10:
        p, rho = fp.atm_parameters(atm_obj, h)[0], fp.atm_parameters(atm_obj, h)[2]
        V   = np.sqrt(2*W*atm_obj.g/(rho*ac_obj.Sw*CL))
        Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD
        Pa  = ac_obj.power * 745.699872 * p/atm_obj.p0 * ac_obj.prop_eff
        if Pr > Pa:
            print("The combination of cruising altitude and cruise speed is not attainable")
            break
        Pbr = Pr/ac_obj.prop_eff
        FF  = Pbr * ac_obj.SFC
        x   += V * dt
        W_F_used += FF * dt
        W_F -= FF * dt
        W   -= FF * dt
        t   += dt
        x_arr = np.append(x_arr, x)
        h_arr = np.append(h_arr, h)
    # Descent
    t_des, x_des, W_F_used_des, h = fp.descentmaneuver(ac_obj, atm_obj, h, 15.0, W)
    t += t_des
    x += x_des
    W_F_used += W_F_used_des
    W_F -= W_F_used_des
    W   -= W_F_used_des
    x_arr = np.append(x_arr, x)
    h_arr = np.append(h_arr, h)
    # Landing, taxi, shutdown
    W_beforelanding = W
    W *= ac_obj.WfinalW10
    W_F_used += W_beforelanding - W
    W_F -= W_beforelanding - W
    W   -= W_beforelanding - W
    if summary:
        print("================ Ferry Summary ================")
        print(f"Ferry flight time: {round(t/3600, 2)} [hrs]")
        print(f"Ferry range: {round(x/1000, 2)} [km]")
        print(f"Fuel used: {round(W_F_used, 2)} [kg]")
        print(f"Final fuel weight: {round(W_F, 2)} [kg]")
        print("===============================================")
    if plot:
        plt.plot(x_arr/1000, h_arr)
        plt.xlabel("Horizontal distance [km]")
        plt.ylabel("Altitude [m]")
        plt.show()
    return

def VV_sortiesim_boundary(ac_obj, atm_obj):
    # Boundary value tests for loading
    sortie_nopayload = fp.fuelusesortie(ac_obj, atm_obj, 0, 0, 3048.0, 50, 54.012, Range = 250)
    sortie_12drops   = fp.fuelusesortie(ac_obj, atm_obj, 12, 12, 3048.0, 55, 54.012, 250)
    sortie_smalldz   = fp.fuelusesortie(ac_obj, atm_obj, 12, 6, 3048.0, 50, 54.012, 250, dropregion=10)
    sortie_largedz   = fp.fuelusesortie(ac_obj, atm_obj, 12, 6, 3048.0, 50, 54.012, 250, dropregion=200)
    print("=====================================================")
    print(f"Sortie 0  boxes 0  drops Range 250 [km] Dropregion None   | W_F_used: {round(sortie_nopayload[0], 2)} [kg], sortie time: {round(sortie_nopayload[2]/3600, 2)} [hr]")
    print(f"Sortie 12 boxes 12 drops Range 250 [km] Dropregion None   | W_F_used: {round(sortie_12drops[0], 2)} [kg], sortie time: {round(sortie_12drops[2]/3600, 2)} [hr]")
    print(f"Sortie 12 boxes 6  drops Range 250 [km] Dropregion 10 km  | W_F_used: {round(sortie_smalldz[0], 2)} [kg], sortie time: {round(sortie_smalldz[2]/3600, 2)} [hr]")
    print(f"Sortie 12 boxes 6  drops Range 250 [km] Dropregion 200 km | W_F_used: {round(sortie_largedz[0], 2)} [kg], sortie time: {round(sortie_largedz[2]/3600, 2)} [hr]")
    print("======================================================")

def Cessnacomparison(ac_obj, atm_obj, h_cruise, V_cruise, result = False):
    h_cruise *= 0.3048
    V_cruise *= 0.5144
    # This function assesses the model output compared to reality for a Cessna 172
    # Max fuel, payload accordingly (no drops as Cessna)
    W_F_TO = ac_obj.fuelcapacity * ac_obj.fueldensity
    W_TO   = ac_obj.W_TO
    W_F_used = 0.0
    t = 0.0
    x = 0.0
    h = 15.0
    # After take-off:
    W_a_TO = W_TO * ac_obj.W1W_TO * ac_obj.W2W1 * ac_obj.W3W2
    W = W_a_TO
    W_F_used += W_TO - W_a_TO
    # Climb
    t_climb, x_climb, W_F_used_climb, h = fp.climbmaneuver(ac_obj, atm_obj, 15.0, h_cruise, W)
    t += t_climb
    x += x_climb
    W_F_used += W_F_used_climb
    W -= W_F_used_climb
    # Cruise
    d_remain = ac_obj.R - x_climb
    CL = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    CD = fp.dragpolar(ac_obj, CL)
    d_descent = CL/CD * h_cruise
    d_cruise = d_remain - d_descent
    t_cruise, W_F_used_cruise, cruiseNAT = fp.cruisecalc(ac_obj, atm_obj, h_cruise, d_cruise, W, V_cruise)
    if cruiseNAT:
        print("This cruise speed - altitude combination is impossible")
    t += t_cruise
    x += d_cruise
    W_F_used += W_F_used_cruise
    # Descent
    t_descent, x_descent, W_F_used_descent, h = fp.descentmaneuver(ac_obj, atm_obj, h_cruise, 15.0, W)
    t += t_descent
    x += x_descent
    W_F_used += W_F_used_descent
    if result:
        print("============================================================")
        print(f"Distance flown: {round(x/1000, 2)} [km]")
        print(f"Fuel used: {round(W_F_used, 2)} [kg]")
        print(f"Fuel remaining: {round(W_F_TO - W_F_used, 2)} [kg]")
        print(f"Sortie time: {round(t/3600, 2)} [hr]")
        print("============================================================")
    return 

Cessnacomparison(Cessna, atm, 10000, 115, result=True)