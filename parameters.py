import numpy as np

class UAV:
    def __init__(self):
        "Physical parameters"
        self.g = 9.80665             # m/s^2

        "Parameters of aircraft"
        # Mission flight plan
        self.n_drops = 2
        self.R = 500000              # m
        self.W_PL = 240
        self.M_res = 0.15
        self.h_cruise = 10000*0.3048     # cruising at 10000 ft

        # Wing characteristics
        self.S = 180 # ft^2
        self.n_ult = 1.5
        self.A = 10
        self.e = 0.9
        self.t_r = 2                 #ft
        self.b = 8                   #ft
        self.sweep_angle = 0.36
        self.CD0 = 0.027
        self.CL_CD = 10              #conservative

        self.prop_eff = 0.82
        self.c_p = 90E-9

        self.W1W_TO = 0.995          # engine start-up
        self.W2W1 = 0.997            # taxi
        self.W3W2 = 0.998            # take_off
        self.W4W3 = 0.992            # climb
        self.W10W9 = 0.993           # descent
        self.WfinalW10 = 0.993       # landing, taxi, shut-down

        # # MTOW vs OEW reduced by pilot weight R2=0.9548
        #lin_par1 = 0.5522 #gradient of the linear trend OEW/MTOW
        #lin_par2 = -40.838 #y axis crossing of the linear trend OEW/MTOW

        # # MTOW for drones, R2=0.9988
        #lin_par1 = 0.4631
        #lin_par2 = 52.058

        # # MTOW vs OEW for general aviation R2=0.9548 , original y=0.5482 x + 486.68
        #lin_par1 = 0.5522
        #lin_par2 = 39.162

        # # MTOW vs OEW ultra-light, R2=0.8704
        # lin_par1 = 0.7134
        # lin_par2 = -52.981

        # MTOW vs OEW ultra-light reduced by pilot weight, R2=0.9704
        self.lin_par1 = 0.7134
        self.lin_par2 = -132.98

        self.b = 14 # m
        self.lambda_mid = 0.36 #rad
        self.n_ult = 1.5 # -
        self.t_c = 0.12 # -
        self.cwr = 2.5 # m
        self.W_loading = 562/9.81 # kg/m^2
        self.W_G = 700 # kg
        self.s_tail = 2 # m^2
        self.W_TO = 700 # kg
        self.P_TO = 62 # hp
        self.V_D = 150*0.514444 # m/s
        self.l_t = 3.5 # m
        self.b_f = 1.1 # m
        self.h_f = 1.1 # m
        self.S_G = 4*np.pi*self.b_f # m^2 #TO CHANGE
        self.N_e = 1 # -
        self.W_e = 62.6 # kg
        self.l_f = 4 # m

        self.W_fuel_estimated = 75*0.82    # estimated value based on 20L/h fuel burn of rotax and 3 hour sortie and 15L reserve and density of fuel
        self.W_LDG   = self.W_TO - self.W_PL - self.W_fuel_estimated
        self.f = self.W_LDG/self.W_TO
        self.cruise_frac = self.W1W_TO*self.W2W1*self.W3W2*self.W4W3*0.85    #assume halfway through the cruise with cruise fuel fraction 0.3

        # Parameters for geometry estimation
        # This file needs to be merged with parameters from Torben, which was not in Git by the time I started making this
        # From ADSEE 1 slides, typical values
        # Taking average of reported range of values  
        self.CL_max_clean    = 1.3   #1.3 - 1.9
        self.CL_max_TO       = 1.3
        self.CL_max_land     = 1.6
        self.CL_TO           = self.CL_max_TO/(1.1**2)
        self.CL_LDG          = self.CL_max_land/(1.1**2)
        #parameter difference in different configurations
        #ASSUMPTION: fixed undercarriage
        self.d_CD0_TO_Flaps  = np.arange(0.010, 0.021, 0.001)    #0.010 - 0.020
        self.d_CDO_LDG_Flaps = np.arange(0.055, 0.076, 0.001)    #0.055 - 0.075
        self.d_e_TO_Flaps    = 0.05
        self.d_e_LDG_Flaps   = 0.10

        #atmospheric properties
        self.p0      = 101325        #Pa
        self.rho0    = 1.225         #kg/m^3
        self.T0      = 288.15        #K
        self.Lambda  = -0.0065       #deg K/m
        #self.R       = 287.05        #J/kgK
        self.g0      = 9.80665       #m/s^2

        #speeds
        self.V_s_max = 61*(1.852/3.6)    #CS23 Vs at take off not allowed to be above 61 kts, *1.852 to get to m/s
        self.V_s_min = 50*(1.852/3.6)    #Dropping speed
        self.V_cruise= 105*(1.852/3.6)   #Cruise speed
        self.V_TO_max= 1.1*self.V_s_max       #Maximum take off speed
        self.V_TO_min= 1.1*self.V_s_min       #Minimum take off speed
        self.V_climb = 70*(1.852/3.6)

        #power
        self.P_max   = 100           #break horse power
        self.eta_p   = 0.82          #propeller efficiency
        self.power_setting = 0.9     #cruise power

        #Take off distance
        self.TO_dist = 750           #m, seems long, might need to change this, STOL gang
        self.LDG_dist= 750
        self.lin_par1 = 0.5482       #gradient of the linear trend OEW/MTOW
        self.lin_par2 = 486.68       #y axis crossing of the linear trend OEW/MTOW

    # def Class_I():
        
