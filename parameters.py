import numpy as np

class UAV:
    def __init__(self, name, engine_pos, boom):
        "==== Aircraft Parameters ===="
        self.name                = name          # Name of the aircraft [-]

        "-CS23 Type"
        self.type                = "normal"      # CS23 aircraft type: "normal" for normal/commuter and "utility" for utility          

        "-Aircraft geometry"
        self.Sw                  = 13            # Surface area [m2]
        self.A                   = 8             # Aspect ratio [-]
        self.e                   = 0.7           # Oswald factor [-]
        self.b                   = 11            # Wing span [m]
        self.sweep_angle         = 0.36          # Sweep angle [rad]
        self.lambda_mid          = 0.36          # Sweep angle at mid-wing [rad]
        self.t_c                 = 0.12          # Thickness over chord ratio [-]
        self.rootchord           = 2.5           # Chord length at root [m]
        self.dihedral            = 1             # Wing dihedral [deg]
        self.braced_wing         = False         # True if wing is braced

        self.h_out               = 1.1           # Fuselage width [m]
        self.w_out               = 1.1           # Fuselage height [m]
        self.d_eff               = 1.421         # meter, ADSEE 1 slides 
        self.l_f                 = 4             # Fuselage length [m]
        self.S_G                 = 21.029        # Gross shell area of fuselage [m^2]

        self.s_tail              = 2             # Tail surface area [m]
        self.l_t                 = 3.5           # Tail arm [m]

        self.boom                = boom          # Boom, true if boom tail is implemented
        self.b_boom              = 0.15          # Boom width [m]
        self.h_boom              = 0.15          # Boom height [m]
        self.d_eff_boom          = np.sqrt(self.b_boom*self.h_boom) # 
        self.l_f_boom            = 2             # Boom length [m]
        self.S_G_boom            = self.l_f_boom*np.pi*self.d_eff_boom + 2*np.pi*(self.d_eff_boom/2)**2

        self.l_t_boom            = self.l_f_boom+0.4*self.l_f           # Boom tail arm [m]

        self.xc_OEW_p            = 0.25          # Center of gravity of OEW as a fraction of the fuselage length [-]

        self.pos_main_carriage   = "fuselage"    # Position of main carriage: "fuselage" or "wing"
        self.main_gear_type      = "fixed"       # Type of main gear: "fixed" or "retractable"
        self.nose_gear_type      = "fixed"       # Type of nose gear: "fixed" or "retractable"

        "-Aerodynamic properties"
        self.CD0                 = 0.027         # Zero lift coefficient [-]
        self.L_D                 = 10            # Lift over drag [-] | ASSUMPTION/NOTES: Conservative

        # ASSUMPTION/NOTES: ADSEE 1 slides mention ranges for CL, the code automatically runs over all the CL's in these lists
        # but this means that CL_max_clean, CL_max_TO and CL_max_land must always be stored in an array. For an array with length 1
        # the code just runs once
        self.CL_max_clean        = np.array([1.3])              # Maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.CL_max_TO           = np.array([1.3])              # Maximum lift coefficient at take-off [-]
        self.CL_max_land         = np.array([1.6])              # Maximum lift coefficient at landing [-]
        self.CL_TO               = self.CL_max_TO / (1.1**2)    # [-]
        self.CL_LDG              = self.CL_max_land / (1.1**2)  # [-]

        # ASSUMPTION: Fixed undercarriage
        self.d_CD0_TO_Flaps      = np.arange(0.010, 0.021, 0.001)    # Change in CD0 at take-off due to flaps [-] | Range: 0.010 - 0.020
        self.d_CDO_LDG_Flaps     = np.arange(0.055, 0.076, 0.001)    # Change in CD0 at landing due to flaps [-] | Range: 0.055 - 0.075
        self.d_e_TO_Flaps        = 0.05                              # Change in oswald factor at take-off due to flaps [-] 
        self.d_e_LDG_Flaps       = 0.10                              # Change in oswald factor at landing due to flaps [-]

        "-Weights"
        self.W_e                 = 62.6          # Definitive weight per engine [kg]
        self.W_G                 = 700           # Gross weight [kg]
        self.W_TO                = 700           # Take-off weight [kg]
        self.W_PL                = 240           # Payload weight [kg]
        self.WS                  = 600           # Wing Loading [N/m^2]

        "-Weight fractions"
        self.W1W_TO              = 0.995         # Engine startup fraction [-]
        self.W2W1                = 0.997         # Taxi fraction [-]
        self.W3W2                = 0.998         # Take_off fraction [-]
        self.W4W3                = 0.992         # Climb fraction [-]
        self.W10W9               = 0.993         # Descent fraction [-]
        self.WfinalW10           = 0.993         # Landing, taxi & shut-down fraction [-]

        # self.W_fuel_estimated    = 75*0.82                                          # Fuel weight [kg] | ASSUMPTION: Estimated value based on 20L/h fuel burn of rotax and 3 hour sortie and 15L reserve and density of fuel
        # self.W_LDG               = self.W_TO - self.W_PL - self.W_fuel_estimated    # Weight at landing [kg] | ASSUMPTION: All packages have been dropped
        # self.fuel_frac           = self.W_LDG/self.W_TO                             # Fuel fraction [-]
        self.cruise_frac         = self.W1W_TO*self.W2W1*self.W3W2*self.W4W3*0.85   # Assume halfway through the cruise with cruise fuel fraction 0.3

        "-Propulsive properties"
        self.engine_pos          = engine_pos     # Engine position
        self.P_max               = 100           # Maximum power [bhp]
        self.P_TO                = 62            # Power at take-off [hp]

        self.prop_eff            = 0.82          # Propulsive efficiency [-]
        self.eta_p               = self.prop_eff # Propulsive efficiency [-]

        self.power_setting       = 0.9           # Power setting in cruise [-]

        self.c_p                 = 72E-9         #
        self.N_e                 = 1             # Number of engines [-]

        "==== Mission profile/Atmospheric properties ===="
        "-Mission characteristics"
        self.n_drops             = 1             # Number of drops [-]
        self.n_boxes             = 12            # Number of boxes [-]
        self.R                   = 500000        # Range [m]
        self.R_ferry             = 1000000       # Ferry range [m]
        self.M_res               = 0.10          # Fraction of remaing fuel at the end of flight/reserve fuel [-]
        self.h_cruise            = 10000*0.3048  # Cruise altitude [m] | NOTES: Conversion
        self.h_TO                = 0             # Take-off Height [m]

        self.TO_dist             = 750           # Take-off distance [m]           
        self.LDG_dist            = 750           # Landing distance [m]

        self.n_ult               = 5.7           # Ultimate load factor [-]

        "-Speeds"
        self.V_s_max             = 61*(1.852/3.6)    # CS23 Vs at take off not allowed to be above 61 kts [m/s] | NOTES: *1.852 to get to m/s
        self.V_s_min             = 50*(1.852/3.6)    # Dropping speed [m/s]
        self.V_cruise            = 105*(1.852/3.6)   # Cruise speed [m/s]
        self.V_TO_max            = 1.1*self.V_s_max  # Maximum take off speed [m/s]
        self.V_TO_min            = 1.1*self.V_s_min  # Minimum take off speed [m/s]
        self.V_climb             = 70*(1.852/3.6)    # Climb speed [m/s]
        self.V_D                 = 150*0.514444      # Dive speed [m/s]

        "-Atmospheric properties"
        self.p0                  = 101325        # [Pa]
        self.rho0                = 1.225         # [kg/m^3]
        self.T0                  = 288.15        # [K]
        self.Lambda              = -0.0065       # [deg K/m]
        self.R_gas               = 287.05        # [J/kgK]
        self.g0                  = 9.80665       # [m/s^2]

        "==== Miscellaneous ===="

        #self.lin_par1            = 0.5482        # Gradient of the linear trend OEW/MTOW
        #self.lin_par2            = 486.68        # Y axis crossing of the linear trend OEW/MTOW

        "-MTOW vs OEW reduced by pilot weight R2=0.9548"
        #self.lin_par1           = 0.5522
        #self.lin_par2           = -40.838

        "-MTOW for drones, R2=0.9988"
        #self.lin_par1           = 0.4631
        #self.lin_par2           = 52.058

        "-MTOW vs OEW for general aviation R2=0.9548 , original y=0.5482 x + 486.68"
        #lin_par1           = 0.5522
        #lin_par2           = 39.162

        "- MTOW vs OEW ultra-light, R2=0.8704"
        #lin_par1           = 0.7134
        #lin_par2           = -52.981

        "- MTOW vs OEW ultra-light reduced by pilot weight, R2=0.9704"
        #self.lin_par1           = 0.7134
        #self.lin_par2           = -132.98

        "MTOW vs OEW GA, "
        self.lin_par1            = 0.5249
        self.lin_par2            = 42.049

