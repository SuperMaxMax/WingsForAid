import numpy as np

class UAV:
    def __init__(self, name, engine_pos, boom):
        "==== Aircraft Parameters ===="
        self.name                = name          # Name of the aircraft [-]

        "-CS23 Type"
        self.type                = "normal"      # CS23 aircraft type: "normal" for normal/commuter and "utility" for utility          

        "-Aircraft geometry"
        self.Sw                  = 13            # Surface area [m2] # Sw,b and MGC to be removed -> then V_n diagrams should be created after iteration instead of separate
        self.A                   = 8             # Aspect ratio [-]
        self.e                   = 0.7           # Oswald factor [-]
        self.b                   = 11            # Wing span [m]
        self.MGC                 = self.Sw / self.b #Mean geometric chord [m]
        self.braced_wing         = False         # True if wing is braced

        self.s_tail              = 2             # Tail surface area [m]
        self.l_t                 = 3.5           # Tail arm [m]

        self.boom                = boom          # Boom, true if boom tail is implemented
        self.W_boom              = 0            # Boom weight [kg]
        self.l_f_boom            = 2             # Boom length [m]

        self.xc_OEW_p            = 0.2          # Center of gravity of OEW as a fraction of the fuselage length [-]

        self.pos_main_carriage   = "fuselage"    # Position of main carriage: "fuselage" or "wing"
        self.main_gear_type      = "fixed"       # Type of main gear: "fixed" or "retractable"
        self.nose_gear_type      = "fixed"       # Type of nose gear: "fixed" or "retractable"

        "-Aerodynamic properties"
        self.CD0                 = 0.027         # Zero lift coefficient [-]
        self.CLa                 = 4.2          # Lift curve slope [] | CHANGE TO ACTUAL VALUE

        # ASSUMPTION/NOTES: ADSEE 1 slides mention ranges for CL, the code automatically runs over all the CL's in these lists
        # but this means that CL_max_clean, CL_max_TO and CL_max_land must always be stored in an array. For an array with length 1
        # the code just runs once
        self.CL_max_clean        = np.array([1.3])              # Maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.CL_max_TO           = np.array([1.3])              # Maximum lift coefficient at take-off [-]
        self.CL_max_land         = np.array([1.9])              # Maximum lift coefficient at landing [-]
        self.CL_TO               = self.CL_max_TO / (1.1**2)    # [-]
        self.CL_LDG              = self.CL_max_land / (1.1**2)  # [-]

        "-Weights"
        self.W_e                 = 62.6          # Definitive weight per engine [kg]
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
        self.cruise_frac         = self.W1W_TO*self.W2W1*self.W3W2*self.W4W3*0.85   # Assume halfway through the cruise with cruise fuel fraction 0.3

        "-Propulsive properties"
        self.engine_pos          = engine_pos     # Engine position
        self.P_TO                = 62            # Power at take-off [hp]

        self.prop_eff            = 0.82          # Propulsive efficiency [-]
        self.eta_p               = self.prop_eff # Propulsive efficiency [-]
        if engine_pos == "pusher":
            self.prop_eff        *= 0.92         # Propulsive efficiency [-]
            self.eta_p           *= 0.92         # Propulsive efficiency [-]

        self.power_setting       = 0.9           # Power setting in cruise [-]

        self.c_p                 = 72E-9         #
        self.N_e                 = 1             # Number of engines [-]

        "==== Mission profile/Atmospheric properties ===="
        "-Mission characteristics"
        self.n_drops             = 1             # Number of drops [-]
        self.n_boxes             = 12            # Number of boxes [-]
        self.R                   = 500000        # Range [m]
        self.M_res               = 0.10          # Fraction of remaing fuel at the end of flight/reserve fuel [-]
        self.h_cruise            = 10000*0.3048  # Cruise altitude [m] | NOTES: Conversion
        self.h_TO                = 0             # Take-off Height [m]
        
        self.LDG_dist            = 750           # Landing distance [m]

        self.n_ult               = 3.8 * 1.5     # Ultimate load factor [-]

        "-Speeds"
        self.V_s_min             = 50*(1.852/3.6)    # Dropping speed [m/s]
        self.V_cruise            = 105*(1.852/3.6)   # Cruise speed [m/s]
        self.V_climb             = 70*(1.852/3.6)    # Climb speed [m/s]
        self.V_D                 = 140*0.514444      # Dive speed [m/s]
        self.V_B                 = 46.01347201449718 # Design speed for maximum gust intensity [m/s] | NOTES: Follow guidelines to choose this speed

        "-Atmospheric properties"
        self.rho0                = 1.225         # [kg/m^3]
        self.T0                  = 288.15        # [K]
        self.Lambda              = -0.0065       # [deg K/m]
        self.R_gas               = 287.05        # [J/kgK]
        self.g0                  = 9.80665       # [m/s^2]

        "==== Miscellaneous ===="
        "MTOW vs OEW GA, "
        self.lin_par1            = 0.5249
        self.lin_par2            = 42.049