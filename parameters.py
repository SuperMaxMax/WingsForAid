import numpy as np
class UAV:
    def __init__(self, name, engine_pos, boom, braced_wing):
        "==== Aircraft Parameters ===="
        self.name                = name         # Name of the aircraft [-]

        "-CS23 Type"
        self.type                = "utility"    # CS23 aircraft type: "normal" for normal/commuter and "utility" for utility          

        "-Aircraft geometry"
        self.A                   = 10            # Aspect ratio [-]
        self.e                   = 0.7           # Oswald factor [-]
        self.braced_wing         = braced_wing   # True if wing is braced
        self.kq                  = 0.95          # Volume factor used to calculate wetted area of the wing

        self.s_tail              = 2             # Tail surface area [m]
        self.l_t                 = 3.5           # Tail arm [m]

        self.boom                = boom          # Boom, true if boom tail is implemented
        self.W_boom              = 20            # Boom weight [kg]
        self.l_f_boom            = 2             # Boom length [m]

        self.xc_OEW_p            = 0.25         # Center of gravity of OEW as a fraction of the MAC [-]

        self.pos_main_carriage   = "fuselage"    # Position of main carriage: "fuselage" or "wing"
        self.main_gear_type      = "fixed"       # Type of main gear: "fixed" or "retractable"
        self.nose_gear_type      = "fixed"       # Type of nose gear: "fixed" or "retractable"

        self.mass_penalty_struts = 7             # The weight of the struts

        "- Fuselage geometry"
        self.side_clearance      = 0.2           # Side clearance [m], this is for both sides
        self.top_clearance       = 0.2   
        self.bot_clearance       = 0.1
        self.structural_thickness= 0.2           # Structural thickness fuselage [m], this is for both sides
        
        "-Aerodynamic properties"
        self.CD0                 = 0.027         # Zero lift coefficient [-]
        self.CLa                 = 4.2           # Lift curve slope [] | CHANGE TO ACTUAL VALUE
        self.Drag_increase       = 1.0           # This is used for the calculations of the strut drag if applicable

        # ASSUMPTION/NOTES: ADSEE 1 slides mention ranges for CL, the code automatically runs over all the CL's in these lists
        # but this means that CL_max_clean, CL_max_TO and CL_max_land must always be stored in an array. For an array with length 1
        # the code just runs once
        self.CL_max_clean        = np.array([1.5])              # Maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.CL_max_TO           = np.array([1.5])              # Maximum lift coefficient at take-off [-]
        self.CL_max_land         = np.array([1.9])              # Maximum lift coefficient at landing [-]
        self.CL_TO               = self.CL_max_TO / (1.1**2)    # [-]
        self.CL_LDG              = self.CL_max_land / (1.1**2)  # [-]

        "-Weights"
        self.W_e                 = 62.6         # Definitive weight per engine [kg]
        self.W_TO                = 700          # Take-off weight [kg]
        self.W_PL                = 240          # Payload weight [kg]
        self.WS                  = 708.05       # Wing Loading [N/m^2]

        "-Weight fractions"
        self.W1W_TO              = 0.995        # Engine startup fraction [-]
        self.W2W1                = 0.997        # Taxi fraction [-]
        self.W3W2                = 0.998        # Take_off fraction [-]
        self.W4W3                = 0.992        # Climb fraction [-]
        self.W10W9               = 0.993        # Descent fraction [-]
        self.WfinalW10           = 0.993        # Landing, taxi & shut-down fraction [-]
        self.cruise_frac         = self.W1W_TO*self.W2W1*self.W3W2*self.W4W3*0.85   # Assume halfway through the cruise with cruise fuel fraction 0.3 [-]

        "-Propulsive properties"
        self.engine_pos          = engine_pos   # Engine position [-]
        self.engine_length       = 0.6651       # Engine length [m] - EASA type certificate data sheet ROTAX 912 series
        self.engine_cg           = 0.327        # Engine cg [m] - for ROTAX 912is engine
        self.engine_fairing      = 0.2          # Engine fairing length [m]
        self.d_engine_boxes      = 0.4          # Distance between engine boxes [m], leaves room for possible firewall
        self.power               = 80           # Power at take-off [hp]

        self.prop_eff            = 0.82         # Propulsive efficiency [-]
        self.eta_p               = self.prop_eff# Propulsive efficiency [-]
        if engine_pos == "pusher":
            self.prop_eff        *= 0.92        # Propulsive efficiency [-]
            self.eta_p           *= 0.92        # Propulsive efficiency [-]

        self.power_setting       = 0.9          # Power setting in cruise [-]

        self.c_p                 = 72E-9        # Specific fuel consumption [kg/J]
        self.N_e                 = 1            # Number of engines [-]

        "==== Mission profile/Atmospheric properties ===="
        "-Mission characteristics"
        self.n_drops             = 1             # Number of drops [-]
        self.n_boxes             = 12            # Number of boxes [-]
        self.n_boxes_abreast     = 2             # Number or boxes side by side [-]
        self.n_rows              = 6
        self.R                   = 500000        # Range [m]
        self.M_res               = 0.10          # Fraction of remaing fuel at the end of flight/reserve fuel [-]
        self.h_cruise            = 10000*0.3048  # Cruise altitude [m] | NOTES: Conversion
        self.h_TO                = 0             # Take-off Height [m]
        
        self.LDG_dist            = 750          # Landing distance [m]

        self.n_ult               = 4.4 * 1.5    # Ultimate load factor [-]

        "-Speeds"
        self.V_s_min             = 50*(1.852/3.6)       # Dropping speed [m/s]
        self.V_cruise            = 116.82*(1.852/3.6)   # Cruise speed [m/s]
        self.V_climb             = 70*(1.852/3.6)       # Climb speed [m/s]
        self.V_D                 = 163.54*0.514444      # Dive speed [m/s]
        self.V_B                 = 42.22                # Design speed for maximum gust intensity [m/s] | NOTES: Follow guidelines to choose this speed

        "-Atmospheric properties"
        self.rho0                = 1.225        # [kg/m^3]
        self.T0                  = 288.15       # [K]
        self.Lambda              = -0.0065      # [deg K/m]
        self.R_gas               = 287.05       # [J/kgK]
        self.g0                  = 9.80665      # [m/s^2]

        "==== Miscellaneous ===="
        "MTOW vs OEW GA, "
        self.lin_par1            = 0.5249       # [-]
        self.lin_par2            = 42.049       # [-]

class Cessna_172:
    def __init__(self, engine_pos, braced_wing, boom):
        "==== Aircraft Parameters ===="
        self.name                = "Cessna 172"         # Name of the aircraft [-]

        "-CS23 Type"
        self.type                = "normal"    # CS23 aircraft type: "normal" for normal/commuter and "utility" for utility          

        "-Aircraft geometry"
        self.A                   = 7.32         # Aspect ratio [-]
        self.e                   = 0.8          # Oswald factor [-]
        self.braced_wing         = braced_wing  # True if wing is braced [-]
        self.kq                  = 0.95         # Volume factor used to calculate wetted area of the wing [-]
        self.boom                = boom
        self.s_tail              = 2.83         # Tail surface area [m^2]
        self.l_t                 = 4.14         # Tail arm [m]

        self.xc_OEW_p            = 0.25         # Center of gravity of OEW as a fraction of the MAC [-]

        self.pos_main_carriage   = "fuselage"    # Position of main carriage: "fuselage" or "wing" [-]
        self.main_gear_type      = "fixed"      # Type of main gear: "fixed" or "retractable" [-]
        self.nose_gear_type      = "fixed"      # Type of nose gear: "fixed" or "retractable" [-]

        "-Aerodynamic properties"
        self.CD0                 = 0.027        # Zero lift coefficient [-]
        self.CLa                 = 5.143          # Lift curve slope [-] | change to actual value (follows from aero analysis later)
        self.Drag_increase       = 1.0          # This is used for the calculations of the strut drag if applicable [-]

        self.CL_max_clean        = np.array([1.6])              # Maximum lift coefficient [-], range: 1.3 - 1.9 | CL_max_clean, CL_max_TO and CL_max_land must always be stored in an array
        self.CL_max_TO           = np.array([1.6])              # Maximum lift coefficient at take-off [-]       | Multiple values can be in the arrays
        self.CL_max_land         = np.array([2.3])              # Maximum lift coefficient at landing [-]
        self.CL_TO               = self.CL_max_TO / (1.1**2)    # [-]
        self.CL_LDG              = self.CL_max_land / (1.1**2)  # [-]

        "-Weights"
        self.W_e                 = 124.74       # Definitive weight per engine [kg]
        self.W_TO                = 1088.62      # Take-off weight [kg]
        self.W_PL                = 415.99       # Payload weight [kg]
        self.WS                  = 658.99       # Wing Loading [N/m^2]

        "-Weight fractions"
        self.W1W_TO              = 0.995        # Engine startup fraction [-]
        self.W2W1                = 0.997        # Taxi fraction [-]
        self.W3W2                = 0.998        # Take_off fraction [-]
        self.W4W3                = 0.992        # Climb fraction [-]
        self.W10W9               = 0.993        # Descent fraction [-]
        self.WfinalW10           = 0.993        # Landing, taxi & shut-down fraction [-]
        self.cruise_frac         = self.W1W_TO*self.W2W1*self.W3W2*self.W4W3*0.85   # Assume halfway through the cruise with cruise fuel fraction 0.3 [-]

        "-Propulsive properties"
        self.engine_pos          = engine_pos   # Engine position [-]
        self.power               = 160          # hp

        self.prop_eff            = 0.8         # Propulsive efficiency [-]
        self.eta_p               = self.prop_eff# Propulsive efficiency [-]
        if engine_pos == "pusher":
            self.prop_eff        *= 0.92        # Propulsive efficiency [-]
            self.eta_p           *= 0.92        # Propulsive efficiency [-]

        self.power_setting       = 0.7          # Power setting in cruise [-]

        self.c_p                 = 7.266E-8     # Specific fuel consumption [kg/J]
        self.N_e                 = 1            # Number of engines [-]

        "==== Mission profile/Atmospheric properties ===="
        "-Mission characteristics"
        self.n_drops             = 0            # Number of drops [-]
        self.n_boxes             = 0            # Number of boxes [-]
        self.R                   = 1185000      # Range [m]
        self.M_res               = 0.07         # Fraction of remaining fuel at the end of flight/reserve fuel [-]
        self.h_cruise            = 10000*0.3048 # Cruise altitude [m] | NOTES: Conversion
        self.h_TO                = 0            # Take-off Height [m]
        
        self.LDG_dist            = 420          # Landing distance [m]

        self.n_ult               = 3.8 * 1.5    # Ultimate load factor [-]

        "-Speeds"
        self.V_s_min             = 50*(1.852/3.6)       # Dropping speed [m/s]
        self.V_cruise            = 122*(1.852/3.6)      # Cruise speed [m/s]
        self.V_climb             = 76*(1.852/3.6)       # Climb speed [m/s]
        self.V_D                 = 175*0.514444         # Dive speed [m/s]
        self.V_B                 = 127*0.514444         # Design speed for maximum gust intensity [m/s] | NOTES: Follow guidelines to choose this speed

        "-Atmospheric properties"
        self.rho0                = 1.225        # [kg/m^3]
        self.T0                  = 288.15       # [K]
        self.Lambda              = -0.0065      # [deg K/m]
        self.R_gas               = 287.05       # [J/kgK]
        self.g0                  = 9.80665      # [m/s^2]

        "==== Miscellaneous ===="
        "MTOW vs OEW GA, "
        self.lin_par1            = 0.5522       # [-]
        self.lin_par2            = 39.162       # [-]