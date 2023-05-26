import numpy as np
class UAV:
    def __init__(self, name):
        "=== Class I / Class II parameters ==="
        self.A = 10                         # Aspect ratio [-]
        self.BHP_cruise = 76.3436
        self.CD0 = 0.0273                   # Zero lift coefficient [-]
        self.CL_LDG = 1.5702                # [-]
        self.CL_TO = 1.2397                 # [-]
        self.CL_max_TO = 1.5                # Maximum lift coefficient at take-off [-]
        self.CL_max_clean = 1.5             # Maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.CL_max_land = 1.9              # Maximum lift coefficient at landing [-]
        self.CLa = 4.2                      # Lift curve slope [-] | CHANGE TO ACTUAL VALUE
        self.Drag_increase = 1.0126         # This is used for the calculations of the strut drag if applicable [-]
        self.RFL = 750.0                    # Required field length [m]
        self.LDG_dist = 750.0               # Landing distance [m]
        self.L_D = 14.1804                  # Lift to drag ratio [-]
        self.Lambda = -0.0065
        self.MAC_length = 1.0161            # Mean aerodynamic chord [m]
        self.MGC = 1.0822                   # Mean geometric chord [m]
        self.M_res = 0.1 
        self.Mff = 0.9254                   # Fuel fraction [-]
        self.N_e = 1                        # Number of engines [-]
        self.R = 500000                     # Range [m]
        self.R_gas = 287.05                 # Gas constant [J/kgK]
        self.S_G = 19.77                    # Gross shell area fuselage [m^2]
        self.Sw = 11.7113                   # Wing area [m^2]
        self.Sw_wetted = 23.4226            # Wetted area of the wing [m^2]
        self.T0 = 288.15                    # Sea level temperature [K]
        self.TOP_req = 250 
        self.V_A = 47.9382                  # [m/s]
        self.V_B = 45.0                     # [m/s]
        self.V_D = 90.7328                  # [m/s]
        self.V_climb = 36.0111              # [m/s]
        self.V_cruise = 60.4885             # Cruise velocity [m/s]
        self.V_s_min = 22.8537              # Minimum stall velocity [m/s]
        self.W10W9 = 0.993                  # Descent fraction [-]
        self.W1W_TO = 0.995                 # Engine startup fraction [-]
        self.W2W1 = 0.997                   # Taxi fraction [-]
        self.W3W2 = 0.998                   # Take-off fraction [-]
        self.W4W3 = 0.992                   # Climb fraction [-]
        self.WP = 0.1215
        self.WS = 607.8751                  # Wing Loading [N/m^2]
        self.W_F = 59.8091                  # Fuel weight [kg]
        self.W_OE = 429.1354                # Operational empty weight [kg]
        self.W_PL = 240                     # Payload weight [kg]
        self.W_TO = 728.9445                # Take-off weight [kg]
        self.W_boom = 20                    # Boom weight [kg]
        self.W_e = 62.6                     # Definitive weight per engine [kg]
        self.W_eq = 58.075                  # Equipment weight [kg]
        self.W_fus = 95.9102                # Fuselage weight [kg]
        self.W_n = 11.1013                  # Nacelle weight [kg]
        self.W_pg = 84.7333                 # Propulsion group weight [kg]
        self.W_sc = 27.2948                 # Control systems weight [kg]
        self.W_t = 17.0451                  # Tail weight [kg]
        self.W_uc = 59.3859                 # Undercarriage weight [kg]
        self.W_w = 55.5897                  # Wing weight [kg]
        self.WfinalW10 = 0.993              # Landing, taxi & shut-down fraction [-]
        self.X_LEMAC = 2.276                # Leading edge mean aerodynamic chord [m]
        self.X_cg_aft = 0.5335              # Aft cg location CG/MAC [-]
        self.X_cg_full = 0.4115             # MTOW cg location CG/MAC [-]
        self.X_cg_fwd = 0.1704              # Forward cg location CG/MAC [-]
        self.X_cg_range = 0.363             # Range of cg location CG/MAC [-]
        self.b = 10.8219                    # Wing span [m]
        self.boom                = True     # Boom, true if boom tail is implemented
        self.bot_clearance = 0.1            # Bottom clearance [m]
        self.braced_wing         = True     # True if wing is braced
        self.boxweight = 23.0               # weight of a single box [kg]
        self.c_p = 72E-9                    # Specific fuel consumption [kg/J]
        self.climb_rate = 2.9889 
        self.cos_lambda_c04 = 1 
        self.cruise_frac = 0.8348           # Assume halfway through the cruise with cruise fuel fraction 0.3 [-]
        self.d_eff = 1.241                  # Effective diameter [m]
        self.d_engine_boxes = 0.4           # Distance between engine and wing box [m]
        self.dihedral = 1
        self.e                   = 0.7      # Oswald factor [-]
        self.engine_cg = 0.327              # Engine cg location [m]
        self.engine_fairing = 0.2           # Engine fairing length [m]
        self.engine_length = 0.6651         # Engine length [m]
        self.engine_pos = 'tractor'         # Engine position: "tractor" or "pusher" or "fuselage"
        self.eta_p = 0.82                   # Propulsive efficiency [-]
        self.g0 = 9.8066                    # Gravitational acceleration [m/s^2]
        self.h_TO = 0                       # Take-off altitude, airport altitude [m]
        self.h_cruise = 3048.0              # Cruise altitude [m]
        self.h_in = 0.9                     # Inner fuselage height [m]
        self.h_out = 1.1                    # Outer fuselage height [m]
        self.kq = 0.95                      # Volume factor used to calculate wetted area of the wing [-]
        self.l_f = 5.4651                   # Fuselage length [m]
        self.l_f_boom = 2                   # Boom length [m]
        self.l_n = 0.8651                   # Nosecone length [m]
        self.l_t = 3.5                      # Tail arm [m]
        self.l_tc = 0.8                     # Tail cone length [m]
        self.lambda_co2 = -0.0428           # Half chord sweep angle [rad]
        self.lambda_co4 = 0.0               # Quarter chord sweep angle [rad]
        self.lin_par1 = 0.5249              # [-]
        self.lin_par2 = 42.049              # [-]
        self.main_gear_type = 'fixed'       # Type of main gear: "fixed" or "retractable"
        self.mass_penalty_struts = 7        # The weight of the struts [kg]
        self.n_boxes = 12                   # [-]
        self.n_boxes_abreast = 2            # [-]
        self.n_drops = 1                    # [-]
        self.n_rows = 6                     # [-]
        self.n_ult = 6.6                    # Ultimate load factor [-]
        self.name = name                    # Name of the aircraft
        self.nose_gear_type = 'fixed'       # Type of nose gear: "fixed" or "retractable"
        self.pos_main_carriage = 'fuselage' # Position of main carriage: "fuselage" or "wing"
        self.power = 95.8347                # Power at takeoff [hp]
        self.power_setting = 0.9            # Power setting in cruise [-]
        self.prop_eff = 0.82                # Propulsive efficiency [-]
        self.rho0 = 1.225                   # Air density at sea level [kg/m^3]
        self.rho_TO = 1.225                 # Take-off air density if airport is at sea level [kg/m^3]
        self.rho_cruise = 0.9046            # Cruise air density [kg/m^3]
        self.rootchord = 1.546              # Root chord [m]
        self.s_tail = 2                     # Tail surface area [m]
        self.side_clearance = 0.2           # Side clearance [m], this is for both sides
        self.sigma_TO = 1.0
        self.sigma_cruise = 0.7385
        self.structural_thickness = 0.2     # Structural thickness fuselage [m], this is for both sides
        self.t_c = 0.12                     # Thickness to chord ratio [-]
        self.taper = 0.4                    # Taper ratio [-]
        self.tipchord = 0.6184              # Tip chord [m]
        self.top_clearance = 0.2            # Top clearance [m]
        self.type = "utility"               # CS23 aircraft type: "normal" for normal/commuter and "utility" for utility    
        self.w_in = 1.2                     # Inner fuselage width [m]
        self.w_out = 1.4                    # Outer fuselage width [m]
        self.x_lemac = 0.2871               # Distance from LE root chord to the leading edge mean aerodynamic chord [m]
        self.xc_OEW_p = 0.25                # Center of gravity of OEW as a fraction of the MAC [-]
        self.y_mac = 2.3182                 # Spanwise location of the MAC [m]

        "Structural parameters"             # NOTE: Add identifier "ST_" before variable names
        self.something = 1 # add units
        self.y_mac = 2                      # Spanwise location of the MAC [m]
        self.ST_SF = 1.5

        "Aerodynamic parameters"            # NOTE: Add identifier "AE_" before variable names
        self.AE_Cl0 = 0.4                   # TODO: Change to real value - Lift coeff of airfoil @ 0 AOA, cruise velocity [-]

        "Flight Performance parameters"     # NOTE: Add identifier "FP_" before variable names
        self.screenheight = 50*0.3048       # screen height of 50 ft (CS23)
        self.rpm_maxcont  = 5500            # rpm
        self.omega_prop   = 237             # rad/s, based on 5500 rpm max continuous power and 2.43 gearbox ratio
        self.prop_radius  = 0.8255          # [m] based on 3 blade rotax 3B0 ground adjustable propeller by sensenich propellers
        self.ceiling      = 18000*0.3048    # [m] 18000 ft service ceiling
        self.th_ceil      = 30000*0.3048
        self.SFC          = 7.91666667E-8   # kg/J specific fuel consumption

        "Control and stability parameters"  # NOTE: Add identifier "CS_" before variable names
        self.CS_eta = 0.95                  # airfoil efficiency factor [-]
        # self.CS_cf = 0.3                    # flap chord [m] TODO: update value
        self.CS_mu1 = 0.24
        self.CS_mu2 = 0.78
        self.CS_mu3 = 0.525
        self.CS_x_ac_w = 0.25               # location of wing ac, divided by MAC [-] - SEAD L7, S34   
        self.CS_l_h = 6.5                   # [m] tail length; length of aerodynamic centre of wing to aerodynamic centre tail. NOTE: This is a design choice, so for now it is a guestimate.
        self.CS_Cm_0_airfoil = -0.083       # TODO: Update value - Moment coefficient of airfoil [-]
        
        self.Vh_V = 0.95                    # Ratio between velocity at tail and wing [-] NOTE: This is a guestimate
        self.A_h = 6                        # Aspect ratio horizontal tail. NOTE: This is a guestimate  
        self.lambda_co2_h = 0               # [rad] Half chord sweep of horizontal tailplane [-] NOTE: This is a guestimate  
        self.dEpsilondA = 0.02              # Downwash [-] TODO: check this value, this is a pure guess
        self.Sh_S = 0.3

        "Operations parameters"             # NOTE: Add identifier "OP_" before variable names
        self.something = 1 # add units

class airport:
    def __init__(self, name):
        self.name       = name
        self.mu_ground  = 0.05                  #buildingspeed.org
        self.rwyslope   = 1.0                   #runway slope in degrees --> CONVERT TO RADIANS

class atmosphere:
    def __init__(self):
        self.rho0   = 1.225     # kg/m^3
        self.lambd  = -0.0065   # troposphere, deg K/m
        self.T0     = 273.15    # K
        self.g      = 9.80665   # m/s^2
        self.R      = 287.05    
        self.p0     = 101325    # Pa 
        self.gamma  = 1.4

