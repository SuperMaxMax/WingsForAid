import numpy as np
class UAV:
    def __init__(self, name):
        "=== Class I / Class II parameters ==="
        # A
        self.A = 7.75                       # Aspect ratio [-]
        self.alpha_f = 0                     # Still to be updated angle of attack of the fuselage [rad]

        # B
        self.BHP_cruise = 76.3436
        self.CD0 =  0.027                    # Zero lift coefficient [-]
        self.CL_land = 1.5702                # [-]
        self.b = 9.527          # Wing span [m]
        self.b_v = 1.60         # Vertical tail span [m] NOTE: Is this value updated during design vertical tail?
        self.boom                = True     # Boom, true if boom tail is implemented
        self.bot_clearance = 0.1            # Bottom clearance [m]
        self.braced_wing         = True     # True if wing is braced
        self.boxweight = 23.0               # weight of a single box [kg]
        # C
        self.CD0 =  0.02578                 # Zero lift coefficient [-]
        self.CL_LDG = 1.5702                # [-]
        self.CL_TO = 1.2397                 # [-]
        self.CL_max_TO = 1.5                # Maximum lift coefficient at take-off [-]
        self.CL_max_clean = 1.5             # Maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.CL_max_land = 1.9                # Maximum lift coefficient at landing [-]
        self.CLa = 4.2                      # Lift curve slope [-] | CHANGE TO ACTUAL VALUE
        self.CL_max_clean = 1.5615             # Maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.CL_max_land = 1.9              # Maximum lift coefficient at landing [-]
        self.CLa = 4.743                      # Lift curve slope [-] | CHANGE TO ACTUAL VALUE
        self.c_p = 72E-9                    # Specific fuel consumption [kg/J]
        self.climb_rate = 2.9889
        self.cos_lambda_c04 = 1
        self.cruise_frac = 0.8348           # Assume halfway through the cruise with cruise fuel fraction 0.3 [-]
        self.CL_a_w = 4.743                  # Updated Lift curve slope [1/rad]


        # D
        self.Drag_increase = 1.0126         # This is used for the calculations of the strut drag if applicable [-]
        self.d_eff = 1.241                  # Effective diameter [m]
        self.d_engine_boxes = 0.4           # Distance between engine and wing box [m]
        self.dihedral = 0
        self.d_fuselage = 1.102                 # Fuselage diameter in [m]
        # E
        self.e = 0.776                        # Oswald factor [-]
        self.engine_cg = 0.267              # Engine cg location [m]
        self.engine_fairing = 0.2           # Engine fairing length [m]
        self.engine_length = 0.534         # Engine length [m]
        self.engine_pos = 'tractor'         # Engine position: "tractor" or "pusher" or "fuselage"
        self.eta_p = 0.82                   # Propulsive efficiency [-]
        # F

        # G
        self.g0 = 9.8066                    # Gravitational acceleration [m/s^2]
        # H
        self.h_TO = 0                       # Take-off altitude, airport altitude [m]
        self.h_cruise = 3048.0              # Cruise altitude [m]
        # self.h_in = 0.9                     # Inner fuselage height [m]
        # self.h_out = 0.737                    # Outer fuselage height [m]
        # I
        self.i_w = 0.935 * np.pi / 180       # Updated incidence angle of wing wrt fuselage [rad]


        # J

        # K
        self.kq = 0.95                      # Volume factor used to calculate wetted area of the wing [-]


        # L
        self.LDG_dist = 750.0               # Landing distance [m]
        self.L_D = 14.1804                  # Lift to drag ratio [-]
        self.Lambda = -0.0065
        self.l_fus_main_cone = 2.9          # Fuselage main cilindrical part length [m]
        self.l_fus_tail_cone = 0.8          # Fuselage tail cone length [m]
        self.l_f = 4.3                      # Fuselage length [m]
        self.l_f_boom = 2.5                   # Boom length [m]
        self.l_n = 0.9342                  # Nosecone length [m]
        self.l_t = 3.5                      # Tail arm [m]
        self.l_tc = 0.8                     # Tail cone length [m]
        self.lambda_co2 = -0.02736364725188103           # Half chord sweep angle [rad]
        self.lambda_co4 = 0.0               # Quarter chord sweep angle [rad]
        self.lin_par1 = 0.5249              # [-]
        self.lin_par2 = 42.049              # [-]
        self.l_h = 4.15                        # [m] tail length; length of aerodynamic centre of wing to aerodynamic centre tail. NOTE: This is a design choice, so for now it is a guestimate.

        # M
        self.MAC_length = 1.2477198640078417            # Mean aerodynamic chord [m]
        self.MAC_ac = 0.24                  # Location of aerodynamic center relative to MAC [-]
        self.MGC = 1.0822                   # Mean geometric chord [m]
        self.M_res = 0.075 
        self.Mff = 0.9254                   # Fuel fraction [-]
        self.main_gear_type = 'fixed'       # Type of main gear: "fixed" or "retractable"
        self.mass_penalty_struts = 7        # The weight of the struts [kg]
        # N
        self.N_e = 1                        # Number of engines [-]
        self.n_boxes = 12                   # [-]
        self.n_boxes_abreast = 2            # [-]
        self.n_drops = 1                    # [-]
        self.n_rows = 6                     # [-]
        self.n_ult = 6.6                    # Ultimate load factor [-]
        self.name = name                    # Name of the aircraft
        self.nose_gear_type = 'fixed'       # Type of nose gear: "fixed" or "retractable"
        # O

        # P
        self.pos_main_carriage = 'fuselage' # Position of main carriage: "fuselage" or "wing"
        self.power = 95.8347                # Power at takeoff [hp]
        self.power_setting = 0.9            # Power setting in cruise [-]
        self.prop_eff = 0.82                 # Propulsive efficiency [-]
        # Q

        # R
        self.R = 500000                     # Range [m]
        self.R_gas = 287.05                 # Gas constant [J/kgK]
        self.RFL = 750.0                    # Required field length [m]
        self.rho0 = 1.225                   # Air density at sea level [kg/m^3]
        self.rho_TO = 1.225                 # Take-off air density if airport is at sea level [kg/m^3]
        self.rho_cruise = 0.9046            # Cruise air density [kg/m^3]
        self.rootchord = 1.490             # Root chord [m]
        # S
        self.S_G = 19.77                    # Gross shell area fuselage [m^2]
        self.Sh_S = 0.174                    # [-] Ratio between horizontal tailplane surface area and surface area wing
        self.Sv_S = 0.083                  # [-] Ratio between vertical tailplane surface area and surface area wing
        self.Sw = 11.7113                   # Wing area [m^2]
        self.Sw_wetted = 23.4226            # Wetted area of the wing [m^2]
        self.s_tail = 2                     # Tail surface area [m]
        self.side_clearance = 0.2           # Side clearance [m], this is for both sides
        self.sigma_TO = 1.0
        self.sigma_cruise = 0.7385
        self.structural_thickness = 0.2     # Structural thickness fuselage [m], this is for both sides
        self.sweep_co2 = -0.02736364725188103 # Updated half chord sweep angle [rad]
        self.sweep_co4 = 0.0                 # Updated half chord sweep [rad]

  
        # T
        self.T0 = 288.15                    # Sea level temperature [K]
        self.TOP_req = 250 
        self.t_c = 0.15                     # Thickness to chord ratio [-]
        self.taper = 0.65                   # Taper ratio [-]
        self.tipchord = 0.9685              # Tip chord [m]
        self.tire_nose_height = 0.35306     # Nose landing gear tire height [m]
        self.tire_nose_width = 0.14478      # Nose landing gear tire width [m]
        self.tire_main_height = 0.381       # Main landing gear tire height [m]
        self.tire_main_width = 0.1524       # Main landing gear tire height [m]
        self.tipchord = 0.9685             # Tip chord [m]
        self.top_clearance = 0.2            # Top clearance [m]
        self.type = "utility"               # CS23 aircraft type: "normal" for normal/commuter and "utility" for utility
        # U
        self.V_A = 47.9382                  # [m/s]
        self.V_B = 45.0                     # [m/s]
        self.V_D = 90.7328                  # [m/s]
        self.V_climb = 36.0111              # [m/s]
        self.V_cruise = 60.4885             # Cruise velocity [m/s]
        self.V_s_min = 22.8537              # Minimum stall velocity [m/s]
        # W
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
        self.W_TO = 752                    # Take-off weight [kg]
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
        self.w_in = 1.2                     # Inner fuselage width [m]
        self.w_out = 1.1                    # Outer fuselage width [m]
        self.wing_twist = -2.0 *np.pi/180   # Updated wing twist (difference root and chord) [rad]

        # X
        self.xc_aft_spar = 0.8              # Aft spar location as fraction of MAC
        self.X_LEMAC = 2.276                # Leading edge mean aerodynamic chord [m]
        self.x_lemac = 0.06057988483270884  # Distance from LE root chord to the leading edge mean aerodynamic chord [m]
        self.xc_OEW_p = 0.25                # Center of gravity of OEW as a fraction of the MAC [-]
        self.X_cg_aft = 0.5335              # Aft cg location CG/MAC [-]
        self.X_cg_full = 0.4115             # MTOW cg location CG/MAC [-]
        self.X_cg_fwd = 0.1704              # Forward cg location CG/MAC [-]
        self.X_cg_range = 0.363             # Range of cg location CG/MAC [-]
        # Y
        self.y_mac = 2.2133293637093265     # Spanwise location of the MAC [m]
        self.ystart_ail = 3.308             # start location of aileron measured from rootchord
        self.yend_ail = 3.975               # end location of aileron measured from rootchord
        self.yend_flap = 1.5629             # end location of the HLD measured from rootchord


        "Structural parameters"             # NOTE: Add identifier "ST_" before variable names
        self.something = 1 # add units
        #self.y_mac = 2                      # Spanwise location of the MAC [m]
        self.ST_SF = 1.5

        "Aerodynamic parameters"            # NOTE: Add identifier "AE_" before variable names
        "Main wing and overall a/c"
        self.AE_A = 7.75                        # Updated aspect ratio [-]
        self.AE_CD0 = 0.02578                   # Zero lift drag [-]
        self.AE_CL_LDG = 1.5702                 # Still to be updated [-]
        self.AE_CL_max_TO = 1.5                 # Still to be updated maximum lift coefficient at take-off [-]
        self.AE_CL_max_clean = 1.5615              # Still to be updated maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.AE_CL_max_land = 1.9               # Still to be updated maximum lift coefficient at landing [-]
        self.AE_CL_a_w = 4.743                  # Updated Lift curve slope [1/rad]
        self.AE_L_D = 14.1804                   # Still to be updated lift to drag ratio [-]
        self.AE_MAC_length = 1.2477198640078417 # Updated mean aerodynamic chord [m]
        self.AE_MAC_ac = 0.24                   # Updated location of aerodynamic center relative to MAC [-]
        self.AE_Sw = 11.7113                    # Updated wing area [m^2]
        self.AE_Sw_wetted = 23.4226             # Updated wetted area of the wing [m^2]
        self.AE_b = 9.527                       # Updated wing span [m]
        self.AE_dihedral = 0                    # Updated wing dihedral angle [rad]
        self.AE_span_eff = 0.9940               # Updated Span eficiency factor (different from oswald) [-]
        self.AE_tau = 1                         # Factor used in prandtl glauert correction for airfoil curve slope to wing [-]
        self.AE_e = 0.7778                       # Updated oswald efficiency factor [-]
        self.AE_i_w = 0.935 * np.pi / 180       # Updated incidence angle of wing wrt fuselage [rad]
        self.AE_wing_twist = -2.0 *np.pi/180    # Updated wing twist (difference root and chord) [rad]
        self.AE_sweep_co2 = -0.02736364725188103 # Updated half chord sweep angle [rad]
        self.AE_sweep_co4 = 0.0                 # Updated half chord sweep [rad]
        self.AE_sweep_LE = 0.02736364725188103  # Updated leading edge sweep [rad]
        self.AE_taper = 0.65                    # Updated taper ratio [-]
        self.AE_rootchord = 1.490               # Updated Root chord [m]
        self.AE_tipchord = 0.9685               # Updated tip chord [m]
        self.AE_x_lemac = 0.06057988483270884   # Updated distance from LE root chord to the leading edge mean aerodynamic chord [m]
        self.AE_y_mac = 2.2133293637093265      # Updated spanwise location of the MAC [m]
        self.AE_alpha_f = 0                     # Still to be updated angle of attack of the fuselage [rad]

        # Horizontal tailplane
        self.AE_l_h = 4                        # [m] tail length; length of aerodynamic centre of wing to aerodynamic centre tail. NOTE: This is a design choice, so for now it is a guestimate.
        self.AE_Vh_V = 0.95                    # Ratio between velocity at tail and wing [-] NOTE: This is a guestimate
        self.AE_A_h = 5.166666                        # Aspect ratio horizontal tail. NOTE: This is a guestimate  
        self.AE_dEpsilondA = 0.02              # Downwash [-] TODO: check this value, this is a pure guess
        self.AE_Sh_S = 0.174                   # [-] Ratio between horizontal tailplane surface area and surface area win
        self.AE_Sh = 2.0377662
        self.AE_CL_a_h = 4.18773706267545    # Lift curve slope horizontal tailplane [1/rad] 


        self.AE_taper_h = 1                
        self.AE_b_h = 3.649                     
        self.AE_i_w_h = -1.8957       
        self.AE_wing_twist_h = 0.0    
        self.AE_sweep_co4_h = 0.0                 # Updated half chord sweep [rad]
        self.AE_sweep_co2_h = 0.0 
        self.AE_sweep_LE_h = 0
        self.AE_rootchord_h = 0.931167            
        self.AE_tipchord_h = 0.931167       
        self.AE_MAC_length_h = 0.7061       
        self.AE_y_mac_h = 0.91225       #Constant chord, so, quarter of span is taken such that MAC is at half the halfspan
        self.AE_x_lemac_h = 0
        self.AE_lambda_co2_h = 0.0
        self.AE_horizontal_airfoil = '0012'      # Airfoil of horizontal tail (NACA)


        # Vertical tailplane
        self.AE_Vv_V = 1                       # [-] Ratio betweeen velocity at vertical tail and free-stream velocity
        self.AE_A_v = 2                     # [-] Aspect ratio vertical tail
        self.AE_lambda_c02_v = None            # [rad] Half chord sweep of vertical tailplane
        self.AE_Sv_S = 0.1095                  # [-] Ratio between vertical tailplane surface area and surface area wing
        self.AE_Sv = 1.2824
        self.AE_b_v = 1.2015
        self.AE_vertical_airfoil = '0009'      # Airfoil of vertical tail (NACA)
        self.AE_rootchord_v = 0.9530
        self.AE_tipchord_v = 0.6671
        self.AE_lambda_LE_v = 35 * np.pi / 180
        self.AE_i_w_v = 0.07

        "-NACA4415"
        self.airfoil = "4415"

        if self.airfoil == "4415":
            self.af_Cl0 = 0.457                                 # TODO: Change to real value - Lift coeff of airfoil @ 0 AOA, cruise velocity [-]
            self.af_clcd_max = 163.5                            # Maximum clcd
            self.af_clcd32_max = 170.1                          # Maximum clcd**(3/2)
            self.af_clcd12_max = 165.7                          #
            self.af_cl_max = 1.735                              # Maximum cl_max
            self.af_alpha_s = 18.0 * np.pi / 180                # Stall angle of attack
            self.af_cd0 = 0.00595                               # Drag coefficient at zero lift
            self.af_cl_alpha = 0.103 * 180 / np.pi              # Lift curve slope [1 / rad]
            self.af_cm_alpha = 0.00748                          # Moment coefficient derivative [1/rad]
            self.af_cm0 = -0.0941                               # Moment coefficient at zero AoA
            self.af_alpha0 = -0.0774    # Angle of attack at zero lift

        if self.airfoil == "clarky":
            self.af_clcd_max = 154.7
            self.af_clcd32_max = 150.9
            self.af_clcd12_max = 161.7
            self.af_cl_max = 1.786
            self.af_alpha_s = 16.5 * np.pi / 180
            self.af_cd0 = 0.00604
            self.af_Cl0 = 0.404
            self.af_cl_alpha = 0.113 * 180 / np.pi
            self.af_cm_alpha = 0.00627
            self.af_cm0 = -0.0844
            self.af_alpha0 = -self.AE_Cl0 / self.AE_cl_alphaaf

        "Flight Performance parameters"         # NOTE: Add identifier "FP_" before variable names
        self.screenheight   = 50*0.3048         # screen height of 50 ft (CS23)
        self.rpm_maxcont    = 5500              # rpm
        self.omega_prop     = 237               # rad/s, based on 5500 rpm max continuous power and 2.43 gearbox ratio
        self.prop_radius    = 0.7237         # [m] based on 3 blade rotax 3B0 ground adjustable propeller by sensenich propellers
        self.ceiling        = 18000*0.3048      # [m] 18000 ft service ceiling
        self.th_ceil        = 30000*0.3048
        self.SFC            = 1.2*7.77777778e-8 # kg/J specific fuel consumption
        self.fuelcapacity   = 100               # L
        self.fueldensity    = 0.7429            # kg/L
        self.turnrate_half  = 1.5               # deg/s
        self.turnrate_1     = 3.0               # deg/s
        self.turnrate_2     = 6.0               # deg/s
        self.accelheight    = 300*0.3048
        self.FP_CL_max_land = 2
        self.FP_CL_max_to = 1.5
        self.FP_CL_land = 0.83
        self.FP_CL_to = 1.1


        "Control and stability parameters"  # NOTE: Add identifier "CS_" before variable names
        self.CS_eta = 0.95                  # airfoil efficiency factor [-]
        # self.CS_cf = 0.3                    # flap chord [m] TODO: update value
        self.CS_mu1 = 0.24
        self.CS_mu2 = 0.78
        self.CS_mu3 = 0.525
        self.CS_x_ac_w = 0.24              # location of wing ac, divided by MAC [-] - SEAD L7, S34
        self.CS_Cm_0_airfoil = -0.083       # TODO: Update value - Moment coefficient of airfoil [-]
        self.CS_n_blades = 3                   # [-] number of propeller blades NOTE: Depends on chosen propeller
        # self.CS_D_prop = 1.75                  # [m] Diameter of propeller NOTE: Depends on chosen propeller
        "Operations parameters"             # NOTE: Add identifier "OP_" before variable names
        # inputs
        self.OP_fuel_energy_density = 44.65E6 # [J/kg]
        self.OP_h_loiter = 500 * 0.3048  # [m]
        self.OP_h_dropzone = 0 # [m]
        self.OP_n_app_max = 3 # [g]

        # requirements
        self.OP_hmin = 15  # [m]              # requirements
        self.OP_b_dropzone = 25 # [m]
        self.OP_l_dropzone = 25  # [m]
        self.OP_V_crosswind = 10  # [m/s]
        self.OP_V_tailwind = 15  # m/s]
        self.OP_V_headwind = 15  # [m/s]
        self.OP_V_wind = max(self.OP_V_headwind, self.OP_V_crosswind,self.OP_V_tailwind)  # [m/s]

        self.OP_Range = 250 # [km]
        self.OP_N_boxes_per_sortie = 12  # [-]
        self.OP_MR_PL = 20000 # [kg/day]
        self.OP_PL_per_box = 20 # [kg]
        self.OP_TTFD = 72 # [h]

        # box drop maneuver
        self.OP_V_boxlim = 100 / 3.6 # [m/s] box drop max speed
        self.OP_Vbox_LDG = 40 / 3.6  # [m/s] 40kph drop limit
        self.boxDX = 0.5  # [m]
        self.boxDY = 0.3  # [m]
        self.boxDZ = 0  # [m]

        # cost breakdown inputs
        self.OP_T_ops = 28 # [days]
        self.OP_N_ops = 659 # [operations]
        self.OP_AC_per_op = 20 # [#AC] available on average
        self.OP_n_drops = 2 # [#] choice!
        self.OP_TTFS = 66.712 # [h] from contract to finished assembly and first sortie starts
        self.OP_T_sortie_gnd = 2.067 # [hr]
        self.OP_T_ground = 2 # [h]
        self.OP_T_pilot = 0.75 # [h]

        self.OP_CST_nofuel = 217584726.5 # [euro]
        self.OP_fuelprice = 1.1935 # [euro/L]

        self.OP_T_taxi = 5/60 # [h]
        self.OP_T_TO = 10 / 60  # [h]
        self.OP_T_LDG = 10 / 60  # [h]
        self.OP_T_clearance = 5 / 60  # [h]

        "Tim's coefficients:"

        self.lift_coefficients = [1.02549983e+03,  3.40489460e+02, -2.62539473e+03,  6.99194179e+03,     #Coefficents of a polynomial fit for the
                            -9.67481071e+03,  7.76224418e+03, -3.81444542e+03,  1.16482541e+03,     #lift distribution over the half span
                            -2.15539510e+02,  2.21217522e+01, -9.66193857e-01]                      #Highest order coefficient first

        "Jan W's coefficinets:"

        self.ST_U_de = 50 #derived gust velocity (ft/s)
        self.ST_n_nw = 2.25 #load factor for nose wheel load
        self.ST_n_imp = 3.0 #impact inertia load factor
        self.ST_n_LW = 2/3 #L/W at bad landing
        self.ST_n_m = 3.8 # positive limit maneuvering load factor (from Vn)
        self.ST_n_ult_pos = 6.6 #positive ultimate load factor
        self.ST_n_ult_neg = -2.78 #negative ultimate load factor

        self.ST_Torque_eng = 128 #Nm, Rotax 912 torque
        self.ST_Thrust_eng = 2800 #N Rotax 912 thrust
        self.ST_W_eng = 65.7 #kg Rotax 912 plus clutch, exhaust, alternator, air guide hood, NO MOUNT
        self.ST_Torque_eng2 = 240 #Nm, UL 260 torque

        self.ST_l_strut = 2.563733275367119 #m strut lenght (single)
        self.ST_d_strut = 0.015 #strut diameter (outer) [m]
        self.ST_l_LG = 0.4 #m strut lenght (single)
        self.ST_d_LG = 0.04 #strut diameter (outer) [m]
        self.ST_d_boom = 0.05 #strut diameter (outer) [m]
        self.ST_y_strut = 2.0829580913074834  #spanwise location of strut attachment

        self.ST_W_fus = 40.224 #kg mass of fuselage structure (only)
        self.ST_W_tb = 11.52 #kg mass of tail boom, given it is 2.8 m long
        self.ST_W_lg = 7.75 #kg mass of landing gear struts

        self.ST_z_ground = 0.5 #m floor height
        self.ST_z_prop = 0.3 #m propeller clearance
        self.ST_h_fus = 0.67 #m fuselage height ##KEEP IN MIND THAT THE ATTACHMENT OF THE WING IS ABOVE THE FUSELAGE CEILING
        self.ST_w_fus = 1.0 #m fuselage width
        self.ST_h_prop_axis = 1.1 #m height of the propeller axis above the ground
        self.ST_top_eng = 1.87 #m ceiling of the engine above the ground
        self.ST_eng_w = 0.56 #m width of engine box
        self.ST_eng_h = 0.335 #m height of engine box (NOT INCLUDING PROP obviously)
        self.ST_eng_l = 0.534 #m lenght of engine (excluding gas exhaust piping)

        self.ST_z_cg_ground =  0.92 #m estiamted center of gravity of boxes, fuselage, engine, boom, wing, fuel
        self.ST_x_nw =0.4 #m x distance nose wheel
        self.ST_x_mw =2.35 #m x distance main whee
        self.ST_ax_g = 0.51 #- maximum horizontal breaking deceleration scaled by g0
        self.ST_x_cg = 2.94 #max aft x_cg for nose wheel loading


        self.ST_tail_boom_2a=0.15 #height of tail boom
        self.ST_tail_boom_2b=0.1 #width of tail boom
        self.ST_tail_boom_t=0.001 #thinckness of tail boom (steel 410)
        self.ST_h_tail_spar_r=0.025 #radius of h tail spar
        self.ST_h_tail_spar_t=0.001 #thickness of h tail spar (steel 410)


class airport:
    def __init__(self, name):
        self.name       = name
        self.mu_ground = 0.08
        self.mu_ground_break  = 0.2                  #buildingspeed.org
        self.rwyslope   = 1.0                   #runway slope in degrees --> CONVERT TO RADIANS

class atmosphere:
    def __init__(self):
        self.rho0   = 1.225     # kg/m^3
        self.lambd  = -0.0065   # troposphere, deg K/m
        self.T0     = 288.15    # K
        self.g      = 9.80665   # m/s^2
        self.R      = 287.05    
        self.p0     = 101325    # Pa 
        self.gamma  = 1.4

