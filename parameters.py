import numpy as np

class UAV:
    def __init__(self, name):
        "==== Aircraft Parameters ===="
        self.name                = name          # Name of the aircraft [-]

        "-CS23 Type"
        self.type                = "normal"      # CS23 aircraft type: "normal" for normal/commuter and "utility" for utility          

        "-Aircraft geometry"
        self.S                   = 180           # Surface area [m]
        self.A                   = 12            # Aspect ratio [-]
        self.e                   = 0.9           # Oswald factor [-]
        self.b                   = 8             # Wing span [m]
        self.sweep_angle         = 0.36          # Sweep angle [rad]
        self.lambda_mid          = 0.36          # Sweep angle at mid-wing [rad]
        self.t_c                 = 0.12          # Thickness over chord ratio [-]
        self.cwr                 = 2.5           # Chord length at root [m]

        self.b_f                 = 1.1           # Fuselage width [m]
        self.h_f                 = 1.1           # Fuselage height [m]
        self.l_f                 = 4             # Fuselage length [m]
        self.S_G                 = self.l_f*np.pi*self.b_f + 2*np.pi*(self.b_f/2)**2   # Gross shell area of fuselage [m] CHANGE THIS!

        self.s_tail              = 2             # Tail surface area [m]
        self.l_t                 = 3.5           # Tail arm [m]

        "-Aerodynamic properties"
        self.CD0                 = 0.027         # Zero lift coefficient [-]
        self.L_D               = 10            # Lift over drag [-] | ASSUMPTION/NOTES: Conservative

        # ASSUMPTION/NOTES: From ADSEE 1 slides, taking average of reported range of values  
        self.CL_max_clean        = 1.3           # Maximum lift coefficient [-] | Range: 1.3 - 1.9
        self.CL_max_TO           = 1.3           # Maximum lift coefficient at take-off [-]
        self.CL_max_land         = 1.6           # Maximum lift coefficietn at landing [-]
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

        self.W_fuel_estimated    = 75*0.82                                          # Fuel weight [kg] | ASSUMPTION: Estimated value based on 20L/h fuel burn of rotax and 3 hour sortie and 15L reserve and density of fuel
        self.W_LDG               = self.W_TO - self.W_PL - self.W_fuel_estimated    # Weight at landing [kg] | ASSUMPTION: All packages have been dropped
        self.fuel_frac           = self.W_LDG/self.W_TO                             # Fuel fraction [-]
        self.cruise_frac         = self.W1W_TO*self.W2W1*self.W3W2*self.W4W3*0.85   # Assume halfway through the cruise with cruise fuel fraction 0.3

        "-Propulsive properties"
        self.P_max               = 100           # Maximum power [bhp]
        self.P_TO                = 62            # Power at take-off [hp]

        self.prop_eff            = 0.82          # Propulsive efficiency [-]
        self.eta_p               = self.prop_eff # Propulsive efficiency [-]

        self.power_setting       = 0.9           # Power setting in cruise [-]

        self.c_p                 = 90E-9         # 
        self.N_e                 = 1             # Number of engines [-]

        "==== Mission profile/Atmospheric properties ===="
        "-Mission characteristics"
        self.n_drops             = 2             # Number of drops [-]
        self.R                   = 500000        # Range [m]
        self.M_res               = 0.15          # Fraction of remaing fuel at the end of flight/reserve fuel [-]
        self.h_cruise            = 10000*0.3048  # Cruise altitude [m] | NOTES: Conversion 

        self.TO_dist             = 750           # Take-off distance [m]           
        self.LDG_dist            = 750           # Landing distance [m]

        self.n_ult               = 1.5           # Ultimate load factor [-]

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

        self.lin_par1            = 0.5482        # Gradient of the linear trend OEW/MTOW
        self.lin_par2            = 486.68        # Y axis crossing of the linear trend OEW/MTOW

        "-MTOW vs OEW reduced by pilot weight R2=0.9548"
        #lin_par1           = 0.5522       
        #lin_par2           = -40.838 

        "-MTOW for drones, R2=0.9988"
        #lin_par1           = 0.4631
        #lin_par2           = 52.058

        "-MTOW vs OEW for general aviation R2=0.9548 , original y=0.5482 x + 486.68"
        #lin_par1           = 0.5522
        #lin_par2           = 39.162

        "- MTOW vs OEW ultra-light, R2=0.8704"
        #lin_par1           = 0.7134
        #lin_par2           = -52.981

        "- MTOW vs OEW ultra-light reduced by pilot weight, R2=0.9704"
        #lin_par1           = 0.7134
        #lin_par2           = -132.98

