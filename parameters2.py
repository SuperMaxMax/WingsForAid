import numpy as np

"==== Aircraft Parameters ===="
"-Aircraft geometry"
S                   = 180           # Surface area [m]
A                   = 12            # Aspect ratio [-]
e                   = 0.9           # Oswald factor [-]
b                   = 8             # Wing span [m]
sweep_angle         = 0.36          # Sweep angle [rad]
lambda_mid          = 0.36          # Sweep angle at mid-wing [rad]
t_c                 = 0.12          # Thickness over chord ratio [-]
cwr                 = 2.5           # Chord length at root [m]

b_f                 = 1.1           # Fuselage width [m]
h_f                 = 1.1           # Fuselage height [m]
l_f                 = 4             # Fuselage length [m]
S_G                 = 4*np.pi*b_f   # Gross shell area of fuselage [m] CHANGE THIS!

s_tail              = 2             # Tail surface area [m]
l_t                 = 3.5           # Tail arm [m]


"-Aerodynamic properties"
CD0                 = 0.027         # Zero lift coefficient [-]
CL_CD               = 10            # Lift over drag [-] | ASSUMPTION/NOTES: Conservative

# ASSUMPTION/NOTES: From ADSEE 1 slides, taking average of reported range of values  
CL_max_clean        = 1.3           # Maximum lift coefficient [-] | Range: 1.3 - 1.9
CL_max_TO           = 1.3           # Maximum lift coefficient at take-off [-]
CL_max_land         = 1.6           # Maximum lift coefficietn at landing [-]
CL_TO               = CL_max_TO / (1.1**2) 
CL_LDG              = CL_max_land / (1.1**2)

# ASSUMPTION: Fixed undercarriage
d_CD0_TO_Flaps      = np.arange(0.010, 0.021, 0.001)    # Change in CD0 at take-off due to flaps [-] | Range: 0.010 - 0.020
d_CDO_LDG_Flaps     = np.arange(0.055, 0.076, 0.001)    # Change in CD0 at landing due to flaps [-] | Range: 0.055 - 0.075
d_e_TO_Flaps        = 0.05                              # Change in oswald factor at take-off due to flaps [-] 
d_e_LDG_Flaps       = 0.10                              # Change in oswald factor at landing due to flaps [-]


"-Weights"
W_e                 = 62.6          # Definitive weight per engine [kg]
W_G                 = 700           # Gross weight [kg]
W_TO                = 700           # Take-off weight [kg]

"-Weight fractions"
W1W_TO              = 0.995         # Engine startup fraction [-]
W2W1                = 0.997         # Taxi fraction [-]
W3W2                = 0.998         # Take_off fraction [-]
W4W3                = 0.992         # Climb fraction [-]
W10W9               = 0.993         # Descent fraction [-]
WfinalW10           = 0.993         # Landing, taxi & shut-down fraction [-]

W_fuel_estimated    = 75*0.82                           # Fuel weight [kg] | ASSUMPTION: Estimated value based on 20L/h fuel burn of rotax and 3 hour sortie and 15L reserve and density of fuel
W_LDG               = W_TO - W_PL - W_fuel_estimated    # Weight at landing [kg] | ASSUMPTION: All packages have been dropped
f                   = W_LDG/W_TO                        # Fuel fraction [-]
cruise_frac         = W1W_TO*W2W1*W3W2*W4W3*0.85    #assume halfway through the cruise with cruise fuel fraction 0.3


"-Propulsive properties"
P_max               = 100           # Maximum power [bhp]
P_TO                = 62            # Power at take-off [hp]

prop_eff            = 0.82          # Propulsive efficiency [-]
eta_p               = prop_eff      # Propulsive efficiency [-]

power_setting       = 0.9           # Power setting in cruise [-]

c_p                 = 90E-9         # 
N_e                 = 1             # Number of engines [-]



"==== Mission profile/Atmospheric properties ===="
"-Mission characteristics"
n_drops             = 2             # Number of drops [-]
R                   = 500000        # Range [m]
W_PL                = 240           # Payload weight [kg]
M_res               = 0.15          # Fraction of remaing fuel at the end of flight/reserve fuel [-]
h_cruise            = 10000*0.3048  # Cruise altitude [m] | NOTES: Conversion 

TO_dist             = 750           # Take-off distance [m]           
LDG_dist            = 750           # Landing distance [m]

n_ult               = 1.5           # Ultimate load factor [-]


"-Speeds"
V_s_max             = 61*(1.852/3.6)    # CS23 Vs at take off not allowed to be above 61 kts [m/s] | NOTES: *1.852 to get to m/s
V_s_min             = 50*(1.852/3.6)    # Dropping speed [m/s]
V_cruise            = 105*(1.852/3.6)   # Cruise speed [m/s]
V_TO_max            = 1.1*V_s_max       # Maximum take off speed [m/s]
V_TO_min            = 1.1*V_s_min       # Minimum take off speed [m/s]
V_climb             = 70*(1.852/3.6)    # Climb speed [m/s]
V_D                 = 150*0.514444      # Dive speed [m/s]


"-Atmospheric properties"
p0                  = 101325        # [Pa]
rho0                = 1.225         # [kg/m^3]
T0                  = 288.15        # [K]
Lambda              = -0.0065       # [deg K/m]
R_atm               = 287.05        # [J/kgK]
g0                  = 9.80665       # [m/s^2]



"==== Miscellanceous or however you spell it ===="

lin_par1            = 0.5482        #gradient of the linear trend OEW/MTOW
lin_par2            = 486.68        #y axis crossing of the linear trend OEW/MTOW

"-MTOW vs OEW reduced by pilot weight R2=0.9548"
#lin_par1           = 0.5522       
#lin_par2           = -40.838 

"-MTOW for drones, R2=0.9988"
#lin_par1           = 0.4631
#lin_par2           = 52.058

"-MTOW vs OEW for general aviation R2=0.9548 , original y=0.5482 x + 486.68"
#lin_par1           = 0.5522
#lin_par2           = 39.162

"- MTOW vs OEW ultra-light, R2=0.8704
#lin_par1           = 0.7134
#lin_par2           = -52.981

"- MTOW vs OEW ultra-light reduced by pilot weight, R2=0.9704
#lin_par1           = 0.7134
#lin_par2           = -132.98

