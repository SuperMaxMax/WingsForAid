"Physical parameters"
g = 9.80665  # m/s^2

"Parameters of aircraft"

# Mission flight plan
n_drops = 1
R = 500  # km
W_PL = 240
M_res = 0.15

# Wing characteristics
S = 180 #ft^2
n_ult = 1.5
A = 3
e = 0.9
t_r = 2 #ft
b = 8 #ft
sweep_angle = 0.36
CD0 = 0.15



prop_eff = 0.82
c_p = 90E-9

W1W_TO = 0.995          # engine start-up
W2W1 = 0.997            # taxi
W3W2 = 0.998            # take_off
W4W3 = 0.992            # climb
W10W9 = 0.993           # descent
WfinalW10 = 0.993       # landing, taxi, shut-down

lin_par1 = 0.5482 #gradient of the linear trend OEW/MTOW
lin_par2 = 486.68  #y axis crossing of the linear trend OEW/MTOW

import numpy as np
"Physical parameters"
g = 9.80665  # m/s^2

"Parameters of aircraft"

# Mission flight plan
n_drops = 1
R = 500  # km
W_PL = 240
M_res = 0.15
h_cruise = 10000*0.3048     #cruising at 10000 ft

# Wing characteristics
S = 180 #ft^2
n_ult = 1.5
A = 8
e = 0.9
t_r = 2 #ft
b = 8 #ft
sweep_angle = 0.36
CD0 = 0.15
CL_CD = 10 #conservative

prop_eff = 0.82
c_p = 90E-9

W_TO   = 900            # kg
W1W_TO = 0.995          # engine start-up
W2W1 = 0.997            # taxi
W3W2 = 0.998            # take_off
W4W3 = 0.992            # climb
W10W9 = 0.993           # descent
WfinalW10 = 0.993       # landing, taxi, shut-down
W_fuel_estimated = 75*0.82       # estimated value based on 20L/h fuel burn of rotax and 3 hour sortie and 15L reserve and density of fuel
W_LDG   = W_TO - W_PL - W_fuel_estimated
f = W_LDG/W_TO
cruise_frac = W1W_TO*W2W1*W3W2*W4W3*0.85       #assume halfway through the cruise with cruise fuel fraction 0.3

# Parameters for geometry estimation
# This file needs to be merged with parameters from Torben, which was not in Git by the time I started making this
# From ADSEE 1 slides, typical values
# Taking average of reported range of values  
CL_max_clean    = np.arange(1.3, 2.0, 0.1)   #1.3 - 1.9
CL_max_TO       = np.arange(1.3, 2.0, 0.1)
CL_max_land     = np.arange(1.6, 2.4, 0.1)
CL_TO           = CL_max_TO/(1.1**2)
CL_LDG          = CL_max_land/(1.1**2)
#parameter difference in different configurations
#ASSUMPTION: fixed undercarriage
d_CD0_TO_Flaps  = np.arange(0.010, 0.021, 0.001) #0.010 - 0.020
d_CDO_LDG_Flaps = np.arange(0.055, 0.076, 0.001) #0.055 - 0.075
d_e_TO_Flaps    = 0.05
d_e_LDG_Flaps   = 0.10

#atmospheric properties
p0      = 101325        #Pa
rho0    = 1.225         #kg/m^3
T0      = 288.15        #K
Lambda  = -0.0065       #deg K/m
R       = 287.05        #J/kgK
g0      = 9.80665       #m/s^2

#speeds
V_s_max = 61*(1.852/3.6)#CS23 Vs at take off not allowed to be above 61 kts, *1.852 to get to m/s
V_s_min = 50*(1.852/3.6)#Dropping speed
V_cruise= 105*(1.852/3.6) #Cruise speed
V_TO_max= 1.1*V_s_max   #Maximum take off speed
V_TO_min= 1.1*V_s_min   #Minimum take off speed
V_climb = 70*(1.852/3.6)


#power
P_max   = 100           #break horse power
eta_p   = 0.7           #propeller efficiency
power_setting = 0.85    #cruise power


#Take off distance
TO_dist = 750           #m, seems long, might need to change this, STOL gang
LDG_dist= 750
lin_par1 = 0.5482 #gradient of the linear trend OEW/MTOW
lin_par2 = 486.68  #y axis crossing of the linear trend OEW/MTOW