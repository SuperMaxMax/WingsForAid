"Parameters of aircraft"

n_drops = 1
R = 500 #km
M_res = 0.20
W_PL = 240
g = 9.80665

# Wing characteristics
W_TO = 2000 #lbs
S = 180 #ft^2
n_ult = 1.5
A = 3
e = 0.9
t_r = 2 #ft
b = 8 #ft
sweep_angle = 0.36
CD0 = 0.15

# Horizontal plane characteristics
S_h = 40
A_h = 4
t_rh = 3

# Vertical plane characteristics
S_v = 60
A_v =  2.5
t_rv = 6
labda_quarter = 0.4

prop_eff = 0.9
c_p = 9*10**-6

W1W_TO = 0.995          # engine start-up
W2W1 = 0.997            # taxi
W3W2 = 0.998            # take_off
W4W3 = 0.992            # climb
W10W9 = 0.993           # descent
WfinalW10 = 0.993       # landing, taxi, shut-down

