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

