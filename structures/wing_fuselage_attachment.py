import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")

# steel 410 parameters
sigma_yield = 1225*10**6 #Pa
sigma_yield_shear = 0.58*sigma_yield
density = 7800 #kg/m^3

# loads
Ay_max=32251.284804626564 #N
Ay_min= - 14930.31452182178 #N
Az_min = -759.779152495221#N
Az_max= 232.43169183082978 #N
R_pin = 0.01

#Pin failure

max_load = max(abs(Az_max),abs(Ay_max),abs(Az_min),abs(Ay_min))
print(max_load)
#R_pin = np.sqrt(max_load/(2*np.pi*sigma_yield_shear))
t_pin = max_load/(4*np.pi*R_pin*sigma_yield_shear)

print('thickness pin in mm',t_pin*1000)

#bearing shear out

alpha = 2
c = R_pin
rho_s = R_pin
Ksc = 1 + alpha*(c/rho_s)**0.5
R_notch_fus = 0.025
R_notch_w = 0.05

A_notch = abs(Ay_min)*Ksc/sigma_yield_shear

t_notch_fus = abs(Ay_min)*Ksc*0.5/(sigma_yield_shear*(R_notch_fus-R_pin))
print('thickness notsh fuselage for shear out mm',t_notch_fus*1000)
t_notch_w = abs(Ay_min)*Ksc/(sigma_yield_shear*(R_notch_w-R_pin))
print('thickness notch wing for shear out mm',t_notch_w*1000)

#net section failure

t_notch_fus_1 = abs(Ay_min)*Ksc*0.5/(2*sigma_yield*(R_notch_fus-R_pin))
t_notch_w_1 = abs(Ay_min)*Ksc/(2*sigma_yield*(R_notch_w-R_pin))
print('thickness notsh fuselage mm for section failure',t_notch_fus_1*1000)
print('thickness notch wing for section failure mm',t_notch_w_1*1000)

# compression failure yield

t_notch_fus_2 = abs(Ay_max)*Ksc*0.5/(2*R_pin*sigma_yield)
t_notch_w_2 = abs(Ay_max)*Ksc/(2*R_pin*sigma_yield)

t_notch_fus = max(t_notch_fus_1,t_notch_fus,t_notch_fus_2)
t_notch_w = max(t_notch_w,t_notch_w_1,t_notch_w_2)

print('t notch fus final mm',t_notch_fus*1000)
print('t notch wing final mm',t_notch_w*1000)
