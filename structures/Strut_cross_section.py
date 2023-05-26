import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
from loading_diagram_wing import Ay, Az, By, Bz, E, mat_density, sigma_yield, aircraft
import matplotlib.pyplot as plt
from scipy.integrate import quad

safety_factor = 1.5
truss_fraction = 0.437277492
L_truss = aircraft.wing_strut_location/np.cos(aircraft.strut_angle)
P_truss = np.sqrt(By**2 + Bz**2)
P_truss_comp = 1539*safety_factor #N
P_truss_tens = 6732*safety_factor #N
I_min = P_truss_comp * L_truss**2 / (np.pi**2 * E)
A_min = P_truss_tens / sigma_yield
a = np.linspace(0.001,0.0125,100) #semi_minor of ellipse [m]
b = 2*a #semi_major of ellipse [m]
t = 4*I_min/(np.pi*a**3 * (1+(3*b/a))) #thickness of strut sheet [m]
A = np.pi*(a+b)*t
I_check = np.pi*a**3*t*(1+(3*b/a))/4
m = A*mat_density*L_truss
print("a is",a[-1])
print("A is", A[-1])
print("t is",t[-1])
print('b is',b[-1])
print('m is', m[-1])
print('I_min is', I_min)
print('I_check is',I_check[-1])
print("strut length", L_truss)

plt.plot(a,[A_min]*100)
plt.plot(a, A)
#plt.show()