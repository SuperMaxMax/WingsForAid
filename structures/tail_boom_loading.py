import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from parameters import atmosphere
at=atmosphere()
from fuselage_truss_loading import dF_v_gust, dF_h_gust, yield_stress, E, density

a_b = dF_h_gust/dF_v_gust #ellipsoid cross sectio a/b optimal ratio
print(a_b)

M_max = 12220 #5569.2
F_max = 65865
shear_stress = 0.5*yield_stress
#shear_stress=165*10**6 #alu
L=2.8

a = 0.15/2
if a_b > 3/2:
    b = 2/3*a

Ixx_min = M_max*a/shear_stress
A_min = F_max/shear_stress
P = 2*np.pi*np.sqrt((a**2+b**2)/2)
print(Ixx_min,A_min)

t = np.linspace(0.0001,0.005,1000) #for multiple thichnesses
Ixx_min = [Ixx_min]*len(t)
A_min = [A_min]*len(t)
a1=a+t/2
a2=a-t/2
b1=b+t/2
b2=b-t/2

Ixx = 0.25*np.pi*(a1**3*b1-a2**3*b2) #define Ixx for given thickness
A = np.pi*(a1*b1-a2*b2) #define A for given thickness
sigma = M_max*a/Ixx #define max stress

plt.plot(t,Ixx, label="Ixx")
plt.plot(t,A, label='A')
#plt.plot(t,sigma, label="sigma")
plt.plot(t,Ixx_min, label='ixx_min')
plt.plot(t,A_min,label='A_min')
plt.legend()
plt.show()

t=2/1000
print("mass",P*t*L*density)

