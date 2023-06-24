import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from fuselage_truss_loading import dF_h_man
from parameters import atmosphere
at=atmosphere()

# # # s-steel 410
# yield_stress = 1225*10**6
# yield_stress_shear = yield_stress*0.58
# density=7800
# E = 200*10**9


# # # steel 4130 parameters
yield_stress = 460*10**6
yield_stress_shear= 0.58*yield_stress
density = 7850
E=205*10**9



v=0.33
G = E/(2*(1+v))


r_h_tail_spar = 0.035
l_h_tail_spar = 3.09*1/3+0.2 #changed for final value
F_h_max = 1987 #abs(dF_h_man)
T_max = 0.25*F_h_max * 0.58 #ac.AE_MAC_length_h
arm = 1/3*0.5*3.09 #ac.b_h #5/2*1/3 # conservative, assuming hor tail up to 5m
M_h_max = abs(F_h_max*arm)
A_m = np.pi*r_h_tail_spar**2


Ixx_min = M_h_max*r_h_tail_spar/yield_stress
A_min = F_h_max/yield_stress_shear
t = np.linspace(0.0001,0.005,1000) #for multiple thichnesses
Ixx_min = [Ixx_min]*len(t)
#Iyy_min = [Iyy_min]*len(t)
A_min = [A_min]*len(t)
r1=r_h_tail_spar+t/2
r2=r_h_tail_spar-t/2
P=2*np.pi*r_h_tail_spar

Ixx = 1/8*np.pi*(r_h_tail_spar*2)**3*t #define Ixx for given thickness
#Iyy = 0.25*np.pi*(a1*b1**3-a2*b2**3)
#Ixx = 0.25*np.pi*(a1**4-a2**4)
A = 2*np.pi*r_h_tail_spar*t #define A for given thickness
#igma = M_h_max*a/Ixx #define max stress

plt.plot(t,Ixx, label="Ixx")
#plt.plot(t,Iyy,label='Iyy')
plt.plot(t,A, label='A')
#plt.plot(t,sigma, label="sigma")
plt.plot(t,Ixx_min, label='ixx_min')
#plt.plot(t,Iyy_min,label='iyy min')
plt.plot(t,A_min,label='A_min')
plt.legend()
plt.show()


t=0.6/1000
print("mass kg",P*t*l_h_tail_spar*density)
print("area", 2*np.pi*r_h_tail_spar*t)
print('min area', A_min[0])
print('Ixx',1/8*np.pi*(r_h_tail_spar*2)**3*t)
print('min Ixx',Ixx_min[0])
print("T_max req",T_max)
print("T_max allowed", 2*t*yield_stress_shear*A_m)
print("Max twist, deg",(T_max/(4*A_m**2*G)*P/t*l_h_tail_spar/2)/(2*np.pi)*360)
print("lenght",l_h_tail_spar )
#print('Iyy mm^4',0.25*np.pi*(a1*b1**3-a2*b2**3)*1000**4)
