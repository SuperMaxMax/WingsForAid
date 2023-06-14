import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from fuselage_truss_loading import dF_h_man
from parameters import atmosphere
at=atmosphere()

# # s-steel 410
yield_stress = 1225*10**6
yield_stress_shear = yield_stress*0.58
density=7800
E = 200*10**9


# # # steel 4130 parameters
# sigma_yield = 460*10**6
# sigma_yield_shear = 0.58*sigma_yield
# density = 7850
# E=205*10**9



v=0.33
G = E/(2*(1+v))


r_h_tail_spar = 0.035
l_h_tail_spar = ac.AE_b_h
F_h_max = abs(dF_h_man)
T_max = 0.25*ac.AE_MAC_length_h*max(0.5*at.rho0*ac.V_A**2*ac.AE_Sh*ac.AE_CL_a_h*(2/360*2*np.pi),0.5*at.rho0*ac.V_cruise**2*ac.AE_Sh*ac.AE_CL_a_h*(1/360*2*np.pi),F_h_max)
arm = ac.AE_b_h/2*1/3 #5/2*1/3 # conservative, assuming hor tail up to 5m
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

t=1./1000
print("mass kg",P*t*1.2*density)
print("area", 2*np.pi*r_h_tail_spar*t)
print('min area', A_min[0])
print('Ixx',1/8*np.pi*(r_h_tail_spar*2)**3*t)
print('min Ixx',Ixx_min[0])
print("T_max req",T_max)
print("T_max allowed", 2*t*yield_stress_shear*A_m)
print("Max twist, deg",(T_max/(4*A_m**2*G)*P/t*1.2/2)/(2*np.pi)*360)
#print('Iyy mm^4',0.25*np.pi*(a1*b1**3-a2*b2**3)*1000**4)
