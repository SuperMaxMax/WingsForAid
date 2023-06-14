import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from parameters import atmosphere
at=atmosphere()
from fuselage_truss_loading import dF_v_gust, dF_h_gust, yield_stress, E, density, dF_h_man

a_b = dF_h_gust/dF_v_gust #ellipsoid cross sectio a/b optimal ratio
#print(a_b)



M_h_max = 5569.2 #5569.2  4364
F_h_max =  65865 #65865 7626
M_v_max = dF_h_man/a_b*1.8
F_v_max = dF_h_man/a_b
shear_stress = 0.5*yield_stress
#shear_stress=165*10**6 #alu
L=3.6658

L_EF = 0.8
L_FG = 1
L_GT=1
angle_1 = 40/360*2*np.pi
angle_2 = 20/360*2*np.pi
Gx = -F_h_max*L_GT
G = Gx/np.sin(angle_2)
Gz= G*np.cos(angle_2)
Fz=-Gz-F_h_max
F=Fz/np.sin(angle_1)
Fx = F*np.cos(angle_1)
Ex=-Fx-Gx

print('Fz',Fz)
print('Gz',Gz)

a = 0.15/2
# if a_b > 3/2:
#     b = 2/3*a
b=a

Ixx_min = M_h_max*a/shear_stress
Iyy_min = M_v_max*b/shear_stress
A_min = max(F_v_max/shear_stress, F_h_max/shear_stress)
P = 2*np.pi*np.sqrt((a**2+b**2)/2)
#print(Ixx_min,A_min)

t = np.linspace(0.0001,0.005,1000) #for multiple thichnesses
Ixx_min = [Ixx_min]*len(t)
Iyy_min = [Iyy_min]*len(t)
A_min = [A_min]*len(t)
a1=a+t/2
a2=a-t/2
b1=b+t/2
b2=b-t/2

Ixx = 0.25*np.pi*(a1**3*b1-a2**3*b2) #define Ixx for given thickness
Iyy = 0.25*np.pi*(a1*b1**3-a2*b2**3)
#Ixx = 0.25*np.pi*(a1**4-a2**4)
A = np.pi*(a1*b1-a2*b2) #define A for given thickness
sigma = M_h_max*a/Ixx #define max stress

plt.plot(t,Ixx, label="Ixx")
plt.plot(t,Iyy,label='Iyy')
plt.plot(t,A, label='A')
#plt.plot(t,sigma, label="sigma")
plt.plot(t,Ixx_min, label='ixx_min')
plt.plot(t,Iyy_min,label='iyy min')
plt.plot(t,A_min,label='A_min')
plt.legend()
plt.show()

t=2.5/1000
a1=a+t/2
a2=a-t/2
b1=b+t/2
b2=b-t/2
print("mass kg",P*t*L*density)
print("area mm^2", 2*np.pi*np.sqrt((a**2+b**2)/2)*t*1000**2)
print('Ixx mm^4',0.25*np.pi*(a1**3*b1-a2**3*b2)*1000**4)
print('Iyy mm^4',0.25*np.pi*(a1*b1**3-a2*b2**3)*1000**4)

#print('mass',t*2*np.pi*a*L*density)

