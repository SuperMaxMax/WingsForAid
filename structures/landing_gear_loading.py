import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from parameters import atmosphere
at=atmosphere()
from fuselage_truss_loading import dF_v_gust, dF_h_gust, yield_stress, E, density, dF_h_man

masses = np.array([12*ac.boxweight, ac.ST_W_eng, ac.ST_W_tb, ac.ST_W_fus, ac.W_w, ac.W_F,ac.ST_W_lg])
z_cg = np.array([0.3,0.5061,0.67,0.67/2,0.67,0.67,0.35])

z_cg_tot = sum(z_cg*masses)/sum(masses)
print(z_cg_tot)
z_cg_ground = z_cg_tot + ac.ST_z_ground
print(z_cg_ground)






a = 0.03
b=a
r=0.03


Ixx_min = 5.92*10**(-8)
# Iyy_min = M_v_max*b/shear_stress
A_min = 1.94*10**(-5)
P = 2*np.pi*np.sqrt((a**2+b**2)/2)
#print(Ixx_min,A_min)

t = np.linspace(0.0001,0.005,1000) #for multiple thichnesses
Ixx_min = [Ixx_min]*len(t)
#Iyy_min = [Iyy_min]*len(t)
A_min = [A_min]*len(t)
a1=a+t/2
a2=a-t/2
b1=b+t/2
b2=b-t/2

#Ixx = 0.25*np.pi*(a1**3*b1-a2**3*b2) #define Ixx for given thickness
Ixx = 0.25*np.pi*(a1*b1**3-a2*b2**3)
#Ixx = 0.25*np.pi*(a1**4-a2**4)
A = np.pi*(a1*b1-a2*b2) #define A for given thickness
#sigma = M_h_max*a/Ixx #define max stress

plt.plot(t,Ixx, label="Ixx")
#plt.plot(t,Iyy,label='Iyy')
plt.plot(t,A, label='A')
#plt.plot(t,sigma, label="sigma")
plt.plot(t,Ixx_min, label='ixx_min')
#plt.plot(t,Iyy_min,label='iyy min')
plt.plot(t,A_min,label='A_min')
plt.legend()
# plt.show()

t=1/1000
a1=a+t/2
a2=a-t/2
b1=b+t/2
b2=b-t/2
#print("mass kg",P*t*L*density)
# print("area mm^2", 2*np.pi*np.sqrt((a**2+b**2)/2)*t*1000**2)
# print('Ixx mm^4',0.25*np.pi*(a1**3*b1-a2**3*b2)*1000**4)
# print('Iyy mm^4',0.25*np.pi*(a1*b1**3-a2*b2**3)*1000**4)

#print('mass',t*2*np.pi*a*L*density)


### tire loading ###
x_cg = ac.X_cg_aft*ac.AE_MAC_length+ac.X_LEMAC
l_m = ac.ST_x_mw-ac.ST_x_cg
l_n = ac.ST_x_cg-ac.ST_x_nw
P_n = ac.W_TO*ac.g0*((l_m+ac.ST_ax_g*ac.ST_z_ground)/(l_m+l_n))
#print('x_cg',ac.X_cg_full*ac.AE_MAC_length+ac.X_LEMAC)
print('nose wheel loading',P_n,l_m,l_n)

w = 0.9/ac.g0*ac.WS**0.25 #vertical speed
E = ac.W_TO*0.2*0.95/2*w**2 #energy to be absorbed by front wheel, if shock absorber used
