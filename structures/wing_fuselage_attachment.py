import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")

# s steel 410 parameters
sigma_yield = 1225*10**6 #Pa
sigma_yield_shear = 0.58*sigma_yield
density = 7800 #kg/m^3

# steel 4130 parameters
sigma_yield = 460*10**6
sigma_yield_shear = 0.58*sigma_yield
density = 7850
E=205*10**9
# check loads from drag

D_load = (ac.W_TO*ac.g0/10*ac.b/2)/(ac.rootchord/2)
print("D",D_load)


# loads
Ay_max=32251.284804626564 #+ D_load#N
Ay_min= - 14930.31452182178 #- D_load#N
Az_min = -759.779152495221#N
Az_max= 232.43169183082978 #N
Bz_min = -23142.254860289046 #n= 6.6
By_min = -32251.284804626564 #n= 6.6
B_min = np.sqrt(Bz_min**2 + By_min**2)
Bz_max = 10713.407106147533 #n=-2.78
By_max= 14930.31452182178
B_max = np.sqrt(Bz_max**2 + By_max**2)



R_pin = 0.02 # 0.005 for wing braces

#Pin failure

max_load = max(abs(Az_max),abs(Ay_max),abs(Az_min),abs(Ay_min))
max_load_B = max(B_min,B_max)
#print(max_load)

t_pin = max_load/(4*np.pi*R_pin*sigma_yield_shear)
#R_pin = np.sqrt(max_load_B/(2*np.pi*sigma_yield_shear))
#print('R_pin mm',R_pin*1000)
R_pin = 0.005

#print('thickness pin in mm',t_pin*1000)

# shear out

alpha = 2
c = R_pin
rho_s = R_pin
Ksc = 1 + alpha*(c/rho_s)**0.5
print("Ksc",Ksc)
R_notch_fus = 0.03
#R_notch_fus = 0.02
R_notch_w = 0.05

A_notch = abs(Ay_min)*Ksc/(2*sigma_yield_shear)

t_notch_fus = abs(Ay_min)*Ksc*0.5/(2*sigma_yield_shear*(R_notch_fus-R_pin))
#t_notch_fus = abs(B_max)*Ksc*0.5/(2*sigma_yield_shear*(R_notch_fus-R_pin))
print('thickness notsh fuselage for shear out mm',t_notch_fus*1000)
t_notch_w = abs(Ay_min)*Ksc/(2*sigma_yield_shear*(R_notch_w-R_pin))
print('thickness notch wing for shear out mm',t_notch_w*1000)

#net section failure

t_notch_fus_1 = abs(Ay_min)*Ksc*0.5/(2*sigma_yield*(R_notch_fus-R_pin))
#t_notch_fus_1 = abs(B_max)*Ksc*0.5/(2*sigma_yield*(R_notch_fus-R_pin))
t_notch_w_1 = abs(Ay_min)*Ksc/(2*sigma_yield*(R_notch_w-R_pin))
print('thickness notsh fuselage mm for section failure',t_notch_fus_1*1000)
print('thickness notch wing for section failure mm',t_notch_w_1*1000)

# compression bearing failure

t_notch_fus_2 = abs(Ay_max)*Ksc*0.5/(2*R_pin*sigma_yield)
#t_notch_fus_2 = abs(B_min)*Ksc*0.5/(2*R_pin*sigma_yield)
t_notch_w_2 = abs(Ay_max)*Ksc/(2*R_pin*sigma_yield)

t_notch_fus = max(t_notch_fus_1,t_notch_fus,t_notch_fus_2)
t_notch_w = max(t_notch_w,t_notch_w_1,t_notch_w_2)

print('t notch fus final mm',t_notch_fus*1000)
print('t notch wing final mm',t_notch_w*1000)


################ Some truss update ################





a_truss = 0.02      # semi_minor of ellipse [m]
b_truss = 2*a_truss          # semi_major of ellipse [m]
P = 2*np.pi*np.sqrt((a_truss**2+b_truss**2)/2) #perimeter

A_tens = B_min / sigma_yield
t_tens = P/A_tens
m_tens = A_tens*density*ac.ST_l_strut

I_min = B_max * ac.ST_l_strut**2 / (np.pi**2 * E)
#t = 4*I_min/(np.pi*a**3 * (1+(3*b/a)))  # thickness of strut sheet [m]
#A_comp = P*t
#I_check = np.pi*a**3*t*(1+(3*b/a))/4
#m_comp = A*mat_density*ac.l_strut

t = np.linspace(0.0001,0.005,1000) #for multiple thichnesses
I_min = [I_min]*len(t)
#Iyy_min = [Iyy_min]*len(t)
A_min = [A_tens]*len(t)
a1=a_truss+t/2
a2=a_truss-t/2
b1=b_truss+t/2
b2=b_truss-t/2

Ixx = 0.25*np.pi*(a1**3*b1-a2**3*b2) #define Ixx for given thickness
#Iyy = 0.25*np.pi*(a1*b1**3-a2*b2**3)
#Ixx = 0.25*np.pi*(a1**4-a2**4)
A = np.pi*(a1*b1-a2*b2) #define A for given thickness
#sigma = M_h_max*a/Ixx #define max stress

plt.plot(t,Ixx, label="Ixx")
#plt.plot(t,Iyy,label='Iyy')
plt.plot(t,A, label='A')
#plt.plot(t,sigma, label="sigma")
plt.plot(t,I_min, label='ixx_min')
#plt.plot(t,Iyy_min,label='iyy min')
plt.plot(t,A_min,label='A_min')
plt.legend()
plt.show()

t=1.5/1000
m_strut = density*ac.ST_l_strut*P*t
print("mass of strut",m_strut)



