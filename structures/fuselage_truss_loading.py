import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from parameters import atmosphere
at=atmosphere()
#from loading_diagram_wing import aircraft

## Define material ##
# steel 4130
yield_stress = 460*10**6
density=7850
E = 205*10**9



# # # s-steel 410
# yield_stress = 1225*10**6
# density=7800
# E = 200*10**9

# # #al 5052
# yield_stress = 255*10**6
# density=2680
# E=70.3*10**9
## Define main forces ##
Thrust = 2800
Lift = ac.W_TO*ac.g0
Strut =4000
# print(Thrust,Lift,Strut)

### Define dimensions ###
L_AE=1.2
L_AB=0.8
L_AA_=1.2

## Define angles ##
theta = np.arctan((L_AE/2)/L_AB)
phy = np.arctan(L_AB/L_AA_)
#### Define loads ####
Ax=0.2*Thrust
Az=0.075*Lift
Ax_=Ax
Ay_=0
Az_=-0.25*Lift

Bx=0.2*Thrust
Bz=0.075*Lift
Bx_=Bx
By_=Strut
Bz_=0.475*Lift

Cx=0.1*Thrust
Cz=-0.05*Lift
Cx_=Cx
Cz_=0.25*Lift

Dx=0.2*Thrust
Dz=0.075*Lift
Dx_=Dx
Dy_=Strut
Dz_=0.475*Lift

Ex=0.2*Thrust
Ez=0.075*Lift
Ex_=Ex
Ey_=0
Ez_=-0.25*Lift

Fx=0.1*Thrust
Fz=0
Fx_=Fx
Fz_=0

### Define loads ###

AA_=Ax #L_AA_=L_AA_
BB_=Bx
L_BB_=L_AA_
CC_=Cx
L_CC_=L_AA_
DD_=Dx
L_DD_=L_AA_
EE_=Ex
L_EE_=L_AA_
FF_=Fx
L_FF_=L_AA_

AF=0
L_AF=L_AA_/2
EF=0
L_EF=L_AF
BC=0
L_BC=L_AF
CD=0
L_CD=L_AF

B_C_=By_
L_B_C_=L_AF
C_D_=Dy_
L_C_D_=L_AF


CF=Fz
L_CF=L_AB
C_F_=Fz_
L_C_F_=L_AB

AB=Bz #L_AB=L_AB
DE=Dz
L_DE=L_AB
A_B_=Bz_
L_A_B_=L_AB
D_E_=Dz_
L_D_E_=L_AB

AC=(-Az-Bz)/np.cos(theta)
L_AC=np.sqrt(L_AF**2 + L_AB**2)
CE=(-Ez-Dz)/np.cos(phy)
L_CE=L_AC

A_C_=(Az_+Bz_)/np.cos(theta)
L_A_C_=L_AC
C_E_=(Ez_+Dz_)/np.cos(theta)
L_C_E_=L_CE

A_F_=Ay_-A_C_*np.sin(theta)
L_A_F_=L_AF
E_F_=Ey_-C_E_*np.sin(theta)
L_E_F_=L_EF

CF_=(Cz-AC-CE-CF)/np.sin(phy)
L_CF_=np.sqrt(L_AB**2 + L_AA_**2)



loads = np.array([AA_,BB_,CC_,DD_,EE_,FF_,AF,EF,BC,CD,B_C_,C_D_,CF,C_F_,AB,DE,A_B_,D_E_,AC,CE,A_C_,C_E_,A_F_,E_F_,CF_])
lenghts = np.array([L_AA_,L_BB_,L_CC_,L_DD_,L_EE_,L_FF_,L_AF,L_EF,L_BC,L_CD,L_B_C_,L_C_D_,L_CF,L_C_F_,L_AB,L_DE,L_A_B_,L_D_E_,L_AC,L_CE,L_A_C_,L_C_E_,L_A_F_,L_E_F_,L_CF_])
mass = np.array([])

index_compression = np.where(loads<0)
index_tension = np.where(loads>=0)

for i in index_tension:
    mass_new = lenghts[i]*abs(loads[i])/yield_stress*density
    #print(mass_new)
    mass = np.append(mass, mass_new)
#print(mass)

for i in index_compression:
    mass_new = lenghts[i]**2*density*np.sqrt(abs(loads[i])/(E*np.pi**2))
    #print(mass_new)
    mass = np.append(mass,mass_new)
#print(mass)
# mass = lenghts*abs(loads)/yield_stress*density #tension
# mass = lenghts**2*density*np.sqrt(abs(loads)/(E*np.pi**2)) #compression
# print(index_tension)
# print(index_compression)

total_mass = sum(mass)
# print("loading",loads)
# print("masses",mass)
# print("Total modular fuselage mass",total_mass)


#### Panel fuselage calculations #####


q_h = (Lift/2)/L_AB
b=L_AA_
v=0.33
c_shear=15
c_comp=8
t_min1 = q_h*b**2*12*(1-v**2)/(c_shear*np.pi**2*E) #minimum thickness of side sheet (clamped both sides) for shear buckling
mass1= t_min1*b*L_AB*density
print('minimum mass of side sheet, shear buckling',mass1)

mass2=density*L_AA_*Thrust/yield_stress #minumim area of cross section for thrust tension
print('minimum mass for thrust tension',mass2)

mass3=density*L_AE*4000/yield_stress
t_min3=(Strut/yield_stress)/L_AA_
print('minimum mass for wing strut tension',mass3)

t_min4=(12*(ac.W_sc+ac.W_F)*ac.g0*L_AE*(1-v**2)/(c_comp*E*np.pi**2))**(1/3)
mass4=t_min4*density*L_AE*L_AB
print('minumum mass for compression due to wing structural weight',mass4)

t_min5=(12*0.5*(ac.W_sc+ac.W_F)*ac.g0*L_AA_*(1-v**2)/(c_comp*E*np.pi**2))**(1/3)
mass5=t_min5*density*L_AA_*L_AB
#print(0.5*(ac.W_sc+ac.W_F)*ac.g0)
print('minumum mass for compression of one side pannel due to wing structural weight',mass5)

t_min6=(12*1500*L_AA_*(1-v**2)/(c_comp*E*np.pi**2))**(1/3)
mass6=t_min5*density*L_AA_*L_AE
print('minumum mass for compression of bottom pannel due to wing structural weight',mass6)

#############################################################
## Some structural loading forces calculations

T_eng = 0

C_L_h = -0.1331864964576476
F_h_cruise=ac.W_TO*ac.g0*(ac.AE_Vh_V)**2*ac.Sh_S*C_L_h/ac.CL_max_clean
print("F_h in cruise",F_h_cruise)
rho = at.rho0
U_de = 50 #derived gust velocity (ft/s)
miu_g = 2*ac.WS*0.2/(rho*0.00194032033*ac.MAC_length*3.2808399*ac.CL_a_w*ac.g0*3.280839895013123)
Kg = 0.88*miu_g/(5.3+miu_g)#gust alleviation factor

V_a = ac.V_A /0.51444 # from m/s to kt
lb_to_N = 4.4482216153
m2_to_ft2 = 10.7639104
S_ht=ac.Sh_S*ac.Sw*m2_to_ft2
S_vt = ac.Sv_S*ac.Sw*m2_to_ft2
n_m = ac.ST_n_m

dF_h_man = n_m*ac.W_TO*ac.g0*((ac.X_cg_aft-ac.MAC_ac)*ac.MAC_length/ac.l_h - ac.Sh_S*(ac.AE_CL_a_h/ac.CL_a_w)*(1-ac.AE_dEpsilondA) - rho/2*(ac.Sh_S*ac.Sw*ac.AE_CL_a_h*ac.l_h/ac.W_TO)) # in N

dF_h_gust = Kg*ac.ST_U_de* ac.V_A/0.51444 *S_ht/498 *(1-ac.AE_dEpsilondA)* lb_to_N #in N

dF_v_gust  = Kg*ac.ST_U_de* ac.V_A/0.51444 *S_vt/498 * lb_to_N #in N

print("horizontal tail man load",dF_h_man)
print("hor tail gust load",dF_h_gust)
print("ver tail gust load",dF_v_gust)

# print('test',(ac.X_cg_aft-ac.MAC_ac)*ac.MAC_length/ac.l_h, ac.Sh_S*(ac.AE_CL_a_h/ac.CL_a_w)*(1-ac.AE_dEpsilondA) , rho/2*(ac.Sh_S*ac.Sw*ac.AE_CL_a_h*ac.l_h/ac.W_TO))
# print((ac.X_cg_aft-ac.MAC_ac)*ac.MAC_length)


