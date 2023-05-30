import numpy as np
import matplotlib.pyplot as plt
#from loading_diagram_wing import aircraft

## Define material ##
# steel 4130
yield_stress = 460*10**6
density=7850
E = 205*10**9

### Define dimensions ###
L_AE=1.2
L_AB=0.8
L_AA_=1.2

## Define angles ##
theta = np.arctan((L_AE/2)/L_AB)
phy = np.arctan(L_AB/L_AA_)
#### Define loads ####
Ax=100
Az=100
Ax_=100
Ay_=100
Az_=100

Bx=100
Bz=100
Bx_=100
By_=100
Bz_=100

Cx=100
Cz=100
Cx_=100
Cz_=100

Dx=100
Dz=100
Dx_=100
Dy_=100
Dz_=100

Ex=100
Ez=100
Ex_=100
Ey_=100
Ez_=100

Fx=0
Fz=0
Fx_=0
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
for i in range(loads[index_tension]):
    mass = lenghts[i]*abs(loads)/yield_stress*density


mass = lenghts*abs(loads)/yield_stress*density #tension
mass = lenghts**2*density*np.sqrt(abs(loads)/(E*np.pi**2)) #compression
print(index_tension)
print(index_compression)
