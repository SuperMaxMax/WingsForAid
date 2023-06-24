import numpy as np
import matplotlib.pyplot as plt
import cmath

r1= 0.0353 #bigger
r=0.035
r2 = 0.0347 #smaler
A_m = np.pi*r**2
P = 1987/2
c = 0.58
e = 0.25*c
E = 200*10**9
v=0.33
G = E/(2*(1+v))
L = 6*1.23/2
twist = (P*e/(4*A_m**2*G)*2*np.pi*r/(0.6/1000)*L) #0.81/360*2*np.pi
# print(P*e)
# print((4*A_m**2*G))
# print(2*np.pi*r/(0.6/1000)*L)
# print("twist", twist/(2*np.pi)*360)
m = 6.65/2
S = 1.79/2
I = 0.25*np.pi*(r1**4-r2**4)
b= c/2
kt = P/(2*twist)
kh = 3*E*I/L**3

Cla = 4.204357
rho = 1.225
It = m*(0.5*np.sqrt(r1**2+r2**2))**2
St = m*0.25*c

u_d = np.sqrt(2*kt/(rho*S*e*c*Cla))
print("div speed",u_d)
u_red = 54/0.58*np.sqrt(m/kh)
print("reduced velocity",u_red)

u = np.linspace(0,1000,1001)

p1 = np.array([])
p2 = np.array([])
p3 = np.array([])
p4 = np.array([])

for i in range(0,1001):

    q = 0.5 * rho * i ** 2
    a4 = m * It - St ** 2
    a2 = m * kt + It * kh - (2 * m * e * b + St) * q * S * Cla
    a0 = kh * (kt - 2 * e * b * q * S * Cla)
    p1_new = cmath.sqrt(1/(2*a4)*(-a2+cmath.sqrt(a2**2-4*a4*a0)))
    p2_new = -cmath.sqrt(1/(2*a4)*(-a2+cmath.sqrt(a2**2-4*a4*a0)))
    p3_new =  cmath.sqrt(1/(2*a4)*(-a2-cmath.sqrt(a2**2-4*a4*a0)))
    p4_new = -cmath.sqrt(1 / (2 * a4) * (-a2 - cmath.sqrt(a2 ** 2 - 4 * a4 * a0)))
    p1 = np.append(p1,p1_new)
    p2 = np.append(p2,p2_new)
    p3 = np.append(p3,p3_new)
    p4 = np.append(p4,p4_new)
    #print(p1)
# p1 = cmath.sqrt(1/(2*a4)*(-a2+cmath.sqrt(a2**2-4*a4*a0)))
# p2 = -cmath.sqrt(1/(2*a4)*(-a2+cmath.sqrt(a2**2-4*a4*a0)))
# p3 = cmath.sqrt(1/(2*a4)*(-a2-cmath.sqrt(a2**2-4*a4*a0)))
# p4 = -cmath.sqrt(1/(2*a4)*(-a2-cmath.sqrt(a2**2-4*a4*a0)))

# print(p1)

#print(a.real, a.imag)
plt.plot(u,p1.real,"b")
plt.plot(u,p2.real,"b")
plt.plot(u,p3.real,"b")
plt.plot(u,p4.real,"b",label="Real component")
plt.grid()
#plt.show()
#plt.ylim([-100, 100])

plt.plot(u,p1.imag,"r")
plt.plot(u,p2.imag,"r")
plt.plot(u,p3.imag,"r")
plt.plot(u,p4.imag,"r",label='Imaginary component')
#plt.ylim([-1, 100])
plt.xlabel('Air speed [m/s]',fontsize=15)
plt.ylabel('Solution of characteristic equation',fontsize=15)
plt.legend()
plt.show()
#plt.savefig('flutter.png',dpi=600)