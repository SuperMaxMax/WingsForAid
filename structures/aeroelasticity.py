import numpy as np
import matplotlib.pyplot as plt
import cmath

r1= 0.036 #bigger
r2 = 0.035 #smaler
P = 4000
c = 0.7
E = 205*10*9
L = 0.6
twist = 1/360*2*np.pi
m = 4
S = 2
I = 0.25*np.pi*(r1**4-r2**4)
b= c/2
kt = P/twist
kh = 3*E*I/L**3
e = 0.25*c
Cla = 2*np.pi
rho = 1.225
It = m*(0.5*np.sqrt(r1**2+r2**2))**2
St = m*0.25*c

u_d = np.sqrt(2*kt/(rho*S*e*Cla))
print("div speed",u_d)

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
    p3_new = cmath.sqrt(1/(2*a4)*(-a2-cmath.sqrt(a2**2-4*a4*a0)))
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
plt.plot(u,p1.real)
plt.plot(u,p2.real)
plt.plot(u,p3.real)
plt.plot(u,p4.real)
plt.show()

plt.plot(u,p1.imag)
plt.plot(u,p2.imag)
plt.plot(u,p3.imag)
plt.plot(u,p4.imag)
plt.show()