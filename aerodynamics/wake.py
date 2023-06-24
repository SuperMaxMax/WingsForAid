import numpy as np
import matplotlib.pyplot as plt

rho = 1.225
v_inf = 60
D = 3
Rd = D/2
T = 10000


def axial(x, rho, v_inf, D, Rd, T):


    Tc = T / (rho*v_inf**2*D**2)

    a = -0.5 + 0.5*(1+8/np.pi*Tc)**0.5
    va = v_inf * a * (1 + x / (x**2+Rd**2)**0.5)

    vd = v_inf + va
    return vd


x = np.arange(0, 7, 0.1)
y = np.arange(-Rd, Rd, 0.1)
#z = np.matrix(len(x), len(y))

u1 = axial(x, rho, v_inf, D, Rd, T)
print(u1)
u = u1
for i in range(len(y)-1):
    u = np.vstack((u,u1))
print(u)

v = np.zeros(len(y))
v1 = v

for i in range(len(x)-1):
    v = np.vstack((v,v1))

#print(len(z))
#z.reshape(len(z), len(z))
print(len(x),len(y),len(u),len(v))


x = np.linspace(-1, 1, 21)
y = np.linspace(-1, 1, 21)


z = np.array([i*i+j*j for j in y for i in x])

print(len(z))

X, Y = np.meshgrid(x, y)
Z = z.reshape(21, 21)
Z = np.array([[1,1],
                        [3,3],
                        [4,4]])


x = np.arange(0, 4, 0.1)
y = np.arange(-Rd, Rd, 0.1)
u1 = axial(x, rho, v_inf, D, Rd, T)
Z = axial(x, rho, v_inf, D, Rd, T)
for i in range(len(y)):
    Z = np.vstack((Z, u1))


print(Z)

#plt.pcolor(X, Y, Z)
plt.imshow(Z, cmap = "jet", origin='lower',interpolation='bilinear', extent = [0,4,-Rd,Rd]) #Set map to jet or turbo

plt.show()


