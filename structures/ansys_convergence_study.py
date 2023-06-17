import numpy as np
import matplotlib.pyplot as plt

mesh = np.array([400,200,100,50,25,12.5,6.25,3.0])
deformation_tail_1 = np.array([72.104,51.437,44.575,41.054,39.939,39.754,39.759,39.848])
deformation_tail = deformation_tail_1/deformation_tail_1[-1]

plt.plot(mesh,deformation_tail,"o-", color = "orange")
plt.grid()
#plt.yscale("log")
plt.xscale("log")
plt.xlim(max(mesh), min(mesh))
plt.xlabel("Mesh element size [mm]",size=15)
plt.ylabel("Nondimentionalised node deformation [-]",size=15)
plt.show()

load = np.array([10,50,100,500,1000,5000,10000,50000])
E = 2*10**11
I = 7.6694*10**(-7)
L=1

anal_deformation =1000*load*L/(3*E*I)
anal_max_moment = load*L
#print(anal_deformation)
ansys_deformation = np.array([0.02307,0.11536,0.23072,1.1536,2.3072,11.536,23.072,115.36])
ansys_max_moment = 0.001*np.array([9750,48750,97500,4.875E5,9.75E5,4.875E6,9.75E6,4.875E7])

plt.plot(load,anal_deformation,"o-", color = "blue",label='Analytical')
plt.plot(load,ansys_deformation,"o-",color="orange", label='FEM')
plt.legend()
plt.grid()
plt.yscale("log")
plt.xscale("log")
#plt.xlim(max(load), min(load))
#plt.ylim(max())
plt.xlabel("Load [N]",size=15)
plt.ylabel("Absolute deformation [mm]",size=15)
plt.show()

plt.plot(load,anal_max_moment,"o-", color = "blue",label='Analytical')
plt.plot(load,ansys_max_moment,"o-",color="orange", label='FEM')
plt.legend()
plt.grid()
plt.yscale("log")
plt.xscale("log")
#plt.xlim(max(load), min(load))
#plt.ylim(max())
plt.xlabel("Load [N]",size=15)
plt.ylabel("Internal bending moment [Nm]",size=15)
plt.show()