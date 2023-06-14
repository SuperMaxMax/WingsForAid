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