from matplotlib import pyplot as plt

label = [1, 5, 6, 7, 8, 10, 11, 12, 13, 17, 20, 21, 22, 23, 24, 25]
cost = [0.284027778, 0.000676448, 0.001503, 4.40059E-05, 0.00001088, 0.046, 0.046, 0.046, 0.018, 0.01, 0.0024345, 0.0060826, 0.0012962, 0.00445056, 0.0044059, 0.3]
payload = [1.8, 15422, 2120, 5000, 30000, 9, 6, 4.5, 100, 150, 1393, 458, 4000, 2286, 973, 240]
range = [80, 3330, 556, 600, 800, 27, 18, 15, 200, 300, 1982, 2424, 519, 439, 629, 500]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(cost, range, payload)

ax.set_xlabel('cost')
ax.set_ylabel('range')
ax.set_zlabel('payload')
plt.show()

from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

fig = plt.figure()