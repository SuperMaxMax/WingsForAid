import numpy as np
import matplotlib.pyplot as plt

# read array from numpy file
weights = np.load('weights.npy')
n_stringers = np.arange(0, 31, 1)
n_ribs = np.arange(0, 31, 1)
n_stringers_g, n_ribs_g = np.meshgrid(n_stringers, n_ribs)

# interpolate the weights, to replace the 0 values
# get data points to interpolate, which are the nonzero weight, together with the corresponding n_stringers and n_ribs
data_points = []
for i in range(len(weights)):
    for j in range(len(weights[i])):
        if weights[i][j] != 0:
            data_points.append([n_stringers[i], n_ribs[j], weights[i][j]])

# interpolate the data points
from scipy.interpolate import griddata
data_points = np.array(data_points)
weights_2 = griddata(data_points[:, :2], data_points[:, 2], (n_stringers_g, n_ribs_g), method='linear')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(n_stringers_g, n_ribs_g, weights_2)
# plot the specific points in red
ax.scatter(data_points[:, 0], data_points[:, 1], data_points[:, 2], c='r', s=50)
ax.set_xlabel('Number of stringers')
ax.set_ylabel('Number of ribs')
ax.set_zlabel('Weight [kg]')
plt.show()

print(weights)