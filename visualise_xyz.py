#!/usr/bin/python3
import sys
import numpy as np
from sklearn.neighbors.kde import KernelDensity

if len(sys.argv) < 2:
    print("%s file.xyz" % sys.argv[0])
    exit(42)

# Read file
file=open(sys.argv[1])

points = []
for line in file:
    vector = line.split()
    points += [[float(x) for x in vector]]

#Convert the points into a NP Array
X = np.array(points)
print (X)
# 3D Visualisation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)
xs = X[:,0]
ys = X[:,1]
zs = X[:,2]
print(xs)
ax.scatter(xs, ys, zs)

plt.show()



