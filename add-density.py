#!/usr/bin/python3
import sys
import numpy as np
from sklearn.neighbors.kde import KernelDensity

if len(sys.argv) < 2:
    print("%s image.png [treshold]" % sys.argv[0])
    exit(42)

# Read file
file=open(sys.argv[1])

points = []
for line in file:
    vector = line.split()
    points += [[float(x) for x in vector]]

#Convert the points into a NP Array
X = np.array(points)
#print(X)

#The kernel aplied
kde = KernelDensity(kernel='gaussian', bandwidth=5.0).fit(X)
# Compute the density on each point of the sample
values = kde.score_samples(X)
vmin = min(values)
vmax = max(values)

i=0
for x, y in X:
   #Positive and renormalised log version
   print(x, y, (values[i] - vmin)/(vmax - vmin))
   #Exponential value
   #print(x, y, (np.exp(values[i]) - np.exp(vmin))/np.exp(vmax))
   i+=1
