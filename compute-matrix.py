#!/usr/bin/python3
import sys
import numpy as np
from scipy import spatial

# Input: The input file should contain on each line
# a n+1 vector containing the coordinatesof a point in dimension n
# and a value between 0 and 1 coresponding to it's filtration by density.
#
# Output: The output are the bondary matrices between chain complex,
# where vector basis in dimension one are the points kept in order.
# In higher dimension, simplex are constructed in lexicographic order.
# For example, 3 points will give rise to the basis :
# (0), (1), (2), (0, 1), (0, 2), (0, 1, 2).
#
if len(sys.argv) < 2:
    print("%s points-and-density" % sys.argv[0])
    exit(42)

# Read file
file=open(sys.argv[1])
points = []
densities = []
distances = []

for line in file:
    vector = [float(x) for x in line.split()]
    coordinates = vector[:-1]
    density = vector[-1]
    points += [coordinates]
    densities += [density]

points = np.array(points)
densities = np.array(densities)
nb_pts = len(points)

#cKDTree is like KDTree but (supposedly) faster
tree = spatial.cKDTree(points)

#Since the points are in the tree, we look for the second neighbor
for p in points:
    ds, ps = tree.query(p, k=2)
    d = ds[1]
    p_idx = ps[1]
    distances += [d]
    #print("distance = %f, point_index = %d" % (d, p_idx))
distances = np.array(distances)
# Normalise distances
max_dist = float(max(distances))

for i in range(nb_pts):
    distances[i] = distances[i] / max_dist
    
#print(points)
#print(densities)
#print(distances)

#Now we compute bondary matrix \delta_1
print("Transposed matrix From C_1 to C_0:")
#matrix_d1 = np.empty([nb_pts, 0])

for i in range(nb_pts):
    for j in range(i + 1, nb_pts):
        col = np.array(np.zeros(nb_pts, dtype=np.int8), dtype='U20')
        #We compute the bifiltration index of the current simplex
        #stored in the tuple: (x, y)
        x = max(abs(distances[i] - distances[j]), distances[i],
                distances[j]) # the two point should be close to each other: rips filtration
        y = max(densities[i], densities[j]) # the two points should be in the set
        #Now we output the collumn of the matrix.
        #We simply compute the bondary opperator applied to this simplex.
        #Then, we output a polynomial coefficient in front of each
        #element of the cycle expression wich induce the right
        #grading degree of the expression ; the simplex connect in (x, y).
        #The output is x^distances * y^density = (distance, density)
        col[i] = "-x^%fy^%f" % (x - distances[i], y - densities[i])
        col[j] = "x^%f+y^%f" % (x - distances[j], y - densities[j])
        print(" ".join(col.tolist()))

        # Compute a numpy matrix with the values
        #col = np.array([col]).transpose()
        #matrix_d1 = np.concatenate((matrix_d1, col), axis=1)

#Print the Numpy Matrix
#np.set_printoptions(threshold=np.nan)
#print(matrix_d1)


#Then the bondary matrix \delta_2
print("Transposed matrix From C_2 to C_1:")

for i in range(nb_pts):
    for j in range(i + 1, nb_pts):
        for k in range(j + 1, nb_pts):
            col = np.array(np.zeros(nb_pts*(nb_pts - 1), dtype=np.int8), dtype='U20')
            #TODO
#            col[i] = "-x^%f-y^%f)" % (distances[i],densities[i])
#            col[j] = "x^%f+y^%f" % (distances[j],densities[j])
            print(" ".join(col.tolist()))
