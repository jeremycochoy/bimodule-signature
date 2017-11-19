#!/usr/bin/python3

import sys
import numpy as np
from collections import OrderedDict
from sortedcontainers import SortedList
from sortedcontainers import SortedDict
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
#distances = []

for line in file:
    vector = [float(x) for x in line.split()]
    coordinates = vector[:-1]
    density = vector[-1]
    points += [coordinates]
    densities += [density]

points = np.array(points)
densities = np.array(densities)
nb_pts = len(points)
max_density = max(densities)

#Reverse densities. This code is not needeed, we use it so that we can take
# the sup-levelsets for density, because our imput file have
# outliers with low densities and interesting points at high densities.
for i in range(0, len(densities)):
    densities[i] = max_density - densities[i]

###
## Utilitary function that allow accessing informations
## about simplexes (mostly edges)
###

#Fast distance computation (euclidean)
def distance(i, j):
    return np.linalg.norm(points[i] - points[j])

#Return the index, in our order, for the segment made of points
# i and j.
# Warning: It suppose i < j!
def seg_index(i, j):
    # This is computed in the first matrix pass...
    # Yes, it's an ugly hack. Should be fixed one day
    # by precalculating it at the begining of the program...
    return seg_index_map[(i, j)]
seg_index_map = {}
counter = 0
for i in range(nb_pts):
    for j in range(i + 1, nb_pts):
        seg_index_map[(i, j)] = counter


#Return the pair of time where the segment appear in the filtration.
def seg_time(i, j):
    #We compute the bifiltration index of the current simplex
    #stored in the tuple: (x, y)
    x = distance(i, j) # the two point should be close to each other: rips filtration
    y = max(densities[i], densities[j]) # the two points should be in the set
    return (x, y)

###
## The main algorithm that compute our matrices
###

### As in "Computing multidimensional persistence",
### we compute over Z2, where -1 = 1, so we forget about the sign.
### We also use the position over term order.

#Now we compute bondary matrix \delta_1
#It is stored in a sparse format (list of non zero coefs)
print("Compute transposed matrix From C_1 to C_0...")

d1 = []
for i in range(nb_pts):
    for j in range(i + 1, nb_pts):
        col = SortedDict()
        (x, y) = seg_time(i, j)
        #We simply compute the bondary opperator applied to this simplex.
        #Then, we store a polynomial coefficient in front of each
        #element of the cycle expression wich induce the right
        #grading degree of the expression ; the simplex connect in (x, y).
        #The output is x^distances * y^density = (distance, density)
        #In a rips filtration, all the points are present at the initial time
        col[i] = (x, y - densities[i])
        col[j] = (x, y - densities[j])
        d1 += [col]

print(d1)

#Then the bondary matrix \delta_2
print("Compute the transposed matrix From C_2 to C_1:")

d2 = []
for i in range(nb_pts):
    for j in range(i + 1, nb_pts):
        for k in range(j + 1, nb_pts):
            col = SortedDict()
            # the two point should be close to each other: rips filtration
            x = max(distance(i, j),
                    distance(i, k),
                    distance(j, k))
            # the two points should be in the set
            y = max(densities[i], densities[j], densities[k])
            #Remember seg_index(x, y) require x < y!
            (sx, sy) = seg_time(i, j)
            col[seg_index(i, j)] = (x - sx, y - sy)
            col[seg_index(j, k)] = (x, y)
            col[seg_index(i, k)] = (x, y)
            d2 += [col]
#print(d2)

########### Butcher implementation
###
### We use a slightly modified version of butcher to compute
### division on polynomial vectors instead of just polynoms.
###

# A vector has type SortedDict{line index : (x power, y power)}
def LM(vec):
    return vec.items()[0]

# In Z2, LT = ML :)
def LT(vec):
    return LM(vec)

def LCM_poly(p, q):
    #print("p: ", p, " q: ", q)
    return (max(p[0], q[0]), max(p[1], q[1]))

# A vector has type SortedDict{line index : (x power, y power)}
def LCM(vec1, vec2):
    l1 = LM(vec1)
    l2 = LM(vec2)
    #print("LM Vec1:",l1)
    #print("LM Vec2:",l2)
    if l1[0] != l2[0]:
        return 0
    else:
        return LCM_poly(l1[1], l2[1])
#Debug: (should return (8, 10)
#print(LCM(SortedDict({0:(8, 8), 1:(10, 10), 2:(5, 5)}),
#          SortedDict({0:(0, 10), 1:(1, 10), 2:(2, 10)})))

# Divides the polynomial vector vec by
# the list of polynomial vectors veclist.
def DIVIDES(vec, veclist):
    p = vec
    r = SortedDict()
    q = {} #We use a dictionary for easy and fast acess to qi
    while len(p) != 0: #While p != 0
        we_did_something = false
        for i in range(0, len(veclist)):
            we_did_something = false
        if we_did_something == false:
            (lt_p_key, lt_p_value) = LT(p)
            #Z2 addition, and we use homogeneousness
            # to know that we won't have two different monomials at same time
            #--------- assert:
            printf("r monom:", r.get(lt_p_key, 0), " p lt monom:", lt_p_value)
            #---------
            if r.has_key(lt_p_key):
                del r[lt_p_key]
            else:
                r[lt_p_key] = lt_p_value
                del p[lt_p_key]
        
