from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import scipy as sp
import scipy.misc
import math
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

np.random.seed(0)

fig = plt.figure()
ax = fig.gca(projection='3d')

#V = np.array([[0.,0.,0.],[0.,1.,0.],[1.,1.,1.],[2.,1.,2.],[3.,2.,1.],[4.,3.,6.],[2.,3.,1.],[5.,3.,1.],[6.,5.,3.]])
V = sp.rand(10,3)

ax.plot(V[:,0], V[:,1], V[:,2], 'o')

hull = ConvexHull(V)

print(hull.points)
print(hull.vertices)
print(hull.simplices)
print(hull.neighbors)
for simplex in hull.simplices:
    simplex = np.append(simplex, [simplex[0]])
    ax.plot(V[simplex, 0], V[simplex, 1], V[simplex, 2], 'k-')

# ax.set_xlim3d(0., 1.)
# ax.set_ylim3d(0., 1.)
# ax.set_zlim3d(0., 1.)
            
plt.show()
