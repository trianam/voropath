from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import scipy as sp
#import scipy.misc
#import math
import matplotlib.pyplot as plt
import scipy.spatial
#from scipy.spatial import Voronoi
#from scipy.spatial import ConvexHull

np.random.seed(0)

fig = plt.figure()
ax = fig.gca(projection='3d')

#V = np.array([[0.,0.,0.],[0.,1.,0.],[1.,1.,1.],[2.,1.,2.],[3.,2.,1.],[4.,3.,6.],[2.,3.,1.],[5.,3.,1.],[6.,5.,3.]])
V = sp.rand(10,3)

ax.plot(V[:,0], V[:,1], V[:,2], 'o')

vor = sp.spatial.Voronoi(V)
vorVer = vor.vertices
vorReg = vor.regions

print('vor.vertices:')
print(vorVer)
print('vor.regions:')
print(vorReg)
print('vor.ridge_vertices:')
print(vor.ridge_vertices)
print('vor.ridge_points:')
print(vor.ridge_points)

ax.plot(vorVer[:,0], vorVer[:,1], vorVer[:,2], 'o')
i = 0
for ver in vorVer:
    ax.text(ver[0], ver[1], ver[2], i, color='red')
    i = i+1

for ridge in vor.ridge_vertices:
    for i in range(len(ridge)):
        if (ridge[i] != -1) and (ridge[(i+1)%len(ridge)] != -1):
            ax.plot([vorVer[ridge[i],0], vorVer[ridge[(i+1)%len(ridge)],0]], [vorVer[ridge[i],1], vorVer[ridge[(i+1)%len(ridge)],1]], [vorVer[ridge[i],2], vorVer[ridge[(i+1)%len(ridge)],2]], 'k--', lw=2)
    
# for region in vorReg:
#     if region and (not -1 in region):
#         regVer = np.array(vorVer)[region]
#         hull = sp.spatial.ConvexHull(regVer)
#         for simplex in hull.simplices:
#             simplex = np.append(simplex, [simplex[0]])
#             ax.plot(regVer[simplex, 0], regVer[simplex, 1], regVer[simplex, 2], 'k-')

ax.set_xlim3d(0., 1.)
ax.set_ylim3d(0., 1.)
ax.set_zlim3d(0., 1.)
            
plt.show()
