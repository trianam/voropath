from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import scipy as sp
import scipy.misc
import math
import matplotlib.pyplot as plt
from pyhull.voronoi import VoronoiTess

V = np.array([[0,0,0],[0,1,0],[1,1,1],[2,1,2],[3,2,1]])
plt.gca(projection='3d')
plt.plot(V[:,0], V[:,1], V[:,2], 'o')

vor = VoronoiTess(V)
vorVer = np.array(vor.vertices)
print vorVer

plt.plot(vorVer[:,0], vorVer[:,1], vorVer[:,2])


plt.show()
