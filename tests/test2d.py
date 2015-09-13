import numpy as np
import scipy as sp
import scipy.misc
import math
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import ConvexHull


V = sp.rand(100,2)
plt.plot(V[:,0], V[:,1], 'o')

vor = Voronoi(V)
vorVer = vor.vertices
vorReg = vor.regions
# print(type(vorVer)
# print(vorVer)
# print(type(vorReg))
# print(type(vorReg[0]))
# print(vorReg)

plt.plot(vorVer[:,0], vorVer[:,1], 'o')

#hull = ConvexHull(vorVer)
#for simplex in hull.simplices:
#    plt.plot(vorVer[simplex, 0], vorVer[simplex, 1], 'k-')


for region in vorReg:
    if region and (not -1 in region):
        regVer = np.array(vorVer)[region]
        hull = ConvexHull(regVer)
        for simplex in hull.simplices:
            plt.plot(regVer[simplex, 0], regVer[simplex, 1], 'k-')

            
voronoi_plot_2d(vor)

            
plt.show()


