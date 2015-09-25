#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import voronizator

voronoi = voronizator.Voronizator()

fig = plt.figure()
ax = fig.gca(projection='3d')

Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])

voronoi.setRandomSites(10,0)
voronoi.makeVoroGraph()
voronoi.calculateShortestPath(Vs, Ve)

voronoi.plotSites(ax)
voronoi.plotShortestPath(ax)
voronoi.plotGraph(ax)

ax.set_xlim3d(0., 1.)
ax.set_ylim3d(0., 1.)
ax.set_zlim3d(0., 1.)
            
plt.show()
