#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import voronizator
import polyhedron

voronoi = voronizator.Voronizator()

fig = plt.figure()
ax = fig.gca(projection='3d')

Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])


voronoi.setCustomSites(np.array([[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],[0.,1.,1.],[1.,0.,1.],[0.,0.,1.],[0.5,0.5,0.5]]))
voronoi.makeVoroGraph()
#voronoi.calculateShortestPath(Vs, Ve, 'near')
voronoi.calculateShortestPath(Vs, Ve, 'all')

voronoi.plotSites(ax)
voronoi.plotGraph(ax, pathExtremes=True)
voronoi.plotShortestPath(ax)

ax.set_xlim3d(0., 1.)
ax.set_ylim3d(0., 1.)
ax.set_zlim3d(0., 1.)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
