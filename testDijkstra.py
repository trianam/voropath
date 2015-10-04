#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import voronizator
import polyhedron

voronoi = voronizator.Voronizator()

fig = plt.figure()
ax = fig.gca(projection='3d')

voronoi.addTestPolyhedron()
voronoi.makeTestGraph()
voronoi.calculateTestShortestPath()

voronoi.plotPolyhedrons(ax)
voronoi.plotGraph(ax, pathExtremes=True)
voronoi.plotShortestPath(ax)

# ax.set_xlim3d(0., 1.)
# ax.set_ylim3d(0., 1.)
# ax.set_zlim3d(0., 1.)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
