#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import voronizator
import polyhedron

voronoi = voronizator.Voronizator()

poly1 = polyhedron.Polyhedron(a = [0.1,0.1,0.1], b = [0.1,0.9,0.1], c = [0.9,0.1,0.1], d = [0.5,0.5,0.9])
poly2 = polyhedron.Polyhedron(a = [-0.5,0.1,0.1], b = [-0.5,0.9,0.1], c = [-0.9,0.5,0.1], d = [-0.5,0.5,1.5])
poly3 = polyhedron.Polyhedron(a = [0.5,1.,0.5], b = [0.7,0.5,0.5], c = [0.7,0.7,1.], d = [0.9,0.7,0.5])


fig = plt.figure()
ax = fig.gca(projection='3d')

Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])

voronoi.addPolyhedron(poly1)
voronoi.addPolyhedron(poly2)
voronoi.addPolyhedron(poly3)

#voronoi.addBoundingBox([-2.,-2.,-2.], [3.,3.,3.])
voronoi.setPolyhedronsSites()
voronoi.makeVoroGraph()
voronoi.calculateShortestPath(Vs, Ve)

#voronoi.plotSites(ax)
voronoi.plotPolyhedrons(ax)
#voronoi.plotGraph(ax, edges=False, labels=True)
voronoi.plotGraph(ax, pathExtremes=True)
#voronoi.plotGraph(ax)
voronoi.plotShortestPath(ax)

ax.set_xlim3d(-1., 2.)
ax.set_ylim3d(-1., 2.)
ax.set_zlim3d(-1., 2.)

# ax.set_xlim3d(0., 1.)
# ax.set_ylim3d(0., 1.)
# ax.set_zlim3d(0., 1.)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
