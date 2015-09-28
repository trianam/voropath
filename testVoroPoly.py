#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import voronizator
import polyhedron

voronoi = voronizator.Voronizator()

poly1 = polyhedron.Polyhedron(a = [0.1,0.1,0.1], b = [0.1,0.9,0.1], c = [0.9,0.1,0.1], d = [0.5,0.5,0.9])


fig = plt.figure()
ax = fig.gca(projection='3d')

Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])

voronoi.addPolyhedron(poly1)
#voronoi.addBoundingBox([-1.,-1.,-1.], [2.,2.,2.])
voronoi.setPolyhedronsSites()
voronoi.makeVoroGraph()
voronoi.calculateShortestPath(Vs, Ve)

#voronoi.plotSites(ax)
voronoi.plotPolyhedrons(ax)
voronoi.plotShortestPath(ax)
#voronoi.plotGraph(ax, edges=False, labels=False)

ax.set_xlim3d(-1., 2.)
ax.set_ylim3d(-1., 2.)
ax.set_zlim3d(-1., 2.)
            
plt.show()
