#!/usr/bin/python

import numpy as np
import voronizator
import tetrahedron
import plotter

voronoi = voronizator.Voronizator()

poly1 = tetrahedron.Tetrahedron(a = [0.1,0.1,0.1], b = [0.1,0.9,0.1], c = [0.9,0.1,0.1], d = [0.5,0.5,0.9])
poly2 = tetrahedron.Tetrahedron(a = [-0.5,0.1,0.1], b = [-0.5,0.9,0.1], c = [-0.9,0.5,0.1], d = [-0.5,0.5,1.5])
poly3 = tetrahedron.Tetrahedron(a = [0.5,1.,0.5], b = [0.7,0.5,0.5], c = [0.7,0.7,1.], d = [0.9,0.7,0.5])


Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])

voronoi.addPolyhedron(poly1)
voronoi.addPolyhedron(poly2)
voronoi.addPolyhedron(poly3)

#voronoi.addBoundingBox([-2.,-2.,-2.], [3.,3.,3.])
voronoi.setPolyhedronsSites()
voronoi.makeVoroGraph()
voronoi.calculateShortestPath(Vs, Ve)

plt = plotter.Plotter()

#voronoi.plotSites(plt)
voronoi.plotPolyhedrons(plt)
voronoi.plotGraph(plt)
voronoi.plotShortestPath(plt)

plt.draw()
