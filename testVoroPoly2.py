#!/usr/bin/python

import numpy as np
import voronizator
import tetrahedron
import plotter

voronoi = voronizator.Voronizator()

maxEmptyArea = 0.1
poly1 = tetrahedron.Tetrahedron(a = [0.1,0.1,0.1], b = [0.1,0.9,0.1], c = [0.9,0.1,0.1], d = [0.5,0.5,0.9], maxEmptyArea=maxEmptyArea)
poly2 = tetrahedron.Tetrahedron(a = [-0.5,0.1,0.1], b = [-0.5,0.9,0.1], c = [-0.9,0.5,0.1], d = [-0.5,0.5,1.5], maxEmptyArea=maxEmptyArea)
poly3 = tetrahedron.Tetrahedron(a = [0.5,1.,0.5], b = [0.7,0.5,0.5], c = [0.7,0.7,1.], d = [0.9,0.7,0.5], maxEmptyArea=maxEmptyArea)
poly4 = tetrahedron.Tetrahedron(a = [0.,0.5,0.], b = [0.,1.,0.], c = [0.,0.7,1.], d = [-0.4,0.7,0.5], maxEmptyArea=maxEmptyArea)


Vs = np.array([0.,0.,0.])
#Ve = np.array([1.,1.,1.])
Ve = np.array([-0.5,2.,2.])

voronoi.addPolyhedron(poly1)
voronoi.addPolyhedron(poly2)
voronoi.addPolyhedron(poly3)
voronoi.addPolyhedron(poly4)

#voronoi.addBoundingBox([-2.,-2.,-2.], [3.,3.,3.], maxEmptyArea=10.)
voronoi.setPolyhedronsSites()
voronoi.makeVoroGraph()
voronoi.calculateShortestPath(Vs, Ve, 'near')
#voronoi.calculateShortestPath(Vs, Ve, 'all')

plt = plotter.Plotter()

#voronoi.plotSites(plt)
voronoi.plotPolyhedrons(plt)
#voronoi.plotGraph(plt, edges=False, labels=True)
#voronoi.plotGraph(plt, pathExtremes=True)
voronoi.plotGraph(plt)
voronoi.plotShortestPath(plt)

plt.draw()

