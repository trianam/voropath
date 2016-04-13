#!/usr/bin/python

import numpy as np
import voronizator
import plotter

voronoi = voronizator.Voronizator()

Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])

voronoi.setCustomSites(np.array([[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],[0.,1.,1.],[1.,0.,1.],[0.,0.,1.],[0.5,0.5,0.5]]))
voronoi.makeVoroGraph()
#voronoi.calculateShortestPath(Vs, Ve, 'near')
voronoi.calculateShortestPath(Vs, Ve, 'all')

plt = plotter.Plotter()

voronoi.plotSites(plt)
voronoi.plotGraph(plt)
voronoi.plotShortestPath(plt)

plt.draw()
