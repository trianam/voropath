#!/usr/bin/python

import numpy as np
import voronizator
import plotter

voronoi = voronizator.Voronizator()

Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])

voronoi.setRandomSites(10,0)
voronoi.makeVoroGraph()
voronoi.calculateShortestPath(Vs, Ve)

plt = plotter.Plotter()

voronoi.plotSites(plt)
voronoi.plotShortestPath(plt)
voronoi.plotGraph(plt)

plt.draw()

