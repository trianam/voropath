#!/bin/python

import sys
import pickle
import plotter

if len(sys.argv) >= 2:

    print('Load file', flush=True)
    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    voronoi = record['voronoi']

    print('Build renderer, window and interactor', flush=True)
    plotter = plotter.Plotter()
    
    voronoi.plotPolyhedrons(plotter, verbose = True)
    voronoi.plotSites(plotter, verbose = True)
    #voronoi.plotGraph(ax, edges=False, labels=True)
    #voronoi.plotGraph(ax, pathExtremes=True)
    #voronoi.plotGraph(ax)

    print('Render', flush=True)
    plotter.draw()
    
else:
    print('use: {} sceneFile'.format(sys.argv[0]))
