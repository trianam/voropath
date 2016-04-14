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
    plt = plotter.Plotter()
    
    voronoi.plotPolyhedrons(plt, verbose = True)
    voronoi.plotSites(plt, verbose = True)
    voronoi.plotGraph(plt, verbose = True)

    print('Render', flush=True)
    plt.draw()
    
else:
    print('use: {} sceneFile'.format(sys.argv[0]))
